#include <RcppArmadillo.h>
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <vector>
#include "hicpotts_pair_weight.h"
using namespace Rcpp;


// -----------------------------------------------------------------------------
// Symmetric crossing-barrier penalty between components 2 and 3, used in place
// of a hard ordering constraint (a1 <= a2 <= a3) on the component intercepts.
//
// A hard ordering constraint forces a single labelling of components 2 and 3
// regardless of the data, because it deterministically resolves every tie in
// the same direction. Removing such a constraint without any replacement
// reintroduces label-switching: the chain can flicker rapidly between the two
// orderings with no stable dwell time in either, which post-hoc relabelling
// cannot cleanly resolve.
//
// The barrier below is a smooth, symmetric penalty centred at the crossing
// point (gap = a3 - a2 = 0):
//   penalty(gap) = -kappa * exp( -gap^2 / (2 * w^2) )
// It is a downward Gaussian bump that is near zero once |gap| is a few
// multiples of w away from zero in either direction, and reaches its maximum
// cost (-kappa) only at gap = 0. Because it saturates, it does not keep
// rewarding ever more extreme separation between the two components, so it
// cannot drag either component's estimate to an unrealistic value. It
// discourages lingering near a tie (the mechanism behind label-switching
// flicker) without forbidding crossing and without permanently penalising
// whichever component is currently smaller.
// -----------------------------------------------------------------------------
static inline double crossing_barrier(double gap, double kappa, double w) {
  return -kappa * std::exp(-(gap * gap) / (2.0 * w * w));
}

// Stable log(sigmoid(x)) helper used by the noise-relationship prior below.
static inline double log_sigmoid(double x) {
  if (x >= 0.0) return -std::log1p(std::exp(-x));
  return x - std::log1p(std::exp(x));
}

// Score an elevated component as the false-signal noise component. Slopes are
// compared with component 1 on the standardised log1p-covariate scale, so a
// covariate's numerical units cannot determine the amount of shrinkage. The
// intercept term is a soft, saturating preference for an elevated baseline;
// it never creates a hard truncation.
static double noise_branch_score(const NumericVector &beta1,
                                 const NumericVector &candidate,
                                 const arma::vec &covariate_sds,
                                 double link_sd,
                                 double order_strength,
                                 double order_width) {
  double score = 0.0;
  for (int j = 0; j < 4; j++) {
    const double difference =
      (candidate[j + 1] - beta1[j + 1]) * covariate_sds[j] / link_sd;
    score -= 0.5 * difference * difference;
  }
  score += order_strength *
    log_sigmoid((candidate[0] - beta1[0]) / order_width);
  return score;
}

// Component 3 is the labelled false-signal noise component and is therefore
// coupled directly to component 1. Component 2 is intentionally absent from
// this term: its slopes and intercept remain unrestricted true-signal
// parameters. This asymmetry removes the component-2/3 posterior permutation
// mode instead of attempting to repair it after sampling.
static double noise_relationship_logprior(const NumericVector &beta1,
                                          const NumericVector &beta2,
                                          const NumericVector &beta3,
                                          const arma::vec &covariate_sds,
                                          double link_sd,
                                          double order_strength,
                                          double order_width) {
  (void) beta2;
  return noise_branch_score(beta1, beta3, covariate_sds, link_sd,
                            order_strength, order_width);
}

// Biological orientation diagnostic. Positive values mean that the labelled
// component 3 is more compatible with the component-1-like elevated-noise
// definition than component 2. It is a diagnostic, not an extra target term:
// the actual target contribution remains noise_relationship_logprior().
static double noise_branch_gap(const NumericVector &beta1,
                               const NumericVector &beta2,
                               const NumericVector &beta3,
                               const arma::vec &covariate_sds,
                               double link_sd,
                               double order_strength,
                               double order_width) {
  return noise_branch_score(beta1, beta3, covariate_sds, link_sd,
                            order_strength, order_width) -
         noise_branch_score(beta1, beta2, covariate_sds, link_sd,
                            order_strength, order_width);
}

// Build one set of empirical-Bayes Normal hyperparameters from SOFT component
// membership weights.  The centre is a ridge-stabilised weighted least-squares
// fit of log1p(y) to the same five-column log1p design used by the count model.
// Using conditional membership probabilities rather than hard z labels avoids
// the self-reinforcing z -> prior -> beta -> z feedback that can empty a weak
// component during warm-up.  The sampler updates these hyperparameters only
// before the retained-sample boundary and then freezes them, so retained draws
// target one fixed posterior.
static List soft_eb_priors_impl(
    const NumericMatrix &y,
    const NumericMatrix &cov1, const NumericMatrix &cov2,
    const NumericMatrix &cov3, const NumericMatrix &cov4,
    const NumericMatrix &weight1, const NumericMatrix &weight2,
    const NumericMatrix &weight3,
    bool covariates_are_logged = false) {
  const int nr = y.nrow(), nc = y.ncol();
  if (weight1.nrow() != nr || weight1.ncol() != nc ||
      weight2.nrow() != nr || weight2.ncol() != nc ||
      weight3.nrow() != nr || weight3.ncol() != nc ||
      cov1.nrow() != nr || cov1.ncol() != nc ||
      cov2.nrow() != nr || cov2.ncol() != nc ||
      cov3.nrow() != nr || cov3.ncol() != nc ||
      cov4.nrow() != nr || cov4.ncol() != nc)
    stop("Empirical-Bayes inputs must have identical matrix dimensions.");

  std::vector<arma::mat> xtx(4, arma::mat(5, 5, arma::fill::zeros));
  std::vector<arma::vec> xty(4, arma::vec(5, arma::fill::zeros));
  std::vector<double> yty(4, 0.0);
  std::vector<double> n(4, 0.0);
  arma::vec x_sum(4, arma::fill::zeros);
  arma::vec x_sum2(4, arma::fill::zeros);

  for (int i = 0; i < nr; ++i) {
    for (int j = 0; j < nc; ++j) {
      const double yy = y(i, j);
      const double raw[4] = {cov1(i, j), cov2(i, j), cov3(i, j), cov4(i, j)};
      if (!(yy >= 0.0) || !std::isfinite(yy))
        stop("Empirical-Bayes counts must be finite and non-negative.");
      double weights[3] = {weight1(i, j), weight2(i, j), weight3(i, j)};
      double weight_sum = 0.0;
      for (int k = 0; k < 3; ++k) {
        if (!std::isfinite(weights[k]) || weights[k] < 0.0)
          stop("Empirical-Bayes membership weights must be finite and non-negative.");
        weight_sum += weights[k];
      }
      if (!(weight_sum > 0.0) || !std::isfinite(weight_sum))
        stop("Empirical-Bayes membership weights must have a positive row sum.");
      for (int k = 0; k < 3; ++k) weights[k] /= weight_sum;
      double design_row[5];
      design_row[0] = 1.0;
      for (int p = 0; p < 4; ++p) {
        if (!std::isfinite(raw[p]) ||
            (!covariates_are_logged && !(raw[p] > -1.0)))
          stop("Empirical-Bayes covariates are invalid.");
        design_row[p + 1] = covariates_are_logged ? raw[p] : std::log1p(raw[p]);
        x_sum[p] += design_row[p + 1];
        x_sum2[p] += design_row[p + 1] * design_row[p + 1];
      }
      const double response = std::log1p(yy);
      for (int group = 0; group < 3; ++group) {
        const double weight = weights[group];
        for (int a = 0; a < 5; ++a) {
          const double xy = design_row[a] * response;
          xty[group][a] += weight * xy;
          for (int b = 0; b < 5; ++b) {
            const double xx = design_row[a] * design_row[b];
            xtx[group](a, b) += weight * xx;
          }
        }
        yty[group] += weight * response * response;
        n[group] += weight;
      }
      for (int a = 0; a < 5; ++a) {
        const double xy = design_row[a] * response;
        xty[3][a] += xy;
        for (int b = 0; b < 5; ++b)
          xtx[3](a, b) += design_row[a] * design_row[b];
      }
      yty[3] += response * response;
      n[3] += 1.0;
    }
  }

  const double total_n = std::max(1.0, n[3]);
  arma::vec x_sd(4, arma::fill::ones);
  for (int p = 0; p < 4; ++p) {
    const double variance = n[3] > 1 ?
      (x_sum2[p] - x_sum[p] * x_sum[p] / total_n) / (total_n - 1.0) : 0.0;
    x_sd[p] = std::sqrt(std::max(variance, 1e-8));
  }

  auto one_component = [&](int component) {
    const bool fallback = n[component] < 8;
    const int source = fallback ? 3 : component;
    const double source_n = std::max(1.0, n[source]);
    arma::mat gram = xtx[source];
    const double ridge = std::max(1e-8,
      1e-6 * arma::trace(gram) / (5.0 * source_n));
    arma::mat regularised = gram + ridge * source_n * arma::eye(5, 5);
    arma::vec centre;
    bool solved = arma::solve(centre, regularised, xty[source],
                              arma::solve_opts::likely_sympd);
    if (!solved || !centre.is_finite())
      centre = arma::pinv(regularised) * xty[source];
    if (!centre.is_finite()) centre.zeros(5);

    double rss = yty[source] - 2.0 * arma::dot(centre, xty[source]) +
      arma::as_scalar(centre.t() * gram * centre);
    const double residual_variance = std::max(
      rss / std::max(1.0, source_n - 5.0), 0.05 * 0.05);
    arma::mat information = gram / source_n + ridge * arma::eye(5, 5);
    arma::mat information_inverse;
    bool inverted = arma::inv_sympd(information_inverse, information);
    if (!inverted || !information_inverse.is_finite())
      information_inverse = arma::pinv(information);
    arma::vec prior_sd = arma::sqrt(
      arma::clamp(residual_variance * information_inverse.diag(), 1e-12,
                  arma::datum::inf));

    prior_sd[0] = std::min(5.0, std::max(0.5, prior_sd[0]));
    for (int p = 0; p < 4; ++p) {
      const double floor_sd = 0.25 / x_sd[p];
      const double ceiling_sd = 2.5 / x_sd[p];
      prior_sd[p + 1] = std::min(ceiling_sd,
        std::max(floor_sd, prior_sd[p + 1]));
    }
    if (fallback) prior_sd *= 1.5;

    return List::create(
      Named("meany") = centre[0], Named("meanx1") = centre[1],
      Named("meanx2") = centre[2], Named("meanx3") = centre[3],
      Named("meanx4") = centre[4], Named("sdy") = prior_sd[0],
      Named("sdx1") = prior_sd[1], Named("sdx2") = prior_sd[2],
      Named("sdx3") = prior_sd[3], Named("sdx4") = prior_sd[4],
      Named("n_cells") = n[component], Named("effective_cells") = n[component],
      Named("allocation_weights") = "soft membership probabilities",
      Named("global_fallback") = fallback);
  };

  return List::create(
    Named("component1") = one_component(0),
    Named("component2") = one_component(1),
    Named("component3") = one_component(2));
}

// Backwards-compatible hard-allocation helper for direct R calls and tests.
// The production sampler uses soft_eb_priors_impl() during warm-up.
static List iterated_eb_priors_impl(
    const NumericMatrix &y,
    const NumericMatrix &cov1, const NumericMatrix &cov2,
    const NumericMatrix &cov3, const NumericMatrix &cov4,
    const NumericMatrix &z) {
  const int nr = y.nrow(), nc = y.ncol();
  if (z.nrow() != nr || z.ncol() != nc)
    stop("Empirical-Bayes y and z values must have identical dimensions.");
  NumericMatrix weight1(nr, nc), weight2(nr, nc), weight3(nr, nc);
  for (int i = 0; i < nr; ++i) {
    for (int j = 0; j < nc; ++j) {
      const double value = z(i, j);
      if (!std::isfinite(value) || value < 1.0 || value > 3.0 ||
          std::fabs(value - std::round(value)) > 1e-8)
        stop("Empirical-Bayes z values must be integer labels 1, 2 or 3.");
      const int label = (int) std::round(value);
      weight1(i, j) = label == 1 ? 1.0 : 0.0;
      weight2(i, j) = label == 2 ? 1.0 : 0.0;
      weight3(i, j) = label == 3 ? 1.0 : 0.0;
    }
  }
  return soft_eb_priors_impl(y, cov1, cov2, cov3, cov4,
                             weight1, weight2, weight3);
}

// Internal deterministic test and public-R helper hook.
// [[Rcpp::export(name = ".hicpotts_iterated_eb_priors_cpp")]]
List hicpotts_iterated_eb_priors_cpp(
    NumericMatrix y, NumericMatrix cov1, NumericMatrix cov2,
    NumericMatrix cov3, NumericMatrix cov4, NumericMatrix z) {
  return iterated_eb_priors_impl(y, cov1, cov2, cov3, cov4, z);
}

// Internal hook used to verify soft-allocation EB and its permutation
// invariance without duplicating the weighted regression in R.
// [[Rcpp::export(name = ".hicpotts_soft_eb_priors_cpp")]]
List hicpotts_soft_eb_priors_cpp(
    NumericMatrix y, NumericMatrix cov1, NumericMatrix cov2,
    NumericMatrix cov3, NumericMatrix cov4,
    NumericMatrix weight1, NumericMatrix weight2, NumericMatrix weight3) {
  return soft_eb_priors_impl(y, cov1, cov2, cov3, cov4,
                             weight1, weight2, weight3);
}

static void swap_signal_labels(NumericMatrix &z) {
  for (int i = 0; i < z.nrow(); ++i) {
    for (int j = 0; j < z.ncol(); ++j) {
      const int label = (int) z(i, j);
      if (label == 2) z(i, j) = 3;
      else if (label == 3) z(i, j) = 2;
    }
  }
}

static inline double logspace_add(double x, double y) {
  if (!std::isfinite(x)) return y;
  if (!std::isfinite(y)) return x;
  const double m = std::max(x, y);
  return m + std::log(std::exp(x - m) + std::exp(y - m));
}

// -----------------------------------------------------------------------------
// Cached native emission/posterior engine.
//
// The public R helpers remain unchanged, but the production sampler evaluates
// its hottest likelihood paths entirely in C++. The four log1p covariate
// matrices are created once per fit and reused in every z, beta, size, theta
// and empirical-Bayes update. This removes repeated transformations, logical
// subsetting and R callbacks without changing the target density.
// -----------------------------------------------------------------------------
static NumericMatrix log1p_matrix(const NumericMatrix &x) {
  NumericMatrix out(x.nrow(), x.ncol());
  for (int j = 0; j < x.ncol(); ++j) {
    for (int i = 0; i < x.nrow(); ++i) {
      const double value = x(i, j);
      if (!std::isfinite(value) || value <= -1.0)
        stop("Covariates must be finite and greater than -1.");
      out(i, j) = std::log1p(value);
    }
  }
  return out;
}

struct HicpottsNormalPrior {
  std::array<double, 5> mean;
  std::array<double, 5> sd;
};

static HicpottsNormalPrior normal_prior_from_list(const RObject &object,
                                                  int component) {
  if (object.isNULL()) stop("Regression-prior hyperparameters are missing.");
  List priors(object);
  const std::string key = "component" + std::to_string(component);
  if (!priors.containsElementNamed(key.c_str()))
    stop("Regression priors are missing %s.", key.c_str());
  List p = priors[key];
  const char *mean_names[5] = {"meany", "meanx1", "meanx2", "meanx3", "meanx4"};
  const char *sd_names[5] = {"sdy", "sdx1", "sdx2", "sdx3", "sdx4"};
  HicpottsNormalPrior out;
  for (int k = 0; k < 5; ++k) {
    if (!p.containsElementNamed(mean_names[k]) ||
        !p.containsElementNamed(sd_names[k]))
      stop("Regression prior %s is incomplete.", key.c_str());
    out.mean[k] = as<double>(p[mean_names[k]]);
    out.sd[k] = std::max(1e-6, as<double>(p[sd_names[k]]));
    if (!std::isfinite(out.mean[k]) || !std::isfinite(out.sd[k]) ||
        !(out.sd[k] > 0.0))
      stop("Regression prior %s contains invalid values.", key.c_str());
  }
  return out;
}

static inline double linear_predictor_cached(
    const NumericVector &beta, int i, int j,
    const NumericMatrix &logcov1, const NumericMatrix &logcov2,
    const NumericMatrix &logcov3, const NumericMatrix &logcov4) {
  double eta = beta[0] + beta[1] * logcov1(i, j) +
    beta[2] * logcov2(i, j) + beta[3] * logcov3(i, j) +
    beta[4] * logcov4(i, j);
  if (eta > 30.0) eta = 30.0;
  else if (eta < -30.0) eta = -30.0;
  return eta;
}

static inline double emission_loglik_cached(
    int component, double count, int i, int j, const NumericVector &beta,
    const NumericMatrix &logcov1, const NumericMatrix &logcov2,
    const NumericMatrix &logcov3, const NumericMatrix &logcov4,
    double size_value, double theta, const std::string &dist,
    bool allocation_floor = false) {
  const double mu = std::exp(linear_predictor_cached(
    beta, i, j, logcov1, logcov2, logcov3, logcov4));
  const bool negative_binomial = dist == "NB" || dist == "ZINB";
  double base = negative_binomial ?
    R::dnbinom_mu(count, size_value, mu, true) :
    R::dpois(count, mu, true);
  double out = base;
  if (component == 1 && (dist == "ZIP" || dist == "ZINB")) {
    if (count == 0.0) {
      out = logspace_add(std::log(theta), std::log1p(-theta) + base);
    } else {
      out = std::log1p(-theta) + base;
    }
  }
  // The historical pz_123 path floors natural-scale allocation densities at
  // 1e-300. Preserve that convention exactly for z updates, while parameter
  // posteriors retain the full log density without flooring.
  if (allocation_floor && (!std::isfinite(out) || out < std::log(1e-300)))
    out = std::log(1e-300);
  return out;
}

static inline double normal_logprior_cached(
    const NumericVector &beta, const HicpottsNormalPrior &prior) {
  double out = 0.0;
  for (int k = 0; k < 5; ++k)
    out += R::dnorm(beta[k], prior.mean[k], prior.sd[k], true);
  return out;
}

static inline double size_logprior_cached(double size_value, int component) {
  const double shape = component == 1 ? 3.0 : 2.0;
  const double rate = component == 1 ? 1.0 : 0.2;
  return R::dgamma(size_value, shape, 1.0 / rate, true);
}

static double component_logposterior_cached(
    const NumericVector &beta, const NumericMatrix &z,
    const NumericMatrix &y, int component, double theta,
    double size_value, const std::string &dist,
    const NumericMatrix &logcov1, const NumericMatrix &logcov2,
    const NumericMatrix &logcov3, const NumericMatrix &logcov4,
    const HicpottsNormalPrior &prior) {
  double out = normal_logprior_cached(beta, prior);
  if (dist == "NB" || dist == "ZINB")
    out += size_logprior_cached(size_value, component);
  for (int j = 0; j < y.ncol(); ++j) {
    for (int i = 0; i < y.nrow(); ++i) {
      if ((int)z(i, j) != component) continue;
      out += emission_loglik_cached(
        component, y(i, j), i, j, beta, logcov1, logcov2,
        logcov3, logcov4, size_value, theta, dist, false);
    }
  }
  return out;
}

// Internal equivalence hook used by tests; inputs use the public raw-covariate
// convention and are transformed once before calling the cached engine.
// [[Rcpp::export(name = ".hicpotts_cached_posterior_cpp")]]
double hicpotts_cached_posterior_cpp(
    NumericVector beta, NumericMatrix z, NumericMatrix y, int component,
    double theta, double size_value, std::string dist,
    NumericMatrix cov1, NumericMatrix cov2, NumericMatrix cov3,
    NumericMatrix cov4, List priors) {
  if (beta.size() != 5 || component < 1 || component > 3)
    stop("Expected five coefficients and component 1, 2 or 3.");
  NumericMatrix l1 = log1p_matrix(cov1), l2 = log1p_matrix(cov2);
  NumericMatrix l3 = log1p_matrix(cov3), l4 = log1p_matrix(cov4);
  return component_logposterior_cached(
    beta, z, y, component, theta, size_value, dist, l1, l2, l3, l4,
    normal_prior_from_list(priors, component));
}

// Density on the gamma scale induced by a Gaussian random walk on logit(gamma).
static double gamma_logit_rw_logdensity(double to, double from, double step) {
  if (!(to > 0.0 && to < 1.0) || !(from > 0.0 && from < 1.0) ||
      !(step > 0.0) || !std::isfinite(step)) return R_NegInf;
  const double eta_to = std::log(to / (1.0 - to));
  const double eta_from = std::log(from / (1.0 - from));
  return R::dnorm(eta_to, eta_from, step, true) -
         std::log(to) - std::log1p(-to);
}

// Full proposal density for the gamma mixture. gamma_large_probability is
// conditional on selecting the random-walk part; gamma_independence_probability
// selects a direct draw from the configured Beta prior. Evaluating the full
// mixture in both directions makes the transition valid regardless of which
// component generated the realised proposal.
static double gamma_mixture_logdensity(
    double to, double from, double local_step,
    double gamma_large_probability, double gamma_large_multiplier,
    double gamma_independence_probability,
    double gamma_prior_shape1, double gamma_prior_shape2) {
  const double rw_weight = 1.0 - gamma_independence_probability;
  const double small_weight = rw_weight * (1.0 - gamma_large_probability);
  const double large_weight = rw_weight * gamma_large_probability;
  double out = R_NegInf;
  if (small_weight > 0.0)
    out = std::log(small_weight) +
      gamma_logit_rw_logdensity(to, from, local_step);
  if (large_weight > 0.0)
    out = logspace_add(out, std::log(large_weight) +
      gamma_logit_rw_logdensity(to, from,
                                local_step * gamma_large_multiplier));
  if (gamma_independence_probability > 0.0)
    out = logspace_add(out, std::log(gamma_independence_probability) +
      R::dbeta(to, gamma_prior_shape1, gamma_prior_shape2, true));
  return out;
}

static double draw_gamma_mixture(
    double current, double local_step,
    double gamma_large_probability, double gamma_large_multiplier,
    double gamma_independence_probability,
    double gamma_prior_shape1, double gamma_prior_shape2,
    int &proposal_type) {
  if (R::runif(0.0, 1.0) < gamma_independence_probability) {
    proposal_type = 2;
    return R::rbeta(gamma_prior_shape1, gamma_prior_shape2);
  }
  const bool wide = R::runif(0.0, 1.0) < gamma_large_probability;
  proposal_type = wide ? 1 : 0;
  const double step = local_step * (wide ? gamma_large_multiplier : 1.0);
  const double eta = std::log(current / (1.0 - current)) +
                     R::rnorm(0.0, step);
  return 1.0 / (1.0 + std::exp(-eta));
}

// Internal test hook: this is intentionally not exported from NAMESPACE.
// [[Rcpp::export(name = ".hicpotts_noise_relationship_logprior_cpp")]]
double hicpotts_noise_relationship_logprior_cpp(NumericVector beta1,
                                                NumericVector beta2,
                                                NumericVector beta3,
                                                NumericVector covariate_sds,
                                                double link_sd = 0.5,
                                                double order_strength = 10.0,
                                                double order_width = 0.5) {
  if (beta1.size() != 5 || beta2.size() != 5 || beta3.size() != 5)
    stop("Each beta vector must contain an intercept and four slopes.");
  if (covariate_sds.size() != 4)
    stop("covariate_sds must contain four values.");
  if (!(link_sd > 0.0) || !std::isfinite(link_sd) ||
      !(order_width > 0.0) || !std::isfinite(order_width) ||
      !(order_strength >= 0.0) || !std::isfinite(order_strength))
    stop("Invalid noise-relationship prior settings.");
  arma::vec sds(4);
  for (int j = 0; j < 4; j++) {
    if (!(covariate_sds[j] > 0.0) || !std::isfinite(covariate_sds[j]))
      stop("covariate_sds must be positive and finite.");
    sds[j] = covariate_sds[j];
  }
  return noise_relationship_logprior(beta1, beta2, beta3, sds, link_sd,
                                     order_strength, order_width);
}

// Internal detailed-balance hook for the whole 2/3 branch transition. The
// returned value is the labelled biological-prior contribution to
// log pi(T(x)) - log pi(x). Calling it again with beta2/beta3 exchanged must
// return its exact negative.
// [[Rcpp::export(name = ".hicpotts_branch_prior_logratio_cpp")]]
double hicpotts_branch_prior_logratio_cpp(
    NumericVector beta1, NumericVector beta2, NumericVector beta3,
    NumericVector covariate_sds, double link_sd = 0.5,
    double order_strength = 10.0, double order_width = 0.5) {
  if (beta1.size() != 5 || beta2.size() != 5 || beta3.size() != 5 ||
      covariate_sds.size() != 4)
    stop("Expected three length-five beta vectors and four covariate SDs.");
  arma::vec sds(4);
  for (int j = 0; j < 4; ++j) {
    if (!(covariate_sds[j] > 0.0) || !std::isfinite(covariate_sds[j]))
      stop("covariate_sds must be positive and finite.");
    sds[j] = covariate_sds[j];
  }
  const double current = noise_relationship_logprior(
    beta1, beta2, beta3, sds, link_sd, order_strength, order_width);
  const double swapped = noise_relationship_logprior(
    beta1, beta3, beta2, sds, link_sd, order_strength, order_width);
  return swapped - current;
}

// Internal test hook for the non-symmetric part of the gamma proposal ratio.
// Swapping current/proposed must negate the result.
// [[Rcpp::export(name = ".hicpotts_gamma_proposal_logratio_cpp")]]
double hicpotts_gamma_proposal_logratio_cpp(
    double current, double proposed, double local_step = 0.5,
    double large_probability = 0.10, double large_multiplier = 4.0,
    double independence_probability = 0.05,
    double prior_shape1 = 1.0, double prior_shape2 = 1.0) {
  return gamma_mixture_logdensity(
    current, proposed, local_step, large_probability, large_multiplier,
    independence_probability, prior_shape1, prior_shape2) -
    gamma_mixture_logdensity(
      proposed, current, local_step, large_probability, large_multiplier,
      independence_probability, prior_shape1, prior_shape2);
}

// Compute centring and scaling from the cached transformed covariates. This is
// proposal preconditioning only; returned coefficients remain on the original
// log1p-covariate scale.
static void cached_mean_sd(const NumericMatrix &x, double &mean, double &sd) {
  const int n = x.nrow() * x.ncol();
  if (n < 1) stop("Covariate matrices must not be empty.");
  double sum = 0.0;
  for (int j = 0; j < x.ncol(); ++j)
    for (int i = 0; i < x.nrow(); ++i) {
      if (!std::isfinite(x(i, j))) stop("Cached covariates must be finite.");
      sum += x(i, j);
    }
  mean = sum / (double)n;
  double ss = 0.0;
  for (int j = 0; j < x.ncol(); ++j)
    for (int i = 0; i < x.nrow(); ++i) {
      const double centred = x(i, j) - mean;
      ss += centred * centred;
    }
  sd = n > 1 ? std::sqrt(ss / (double)(n - 1)) : 1.0;
  if (!std::isfinite(sd) || sd < 1e-8) sd = 1.0;
}

static arma::vec beta_to_standardized(const arma::vec &beta,
                                      const arma::vec &means,
                                      const arma::vec &sds) {
  arma::vec q(5, arma::fill::zeros);
  q[0] = beta[0];
  for (int j = 0; j < 4; j++) {
    q[0] += beta[j + 1] * means[j];
    q[j + 1] = beta[j + 1] * sds[j];
  }
  return q;
}

static arma::vec standardized_to_beta(const arma::vec &q,
                                      const arma::vec &means,
                                      const arma::vec &sds) {
  arma::vec beta(5, arma::fill::zeros);
  beta[0] = q[0];
  for (int j = 0; j < 4; j++) {
    beta[j + 1] = q[j + 1] / sds[j];
    beta[0] -= beta[j + 1] * means[j];
  }
  return beta;
}

// Orthogonal proposal coordinates. If Xs is the centred/standardised design
// and Xs = Q R, use u = (R / sqrt(n)) q. The constant linear transformation
// cancels from every symmetric MH proposal ratio; posterior evaluation and
// priors remain on the original beta scale.
static arma::vec beta_to_qr(const arma::vec &beta,
                            const arma::vec &means,
                            const arma::vec &sds,
                            const arma::mat &R_scaled) {
  return R_scaled * beta_to_standardized(beta, means, sds);
}

static arma::vec qr_to_beta(const arma::vec &u,
                            const arma::vec &means,
                            const arma::vec &sds,
                            const arma::mat &R_scaled) {
  arma::vec q;
  bool solved = arma::solve(q, arma::trimatu(R_scaled), u,
                            arma::solve_opts::fast);
  if (!solved || q.n_elem != 5 || !q.is_finite())
    stop("QR proposal back-transformation failed.");
  return standardized_to_beta(q, means, sds);
}

// Batch-means Monte Carlo standard error divided by the ordinary sample SD.
// This is an internal stopping diagnostic only and never enters an acceptance
// probability. Returning Inf for short/invalid series prevents early stop.
static double relative_mcse(const std::vector<double> &values) {
  const int n = (int)values.size();
  if (n < 100) return R_PosInf;
  double mean = 0.0;
  for (int i = 0; i < n; i++) mean += values[i];
  mean /= (double)n;
  double ss = 0.0;
  for (int i = 0; i < n; i++) {
    double d = values[i] - mean;
    ss += d * d;
  }
  double variance = ss / (double)(n - 1);
  if (!std::isfinite(variance)) return R_PosInf;
  // A near-constant series is a frozen chain, not a perfectly precise one.
  // Returning 0 here declared maximal precision and allowed the sampler to
  // stop early on a stuck parameter. Inf prevents early stopping instead.
  if (variance <= 1e-16) return R_PosInf;
  int batch_size = std::max(10, (int)std::floor(std::sqrt((double)n)));
  int batches = n / batch_size;
  if (batches < 5) return R_PosInf;
  std::vector<double> batch_means(batches, 0.0);
  for (int b = 0; b < batches; b++) {
    for (int j = 0; j < batch_size; j++)
      batch_means[b] += values[b * batch_size + j];
    batch_means[b] /= (double)batch_size;
  }
  double bm_mean = 0.0;
  for (int b = 0; b < batches; b++) bm_mean += batch_means[b];
  bm_mean /= (double)batches;
  double bm_ss = 0.0;
  for (int b = 0; b < batches; b++) {
    double d = batch_means[b] - bm_mean;
    bm_ss += d * d;
  }
  double spectral_variance = batch_size * bm_ss / (double)(batches - 1);
  double mcse = std::sqrt(std::max(0.0, spectral_variance) / (double)n);
  return mcse / std::sqrt(variance);
}

static inline bool is_signal_component(int label) {
  return label == 2 || label == 3;
}

static inline double signal_component_loglik_cell(int comp, int i, int j,
                                                  const NumericMatrix &y,
                                                  const NumericVector &params,
                                                  const NumericMatrix &logcov1,
                                                  const NumericMatrix &logcov2,
                                                  const NumericMatrix &logcov3,
                                                  const NumericMatrix &logcov4,
                                                  double size_value,
                                                  const std::string &dist) {
  (void) comp;  // signal components share the same emission family.
  double eta = params[0] +
               params[1] * logcov1(i, j) +
               params[2] * logcov2(i, j) +
               params[3] * logcov3(i, j) +
               params[4] * logcov4(i, j);
  if (eta > 30.0) eta = 30.0; else if (eta < -30.0) eta = -30.0;
  double mu = std::exp(eta);
  double yy = y(i, j);

  if (dist == "Poisson" || dist == "ZIP") {
    return R::dpois(yy, mu, true);
  } else if (dist == "NB" || dist == "ZINB") {
    return R::dnbinom_mu(yy, size_value, mu, true);
  }
  stop("Invalid distribution specified for signal-cluster block move.");
  return R_NegInf;
}

// Non-local cluster move that swaps a spatially-contiguous cluster of cells
// between the two signal components (2 and 3), in the style of a Wolff/
// Swendsen-Wang cluster update for the Potts model. It is a z-only move,
// layered on top of the crossing-barrier penalty above; it never touches
// component 1, and zero counts allocated to a signal component remain
// eligible for the cluster.
static bool signal_cluster_swap_23(NumericMatrix &z,
                                   const NumericMatrix &y,
                                   const NumericMatrix &logcov1,
                                   const NumericMatrix &logcov2,
                                   const NumericMatrix &logcov3,
                                   const NumericMatrix &logcov4,
                                   List chains,
                                   const NumericMatrix &size_chain,
                                   int N, int iter, double gamma_iter,
                                   const std::string &dist) {
  std::vector<int> eligible;
  eligible.reserve(N * N);
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      int label = (int) z(i, j);
      if (is_signal_component(label)) {
        eligible.push_back(i * N + j);
      }
    }
  }
  if (eligible.empty()) return false;

  int seed_pos = (int) std::floor(R::runif(0.0, (double) eligible.size()));
  if (seed_pos >= (int) eligible.size()) seed_pos = (int) eligible.size() - 1;
  int seed = eligible[seed_pos];
  int old_label = (int) z(seed / N, seed % N);
  if (!is_signal_component(old_label)) return false;
  int new_label = (old_label == 2) ? 3 : 2;

  double bond_prob = 1.0 - std::exp(-gamma_iter);
  if (!std::isfinite(bond_prob) || bond_prob < 0.0) bond_prob = 0.0;
  if (bond_prob > 1.0) bond_prob = 1.0;

  std::vector<unsigned char> in_cluster(N * N, 0);
  std::vector<int> cluster;
  std::vector<int> stack;
  cluster.reserve(eligible.size());
  stack.push_back(seed);
  in_cluster[seed] = 1;
  cluster.push_back(seed);

  const int di[4] = {-1, 1, 0, 0};
  const int dj[4] = {0, 0, -1, 1};
  while (!stack.empty()) {
    int idx = stack.back();
    stack.pop_back();
    int i = idx / N, j = idx % N;
    for (int d = 0; d < 4; d++) {
      int ni = i + di[d], nj = j + dj[d];
      if (ni < 0 || ni >= N || nj < 0 || nj >= N) continue;
      int nidx = ni * N + nj;
      if (in_cluster[nidx]) continue;
      // Fortuin-Kasteleyn bonds depend only on the Potts labels.  Excluding
      // zero-count neighbours here made the forward cluster law depend on the
      // observations while the MH boundary correction below still assumed
      // the same bond probability on every equal-label edge.  That mismatch
      // violates detailed balance.  A zero count remains a perfectly valid
      // draw from either signal emission; its likelihood contribution belongs
      // in loglik_delta, not in the cluster-construction rule.
      if ((int) z(ni, nj) != old_label) continue;
      if (R::runif(0.0, 1.0) < bond_prob) {
        in_cluster[nidx] = 1;
        cluster.push_back(nidx);
        stack.push_back(nidx);
      }
    }
  }

  int boundary_old = 0, boundary_new = 0;
  for (size_t c = 0; c < cluster.size(); c++) {
    int idx = cluster[c];
    int i = idx / N, j = idx % N;
    for (int d = 0; d < 4; d++) {
      int ni = i + di[d], nj = j + dj[d];
      if (ni < 0 || ni >= N || nj < 0 || nj >= N) continue;
      int nidx = ni * N + nj;
      if (in_cluster[nidx]) continue;
      int outside_label = (int) z(ni, nj);
      if (outside_label == old_label) boundary_old++;
      else if (outside_label == new_label) boundary_new++;
    }
  }

  NumericMatrix chain_old = chains[old_label - 1];
  NumericMatrix chain_new = chains[new_label - 1];
  NumericVector beta_old = chain_old(iter, _);
  NumericVector beta_new = chain_new(iter, _);
  double size_old = size_chain(old_label - 1, iter);
  double size_new = size_chain(new_label - 1, iter);
  if ((dist == "NB" || dist == "ZINB") &&
      (!(size_old > 0.0) || !(size_new > 0.0) ||
       !std::isfinite(size_old) || !std::isfinite(size_new))) {
    return false;
  }

  // Full-matrix weighting, identical to the beta/dispersion and theta paths,
  // so this cluster move targets the same posterior as every other transition.
  // See src/hicpotts_pair_weight.h.
  const bool y_pairs_symmetric = hicpotts_symmetric_pairs(y);

  double loglik_delta = 0.0;
  for (size_t c = 0; c < cluster.size(); c++) {
    int idx = cluster[c];
    int i = idx / N, j = idx % N;
    const double w = hicpotts_pair_weight(i, j, y_pairs_symmetric);
    loglik_delta += w * (signal_component_loglik_cell(new_label, i, j, y, beta_new,
                                                      logcov1, logcov2, logcov3, logcov4,
                                                      size_new, dist) -
                         signal_component_loglik_cell(old_label, i, j, y, beta_old,
                                                      logcov1, logcov2, logcov3, logcov4,
                                                      size_old, dist));
  }
  if (!std::isfinite(loglik_delta)) return false;

  int boundary_delta = boundary_new - boundary_old;
  double log_potts_ratio = gamma_iter * (double) boundary_delta;
  double log_proposal_ratio = 0.0;
  if (bond_prob > 0.0 && bond_prob < 1.0) {
    log_proposal_ratio = std::log1p(-bond_prob) * (double) boundary_delta;
  }
  double log_alpha = loglik_delta + log_potts_ratio + log_proposal_ratio;

  if (std::log(R::runif(0.0, 1.0)) < log_alpha) {
    for (size_t c = 0; c < cluster.size(); c++) {
      int idx = cluster[c];
      z(idx / N, idx % N) = new_label;
    }
    return true;
  }
  return false;
}

// Internal fixed-parameter test hook.  It repeatedly applies the production
// component-2/3 cluster transition to the supplied state and counts the
// resulting binary 2/3 states.  Small-lattice tests can therefore verify
// detailed balance without maintaining a second implementation in R.
// [[Rcpp::export(.hicpotts_cluster_transition_counts)]]
IntegerVector hicpotts_cluster_transition_counts(
    NumericMatrix z, NumericMatrix y,
    NumericMatrix cov1, NumericMatrix cov2,
    NumericMatrix cov3, NumericMatrix cov4,
    NumericMatrix beta, NumericVector sizes,
    double gamma, std::string dist, int repetitions = 10000) {
  int N = z.nrow();
  if (N != z.ncol() || N < 1 || N > 3)
    stop("The cluster diagnostic requires a square lattice with N in 1..3.");
  if (y.nrow() != N || y.ncol() != N ||
      cov1.nrow() != N || cov1.ncol() != N ||
      cov2.nrow() != N || cov2.ncol() != N ||
      cov3.nrow() != N || cov3.ncol() != N ||
      cov4.nrow() != N || cov4.ncol() != N)
    stop("All diagnostic matrices must have the same dimensions as z.");
  if (beta.nrow() != 3 || beta.ncol() != 5 || sizes.size() != 3)
    stop("beta must be 3 by 5 and sizes must have length three.");
  if (repetitions < 1)
    stop("repetitions must be positive.");
  for (int i = 0; i < N; ++i)
    for (int j = 0; j < N; ++j)
      if ((int)z(i, j) != 2 && (int)z(i, j) != 3)
        stop("The cluster diagnostic state may contain only labels 2 and 3.");

  List chains(3);
  for (int comp = 0; comp < 3; ++comp) {
    NumericMatrix chain(1, 5);
    for (int col = 0; col < 5; ++col)
      chain(0, col) = beta(comp, col);
    chains[comp] = chain;
  }
  NumericMatrix size_chain(3, 1);
  for (int comp = 0; comp < 3; ++comp)
    size_chain(comp, 0) = sizes[comp];
  NumericMatrix logcov1 = log1p_matrix(cov1);
  NumericMatrix logcov2 = log1p_matrix(cov2);
  NumericMatrix logcov3 = log1p_matrix(cov3);
  NumericMatrix logcov4 = log1p_matrix(cov4);

  int state_count = 1 << (N * N);
  IntegerVector counts(state_count, 0);
  for (int rep = 0; rep < repetitions; ++rep) {
    NumericMatrix proposal = clone(z);
    signal_cluster_swap_23(
      proposal, y, logcov1, logcov2, logcov3, logcov4, chains, size_chain,
      N, 0, gamma, dist);
    int code = 0;
    for (int i = 0; i < N; ++i)
      for (int j = 0; j < N; ++j)
        if ((int)proposal(i, j) == 3)
          code |= 1 << (i * N + j);
    counts[code]++;
  }
  return counts;
}

// -----------------------------------------------------------------------------
// Helper: label_neighbour_counts(z, i, j, N, n1, n2, n3)
// For cell (i,j), counts how many of its up/down/left/right neighbours
// (border-clipped, matching Neighbours_combined()'s convention) currently
// hold each of labels 1, 2, 3. This is the same quantity Neighbours_combined
// computes for a cell's own actual label, generalised here to all three
// hypothetical labels at once. Used by the theta update and by the gamma
// ABC's Potts field simulator below.
// -----------------------------------------------------------------------------
static inline void label_neighbour_counts(const NumericMatrix &z, int i, int j, int N,
                                          int &n1, int &n2, int &n3) {
  n1 = n2 = n3 = 0;
  auto tally = [&](int ni, int nj) {
    int lbl = (int) z(ni, nj);
    if (lbl == 1) n1++; else if (lbl == 2) n2++; else if (lbl == 3) n3++;
  };
  if (i > 0)     tally(i - 1, j);
  if (i < N - 1) tally(i + 1, j);
  if (j > 0)     tally(i, j - 1);
  if (j < N - 1) tally(i, j + 1);
}

static inline void cell_logweights_cached(
    double out[3], int i, int j, const NumericMatrix &z,
    const NumericMatrix &y, const NumericVector betas[3],
    const double sizes[3], double theta, double gamma,
    const std::string &dist, const NumericMatrix &logcov1,
    const NumericMatrix &logcov2, const NumericMatrix &logcov3,
    const NumericMatrix &logcov4) {
  int n1, n2, n3;
  label_neighbour_counts(z, i, j, z.nrow(), n1, n2, n3);
  const int neighbours[3] = {n1, n2, n3};
  for (int k = 0; k < 3; ++k) {
    out[k] = gamma * neighbours[k] + emission_loglik_cached(
      k + 1, y(i, j), i, j, betas[k], logcov1, logcov2, logcov3,
      logcov4, sizes[k], theta, dist, true);
  }
}

static void checkerboard_z_sweep_cached(
    NumericMatrix &z, const NumericMatrix &y, List chains,
    const NumericMatrix &size_chain, int iter, double theta, double gamma,
    const std::string &dist, double inverse_temperature,
    const NumericMatrix &logcov1, const NumericMatrix &logcov2,
    const NumericMatrix &logcov3, const NumericMatrix &logcov4,
    NumericMatrix &weight1, NumericMatrix &weight2, NumericMatrix &weight3) {
  NumericMatrix c1 = chains[0], c2 = chains[1], c3 = chains[2];
  NumericVector betas[3] = {c1(iter, _), c2(iter, _), c3(iter, _)};
  const double sizes[3] = {
    size_chain(0, iter), size_chain(1, iter), size_chain(2, iter)
  };
  const int N = z.nrow();
  for (int colour = 0; colour <= 1; ++colour) {
    for (int j = 0; j < N; ++j) {
      for (int i = 0; i < N; ++i) {
        if (((i + j) & 1) != colour) continue;
        double logw[3];
        cell_logweights_cached(logw, i, j, z, y, betas, sizes, theta,
                               gamma, dist, logcov1, logcov2, logcov3,
                               logcov4);
        const double m = std::max(logw[0], std::max(logw[1], logw[2]));
        double probabilities[3] = {
          std::exp(inverse_temperature * (logw[0] - m)),
          std::exp(inverse_temperature * (logw[1] - m)),
          std::exp(inverse_temperature * (logw[2] - m))
        };
        const double total = probabilities[0] + probabilities[1] + probabilities[2];
        for (int k = 0; k < 3; ++k) probabilities[k] /= total;
        weight1(i, j) = probabilities[0];
        weight2(i, j) = probabilities[1];
        weight3(i, j) = probabilities[2];
        const double u = R::runif(0.0, 1.0);
        z(i, j) = u < probabilities[0] ? 1.0 :
          (u < probabilities[0] + probabilities[1] ? 2.0 : 3.0);
      }
    }
  }
}

// Internal test hook for the direct checkerboard conditional. Row and column
// are one-based to match R indexing.
// [[Rcpp::export(name = ".hicpotts_cell_logweights_cpp")]]
NumericVector hicpotts_cell_logweights_cpp(
    NumericMatrix z, NumericMatrix y, NumericMatrix betas,
    NumericVector sizes, double theta, double gamma, std::string dist,
    NumericMatrix cov1, NumericMatrix cov2, NumericMatrix cov3,
    NumericMatrix cov4, int row, int column) {
  if (betas.nrow() != 3 || betas.ncol() != 5 || sizes.size() != 3)
    stop("betas must be 3 by 5 and sizes must have length three.");
  if (row < 1 || row > z.nrow() || column < 1 || column > z.ncol())
    stop("row and column are outside the lattice.");
  NumericMatrix l1 = log1p_matrix(cov1), l2 = log1p_matrix(cov2);
  NumericMatrix l3 = log1p_matrix(cov3), l4 = log1p_matrix(cov4);
  NumericVector beta_values[3] = {
    betas(0, _), betas(1, _), betas(2, _)
  };
  const double size_values[3] = {sizes[0], sizes[1], sizes[2]};
  double out[3];
  cell_logweights_cached(out, row - 1, column - 1, z, y, beta_values,
                         size_values, theta, gamma, dist, l1, l2, l3, l4);
  NumericVector result = NumericVector::create(out[0], out[1], out[2]);
  result.attr("names") = CharacterVector::create(
    "component1", "component2", "component3");
  return result;
}

// -----------------------------------------------------------------------------
// Theta (zero-inflation) update.
//
// Theta is updated through the standard latent-indicator augmentation of a
// ZIP/ZINB model. For every component-1 zero, u=1 denotes a structural zero:
//
// The model puts structural zero inflation on component 1 only, so conditional
// on z the only cells carrying information about theta are those with z == 1:
//
//   Pr(u=1 | y=0,z=1,theta) = theta / {theta + (1-theta) f1(0)}.
//
// Positive component-1 counts necessarily have u=0. With a Beta(1,1) prior,
// the next step is conjugate:
//   theta | u,z ~ Beta(1 + sum(u), 1 + n1 - sum(u)).
// Cells outside component 1 do not contain a zero-inflation indicator. The
// complete N x N field is used, consistently with every other emission term.
// -----------------------------------------------------------------------------
double theta_component1_loglik(const NumericMatrix &z_current, const NumericMatrix &y,
                               const NumericVector &b1,
                               const NumericMatrix &cov1, const NumericMatrix &cov2,
                               const NumericMatrix &cov3, const NumericMatrix &cov4,
                               double size1,
                               double theta_eval, const std::string &dist) {
  int N = z_current.nrow();
  double ll = 0.0;
  double log1m_theta = std::log1p(-theta_eval);

  // Full-matrix weighting, identical to the beta/dispersion path.
  const bool y_pairs_symmetric = hicpotts_symmetric_pairs(y);

  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      // Conditional on z, only component-1 cells inform theta.
      if ((int) z_current(i, j) != 1) continue;

      const double w = hicpotts_pair_weight(i, j, y_pairs_symmetric);

      double eta1 = b1[0] + b1[1] * std::log1p(cov1(i, j)) + b1[2] * std::log1p(cov2(i, j)) +
                    b1[3] * std::log1p(cov3(i, j)) + b1[4] * std::log1p(cov4(i, j));
      if (eta1 > 30.0) eta1 = 30.0; else if (eta1 < -30.0) eta1 = -30.0;
      double mu1 = std::exp(eta1);
      double yy = y(i, j);

      if (yy == 0.0) {
        double p0_base = (dist == "ZINB")
                           ? R::dnbinom_mu(0.0, size1, mu1, false)
                           : R::dpois(0.0, mu1, false);
        double mix = theta_eval + (1.0 - theta_eval) * p0_base;
        ll += w * std::log(mix > 1e-300 ? mix : 1e-300);
      } else {
        double log_dens = (dist == "ZINB")
                            ? R::dnbinom_mu(yy, size1, mu1, true)
                            : R::dpois(yy, mu1, true);
        ll += w * (log1m_theta + log_dens);
      }
    }
  }
  return ll;
}

// Internal testing hook: exposes the theta conditional log-likelihood so the
// test suite can verify detailed balance against the ACTUAL kernel rather than
// an R reimplementation of it. Not part of the user-facing API.
// [[Rcpp::export(.hicpotts_theta_loglik)]]
double hicpotts_theta_loglik(NumericMatrix z_current, NumericMatrix y,
                             NumericVector b1,
                             NumericMatrix cov1, NumericMatrix cov2,
                             NumericMatrix cov3, NumericMatrix cov4,
                             double size1, double theta_eval,
                             std::string dist) {
  return theta_component1_loglik(z_current, y, b1, cov1, cov2, cov3, cov4,
                                 size1, theta_eval, dist);
}

static double theta_gibbs_draw(
    const NumericMatrix &z_current, const NumericMatrix &y,
    const NumericVector &b1, const NumericMatrix &cov1,
    const NumericMatrix &cov2, const NumericMatrix &cov3,
    const NumericMatrix &cov4, double size1, double theta_current,
    const std::string &dist, bool covariates_are_logged = false) {
  int n_component1 = 0;
  int n_structural = 0;
  for (int i = 0; i < z_current.nrow(); ++i) {
    for (int j = 0; j < z_current.ncol(); ++j) {
      if ((int)z_current(i, j) != 1) continue;
      n_component1++;
      if (y(i, j) != 0.0) continue;
      const double x1 = covariates_are_logged ? cov1(i, j) : std::log1p(cov1(i, j));
      const double x2 = covariates_are_logged ? cov2(i, j) : std::log1p(cov2(i, j));
      const double x3 = covariates_are_logged ? cov3(i, j) : std::log1p(cov3(i, j));
      const double x4 = covariates_are_logged ? cov4(i, j) : std::log1p(cov4(i, j));
      double eta = b1[0] + b1[1] * x1 + b1[2] * x2 +
        b1[3] * x3 + b1[4] * x4;
      if (eta > 30.0) eta = 30.0;
      else if (eta < -30.0) eta = -30.0;
      const double mu = std::exp(eta);
      const double sampling_zero = dist == "ZINB" ?
        R::dnbinom_mu(0.0, size1, mu, false) : R::dpois(0.0, mu, false);
      const double denominator = theta_current +
        (1.0 - theta_current) * sampling_zero;
      const double structural_probability = denominator > 0.0 ?
        theta_current / denominator : 0.0;
      if (R::runif(0.0, 1.0) < structural_probability) n_structural++;
    }
  }
  return R::rbeta(1.0 + n_structural,
                  1.0 + n_component1 - n_structural);
}

// [[Rcpp::export(name = ".hicpotts_theta_gibbs_step_cpp")]]
double hicpotts_theta_gibbs_step_cpp(
    NumericMatrix z_current, NumericMatrix y, NumericVector b1,
    NumericMatrix cov1, NumericMatrix cov2, NumericMatrix cov3,
    NumericMatrix cov4, double size1, double theta_current,
    std::string dist) {
  if (b1.size() != 5) stop("b1 must contain five regression coefficients.");
  if (!(theta_current > 0.0 && theta_current < 1.0))
    stop("theta_current must lie strictly inside (0, 1).");
  if (dist != "ZIP" && dist != "ZINB")
    stop("The conjugate theta update requires ZIP or ZINB.");
  return theta_gibbs_draw(z_current, y, b1, cov1, cov2, cov3, cov4,
                          size1, theta_current, dist);
}

double update_theta_gibbs(const NumericMatrix &z_current, const NumericMatrix &y,
                          List chains, int iter,
                          const NumericMatrix &size_chain,
                          const NumericMatrix &cov1,
                          const NumericMatrix &cov2,
                          const NumericMatrix &cov3,
                          const NumericMatrix &cov4,
                          double theta_current, const std::string &dist) {
  NumericMatrix c1 = chains[0];
  NumericVector b1 = c1(iter, _);
  return theta_gibbs_draw(z_current, y, b1, cov1, cov2, cov3, cov4,
                          size_chain(0, iter), theta_current, dist, true);
}

// -----------------------------------------------------------------------------
// Gamma auxiliary-field simulation for the retained-state Gaussian-kernel ABC
// transition.
// -----------------------------------------------------------------------------
#include "abc_s1s2_gamma.h"

// Internal engine used by simulate_hicpotts_potts_truth(). It deliberately
// exposes the same equilibrated three-state field simulator used by the gamma
// transitions so recovery experiments have a known generating gamma.
// [[Rcpp::export(name = ".hicpotts_simulate_potts_labels_cpp")]]
NumericMatrix hicpotts_simulate_potts_labels_cpp(int N, double gamma,
                                                 int sweeps = 500) {
  if (N < 2) stop("N must be at least two.");
  if (!(gamma > 0.0 && gamma < 1.0) || !std::isfinite(gamma))
    stop("gamma must lie strictly between zero and one.");
  if (sweeps < 1) stop("sweeps must be positive.");
  return simulate_potts_abc(N, gamma, sweeps);
}

// Deterministic-calibration engine used once by repeated-chain R wrappers.
// It is intentionally internal: users control the public calibration through
// epsilon/abc_epsilon_quantile/abc_potts_sweeps/abc_sim_reps. The wrapper
// isolates this helper in the same fixed seed used by the standalone sampler.
// [[Rcpp::export(name = ".hicpotts_calibrate_abc_epsilon_cpp")]]
List hicpotts_calibrate_abc_epsilon_cpp(
    int N, double gamma_prior_shape1 = 1.0,
    double gamma_prior_shape2 = 1.0,
    double abc_epsilon_quantile = 0.10,
    int abc_potts_sweeps_arg = 0,
    int abc_sim_reps = 4,
    int abc_calibration_reps = 60) {
  if (N < 2) stop("N must be at least two.");
  if (!(gamma_prior_shape1 > 0.0) || !std::isfinite(gamma_prior_shape1) ||
      !(gamma_prior_shape2 > 0.0) || !std::isfinite(gamma_prior_shape2))
    stop("Gamma prior shapes must be finite and positive.");
  if (!(abc_epsilon_quantile > 0.0) ||
      !(abc_epsilon_quantile < 1.0) ||
      !std::isfinite(abc_epsilon_quantile))
    stop("abc_epsilon_quantile must lie strictly between zero and one.");
  if (abc_potts_sweeps_arg < 0)
    stop("abc_potts_sweeps_arg must be non-negative.");
  if (abc_sim_reps < 1 || abc_calibration_reps < 1)
    stop("ABC replicate counts must be positive.");
  const int sweeps = abc_potts_sweeps_arg > 0 ? abc_potts_sweeps_arg :
    std::max(100, (int) std::llround(
      100.0 * static_cast<double>(N * N) / 400.0));
  AbcPottsWorkspace workspace(N);
  const AbcToleranceCalibration calibration = calibrate_abc_tolerance(
    workspace, gamma_prior_shape1, gamma_prior_shape2,
    abc_epsilon_quantile, abc_calibration_reps, abc_sim_reps, sweeps);
  return List::create(
    Named("epsilon") = calibration.epsilon,
    Named("scale") = calibration.scale,
    Named("N") = N,
    Named("abc_potts_sweeps") = sweeps,
    Named("abc_sim_reps") = abc_sim_reps,
    Named("abc_calibration_reps") = abc_calibration_reps,
    Named("abc_epsilon_quantile") = abc_epsilon_quantile,
    Named("gamma_prior_shape1") = gamma_prior_shape1,
    Named("gamma_prior_shape2") = gamma_prior_shape2,
    Named("simulator_engine") = "compact_lookup_v1");
}


//' @title MCMC Chain Function for the HMRFHiC Model
//'
//' @description
//' The function integrates multiple steps, including proposal distributions, acceptance/rejection steps based on posterior calculations,
//' and retained-state ABC updates for the \eqn{\gamma} parameter. Covariates and prior specifications can be either
//' data-driven or user-defined. The output includes the chains of sampled parameters across the specified number of iterations.
//'
//' @details
//' Gamma uses a retained-state ABC-MCMC transition on a
//' composition-corrected latent-label adjacency summary, a Gaussian kernel,
//' and a default 10-percent prior-predictive tolerance quantile.
//'
//'
//' @param N Integer specifying the dimension of the lattice (\eqn{N \times N}).
//'
//' @param gamma_prior Numeric. The STARTING value for the interaction parameter
//'   \eqn{\gamma} in the Potts model. Note this sets only the initial state of
//'   the chain; the prior density used in the gamma acceptance ratio is the Beta
//'   distribution given by \code{gamma_prior_shape1} and
//'   \code{gamma_prior_shape2}.
//'
//' @param iterations Integer. The number of MCMC iterations to run.
//'
//' @param x_vars A list of covariates used as predictors in the model. This list should contain named elements
//'   corresponding to "distance", "GC", "TES", and "ACC", each containing a list of matrices of values for all \eqn{N \times N}.
//'
//' @param y An \eqn{N \times N} numeric matrix of observed interaction counts. Each element corresponds to interacting locus (i,j).
//'
//' @param use_data_priors Logical. If \code{TRUE}, data-driven priors will be used for each component.
//'   If \code{FALSE}, then \code{user_fixed_priors} must be provided.
//'
//' @param user_fixed_priors (Optional) A list of user-specified priors for the model components if
//'   \code{use_data_priors = FALSE}. Each components priors should be specified as a list of means and standard
//'   deviations for the parameters.
//'
//' @param dist A character string specifying the distribution family to use:
//'   \itemize{
//'     \item \code{"Poisson"}: Poisson distribution
//'     \item \code{"NB"}: Negative Binomial distribution
//'     \item \code{"ZIP"}: Zero-Inflated Poisson distribution
//'     \item \code{"ZINB"}: Zero-Inflated Negative Binomial distribution
//'   }
//'   The default is \code{"ZIP"}.
//'
//' @param epsilon Optional positive ABC bandwidth. If NULL, a deterministic
//'   prior-predictive pilot quantile is used. Supplying a value skips the pilot
//'   simulations entirely.
//'
//' @param distance_metric Retained for source compatibility and ignored.
//'
//' @param size_start (Optional) Required if \code{dist} is "NB" or "ZINB". A numeric vector of length 3 specifying
//'   initial values for the size (overdispersion) parameter for each of the 3 components.
//'
//' @param theta_start (Optional) A numeric value providing an initial value for \eqn{\theta}, the zero-inflation
//'   parameter in ZIP/ZINB models. If not provided and \code{dist} is ZIP or ZINB, \eqn{\theta} is initialized to 0.5.
//'
//' @param z_start Optional N-by-N starting allocation with labels 1, 2 and 3.
//' @param beta_start Optional 3-by-5 matrix of starting regression
//'   coefficients. This is primarily useful for known-truth validation.
//' @param validation_freeze_z,validation_freeze_beta,validation_freeze_gamma
//'   Simulation-diagnostic controls that hold the named state block at its
//'   supplied starting value. All default to false and therefore do not
//'   alter production inference. A frozen beta block requires
//'   \code{beta_start}.
//' @param validation_freeze_size,validation_freeze_theta
//'   Simulation-diagnostic controls for the dispersion and zero-inflation
//'   blocks. Frozen dispersion requires \code{size_start} and frozen theta
//'   requires \code{theta_start} under the corresponding emission family.
//' @param mcse_stop Whether to stop when every monitored parameter satisfies
//'   the batch-means relative MCSE rule.
//' @param mcse_min_iterations Minimum updates before an MCSE check.
//' @param mcse_check_interval Updates between MCSE checks.
//' @param mcse_relative_threshold Maximum MCSE/posterior-SD ratio required
//'   for stopping before the user-requested iteration count.
//' @param tempering_warmup Optional cyclical heating updates for latent-state
//'   warm-up. Zero is the tested default; retained sampling uses temperature 1.
//' @param tempering_beta_min Lowest optional warm-up inverse temperature.
//' @param tempering_cycle Positive optional heating-cycle length.
//'
//' @param gamma_prior_shape1,gamma_prior_shape2 Positive shape parameters of the
//'   Beta prior on \eqn{\gamma} used in either gamma transition, and for the
//'   legacy ABC tolerance calibration. Default \code{1, 1}, i.e.
//'   Beta(1,1) = Uniform(0, 1): a flat prior whose density is constant, so it
//'   contributes exactly zero to the acceptance ratio and pulls \eqn{\gamma}
//'   toward no particular value.
//'
//'   The previous behaviour was a hardcoded Beta(2,2), whose mean of 0.5 sits
//'   above any realistic Potts coupling and biased \eqn{\gamma} upward whenever
//'   the tolerance was loose.
//'
//'   In a screening study over nine simulated datasets at
//'   \eqn{\gamma \in \{0.15, 0.30, 0.45\}}, the flat prior had the lowest
//'   worst-case absolute bias across \eqn{\gamma} at the default tolerance
//'   (0.032, against 0.045 for Beta(2,2) and 0.066 for Beta(1,3)), and was by
//'   far the most stable across tolerance choices: its worst-case bias ranged
//'   only 0.032-0.049 over quantiles 0.05-0.30, where Beta(2,2) degraded to
//'   0.155 and Beta(1,3) to 0.120. That stability means results depend far
//'   less on the tolerance setting. The cost is wider intervals (0.268 versus
//'   0.200 for Beta(1,3) at quantile 0.10), which is the honest price of not
//'   injecting prior information.
//'
//'   Informative alternatives remain available: Beta(1,3) (mean 0.25) if weak
//'   coupling is genuinely expected, or Beta(2,2) to reproduce results from
//'   versions before 1.3.7. Note the study covered \eqn{\gamma} only up to
//'   0.45 and used three datasets per cell, so per-cell estimates are noisy;
//'   it is a screening study, not simulation-based calibration.
//'
//' @param abc_epsilon_quantile Numeric in (0, 1). Quantile of the
//'   prior-predictive summary distances used as the ABC tolerance when
//'   \code{epsilon} is not supplied. Default \code{0.10}; the previous
//'   hardcoded value was \code{0.30}, which produced a tolerance comparable to
//'   the between-\eqn{\gamma} signal so that the prior dominated the posterior.
//'   Smaller values sharpen the ABC approximation at the cost of a lower
//'   acceptance rate. Ignored when \code{epsilon} is supplied.
//'
//' @param gamma_update_interval Positive integer. The selected gamma step - the
//'   dominant per-iteration cost - runs only every \code{gamma_update_interval}
//'   iterations; on the others gamma is held fixed while beta/size/theta/z
//'   continue to update. Default \code{5}. Setting it to
//'   \code{k} divides the gamma cost by ~\code{k} and thins the gamma chain by
//'   \code{k}. Inspect retained gamma movement when choosing this interval.
//'   Per-iteration cost of the gamma step is
//'   approximately \code{abc_sim_reps * abc_potts_sweeps * N^2}.
//'
//' @param abc_potts_sweeps_arg Integer. Number of Gibbs sweeps per simulated
//'   Potts field. \code{0} (default) uses the N-dependent rule
//'   \code{max(100, 100 * N^2 / 400)}. A positive value overrides it. The
//'   default is conservative for weak coupling (gamma well below the critical
//'   value ~1); reducing it speeds simulation but, if set too low, the
//'   simulated field is not stationary and the summary-vs-gamma relationship
//'   shifts, re-biasing gamma. Validate against the estimating equation before
//'   reducing.
//'
//' @param abc_sim_reps Positive integer; number of simulated fields averaged
//'   per ABC candidate.
//' @param z_probability_burnin_arg Integer iteration after which latent-state
//'   membership frequencies are accumulated for classification. The default
//'   \code{-1} selects an automatic value. This controls returned summaries
//'   only and does not change an MCMC transition or target distribution.
//'
//' @param comp23_barrier_kappa Strength of the smooth crossing-barrier penalty
//'   that discourages the component 2 and 3 baseline intensities from sitting
//'   on top of one another (the mechanism behind label-switching flicker).
//'   The default is \code{10}, so this separation term and the component-1/3
//'   biological relationship prior are both active.
//' @param comp23_barrier_w Width of that penalty, i.e. how close the two
//'   intercepts must be before it bites. The default \code{0.3} reproduces the
//'   historical hardcoded value. Must be positive.
//' @param use_noise_relationship_prior Logical. If \code{TRUE} (default),
//'   component 3 is coupled to approximately the component-1 slopes and softly
//'   preferred above component 1 in intercept. Component 2 is absent from this
//'   coupled term and remains the unrestricted true-signal component.
//' @param noise_link_sd Positive standard deviation controlling shrinkage of
//'   component-1 versus false-signal slope differences on the standardised
//'   log1p-covariate scale. Smaller values enforce closer slopes.
//' @param noise_order_strength Non-negative strength of the soft preference
//'   that false-signal noise has a higher intercept than component 1.
//' @param noise_order_width Positive transition width of that soft intercept
//'   preference on the log-mean scale.
//' @param branch_swap_interval Non-negative integer. When positive, attempt
//'   the reversible whole component-2/component-3 branch transition at this
//'   interval. The default is every iteration; zero disables it. The move is active only with the labelled
//'   component-1/3 relationship prior.
//' @param signal_block_move_interval Non-negative integer interval for the
//'   connected component-2/component-3 cluster move. Zero disables it.
//' @param gamma_large_jump_probability Conditional probability of using the
//'   wider logit-normal gamma proposal rather than the adaptive local proposal.
//' @param gamma_large_jump_multiplier Positive multiplier applied to the local
//'   logit proposal scale for a wide jump.
//' @param gamma_independence_probability Probability of proposing gamma
//'   independently from its configured Beta prior. The complete mixture
//'   density is used in the Metropolis-Hastings ratio.
//' @param gamma_method Retained for source compatibility. The only supported
//'   value is \code{"abc"}; the exchange-style transition has been removed.
//' @param verbose Logical. Print periodic sampler progress when true. The
//'   default is false so parallel and batch execution remain quiet.
//' @param progress_interval Positive number of iterations between progress
//'   messages when \code{verbose = TRUE}.
//'
//' @return
//' A \code{list} containing:
//' \itemize{
//'   \item \code{chains}: A list of three \eqn{(iterations+1) \times 5} matrices, each storing the parameter chains
//'         for one of the three components.
//'   \item \code{gamma}: A numeric vector of length \eqn{iterations+1}, the chain of \eqn{\gamma} values.
//'   \item \code{theta}: A numeric vector of length \eqn{iterations+1}, the chain of \eqn{\theta} values
//'         (or initialized and unchanged if not ZIP/ZINB).
//'   \item \code{size}: A \eqn{3 \times (iterations+1)} matrix of size (overdispersion) parameters
//'         if \code{dist} is NB or ZINB; otherwise, this matrix might be unused.
//'   \item \code{z_final}: Final unrestricted latent-state draw.
//'   \item \code{z_checkpoints}: Unrestricted latent-state snapshots retained
//'         for trajectory diagnostics.
//'   \item \code{z_probabilities}: Three \eqn{N \times N} matrices containing
//'         post-burn-in MCMC membership frequencies. These directly estimate
//'         the posterior probability of each component at every cell without
//'         changing the sampler or applying a mean-field approximation.
//'   \item \code{z_probability_draws}: Number of latent-state draws accumulated
//'         in \code{z_probabilities}.
//'   \item \code{performance}: Cached-native engine flags plus total and
//'         stage-specific elapsed times.
//'   \item \code{sampler_settings}: Dimensions, distribution, completed
//'         iterations and computational settings needed to audit the fit.
//' }
//'
//' @examples
//'
//' N <- 10
//' gamma_prior <- 0.5
//' iterations <- 10
//' x_vars <- list(
//'   distance = list(matrix(runif(N*N, 0, 10), nrow=N)),
//'   GC = list(matrix(runif(N*N, 0, 1), nrow=N)),
//'   TES = list(matrix(runif(N*N, 0, 2), nrow=N)),
//'   ACC = list(matrix(runif(N*N, 0, 5), nrow=N))
//' )
//' y <- matrix(rpois(N*N, lambda=5), nrow=N)
//'
//' # Using data-driven priors, ZIP distribution, no user_fixed_priors:
//' results <- run_metropolis_MCMC_betas(
//'   N = N,
//'   gamma_prior = gamma_prior,
//'   iterations = iterations,
//'   x_vars = x_vars,
//'   y = y,
//'   use_data_priors = TRUE,
//'   dist = "ZIP"
//' )
//'
//' # Inspect gamma chain:
//' plot(results$gamma, type='l')
//'
//' @aliases run_metropolis_MCMC_betas
//'
//' @noRd
// [[Rcpp::export]]
List run_metropolis_MCMC_betas(int N, double gamma_prior, int iterations,
                               List x_vars, NumericMatrix y,
                               bool use_data_priors, Nullable<List> user_fixed_priors = R_NilValue,
                               std::string dist = "ZIP", Nullable<double> epsilon = R_NilValue,
                               std::string distance_metric = "manhattan",
                               Nullable<NumericVector> size_start = R_NilValue, Nullable<double> theta_start = R_NilValue,
                               Nullable<NumericMatrix> z_start = R_NilValue,
                               bool mcse_stop = true,
                               int mcse_min_iterations = 10000,
                               int mcse_check_interval = 500,
                               double mcse_relative_threshold = 0.05,
                               int tempering_warmup = 0,
                               double tempering_beta_min = 0.30,
                               int tempering_cycle = 500,
                               double gamma_prior_shape1 = 1.0,
                               double gamma_prior_shape2 = 1.0,
                               double abc_epsilon_quantile = 0.10,
                               int gamma_update_interval = 5,
                               int abc_potts_sweeps_arg = 0,
                               int abc_sim_reps = 4,
                               int z_probability_burnin_arg = -1,
                               double comp23_barrier_kappa = 10.0,
                               double comp23_barrier_w = 0.3,
                               bool use_noise_relationship_prior = true,
                               double noise_link_sd = 0.5,
                               double noise_order_strength = 10.0,
                               double noise_order_width = 0.5,
                               int branch_swap_interval = 1,
                               int signal_block_move_interval = 5,
                               double gamma_large_jump_probability = 0.10,
                               double gamma_large_jump_multiplier = 4.0,
                               double gamma_independence_probability = 0.05,
                               std::string gamma_method = "abc",
                               bool verbose = false,
                               int progress_interval = 50,
                               Nullable<NumericMatrix> beta_start = R_NilValue,
                               bool validation_freeze_z = false,
                               bool validation_freeze_beta = false,
                               bool validation_freeze_gamma = false,
                               bool validation_freeze_size = false,
                               bool validation_freeze_theta = false) {
  if (N < 2) stop("N must be at least two.");
  if (iterations < 1)
    stop("iterations must be a positive integer.");
  if (!(gamma_prior > 0.0 && gamma_prior < 1.0) || !std::isfinite(gamma_prior))
    stop("gamma_prior must lie strictly between zero and one.");
  if (dist != "Poisson" && dist != "NB" && dist != "ZIP" && dist != "ZINB")
    stop("dist must be Poisson, NB, ZIP or ZINB.");
  if (progress_interval < 1)
    stop("progress_interval must be a positive integer.");
  // Barrier controls (see the crossing-barrier note at the top of this file).
  // kappa = 0 disables the penalty entirely and is an explicitly supported
  // sensitivity setting, so only negative/non-finite values are rejected.
  if (!(comp23_barrier_kappa >= 0.0) || !std::isfinite(comp23_barrier_kappa))
    stop("comp23_barrier_kappa must be a non-negative finite number (0 disables the barrier).");
  if (!(comp23_barrier_w > 0.0) || !std::isfinite(comp23_barrier_w))
    stop("comp23_barrier_w must be a positive finite number.");
  if (!(noise_link_sd > 0.0) || !std::isfinite(noise_link_sd))
    stop("noise_link_sd must be a positive finite number.");
  if (!(noise_order_strength >= 0.0) || !std::isfinite(noise_order_strength))
    stop("noise_order_strength must be a non-negative finite number.");
  if (!(noise_order_width > 0.0) || !std::isfinite(noise_order_width))
    stop("noise_order_width must be a positive finite number.");
  if (branch_swap_interval < 0 || signal_block_move_interval < 0)
    stop("branch and signal-block move intervals must be non-negative.");
  if (!(gamma_large_jump_probability >= 0.0 &&
        gamma_large_jump_probability <= 1.0) ||
      !std::isfinite(gamma_large_jump_probability))
    stop("gamma_large_jump_probability must lie in [0, 1].");
  if (!(gamma_large_jump_multiplier > 1.0) ||
      !std::isfinite(gamma_large_jump_multiplier))
    stop("gamma_large_jump_multiplier must be finite and greater than one.");
  if (!(gamma_independence_probability >= 0.0 &&
        gamma_independence_probability < 1.0) ||
      !std::isfinite(gamma_independence_probability))
    stop("gamma_independence_probability must lie in [0, 1).");
  if (gamma_method != "abc")
    stop("gamma_method='exchange' has been removed; use gamma_method='abc'.");
  if (!(gamma_prior_shape1 > 0.0) || !std::isfinite(gamma_prior_shape1) ||
      !(gamma_prior_shape2 > 0.0) || !std::isfinite(gamma_prior_shape2))
    stop("gamma_prior_shape1 and gamma_prior_shape2 must be positive and finite.");
  if (!(abc_epsilon_quantile > 0.0) || !(abc_epsilon_quantile < 1.0) ||
      !std::isfinite(abc_epsilon_quantile))
    stop("abc_epsilon_quantile must lie strictly between 0 and 1.");
  if (gamma_update_interval < 1)
    stop("gamma_update_interval must be a positive integer (1 = every iteration).");
  if (abc_potts_sweeps_arg < 0)
    stop("abc_potts_sweeps must be >= 1, or 0 to use the N-dependent default.");
  if (abc_sim_reps < 1)
    stop("abc_sim_reps must be a positive integer.");
  if (z_probability_burnin_arg < -1 ||
      z_probability_burnin_arg >= iterations)
    stop("z_probability_burnin_arg must be -1 (automatic) or an integer in [0, iterations - 1].");
  if (mcse_min_iterations < 1 || mcse_check_interval < 1)
    stop("MCSE minimum iterations and check interval must be positive.");
  if (!(mcse_relative_threshold > 0.0) || !std::isfinite(mcse_relative_threshold))
    stop("mcse_relative_threshold must be finite and positive.");
  if (tempering_warmup < 0 || tempering_warmup > iterations)
    stop("tempering_warmup must be between zero and iterations.");
  if (!(tempering_beta_min > 0.0 && tempering_beta_min <= 1.0) ||
      !std::isfinite(tempering_beta_min) || tempering_cycle < 1)
    stop("tempering_beta_min must be in (0,1] and tempering_cycle positive.");
  if (y.nrow() != N || y.ncol() != N)
    stop("y must be an N-by-N matrix.");
  for (int j = 0; j < N; ++j) {
    for (int i = 0; i < N; ++i) {
      const double value = y(i, j);
      if (!std::isfinite(value) || value < 0.0 ||
          std::fabs(value - std::round(value)) > 1e-8)
        stop("y must contain finite non-negative integer counts.");
    }
  }
  if (!x_vars.containsElementNamed("distance") ||
      !x_vars.containsElementNamed("GC") ||
      !x_vars.containsElementNamed("TES") ||
      !x_vars.containsElementNamed("ACC"))
    stop("x_vars must contain distance, GC, TES and ACC.");

   // Map x1, x2, x3, and x4 to distance, GC, TES, and ACC
  List x1 = x_vars["distance"];
  List x2 = x_vars["GC"];
  List x3 = x_vars["TES"];
  List x4 = x_vars["ACC"];

  // Covariate matrices, hoisted out of the main loop since x_vars is fixed
  // input data that does not change across iterations. Reused by the theta
  // update and the gamma ABC update.
  NumericMatrix x11 = as<NumericMatrix>(x1[0]);
  NumericMatrix x22 = as<NumericMatrix>(x2[0]);
  NumericMatrix x33 = as<NumericMatrix>(x3[0]);
  NumericMatrix x44 = as<NumericMatrix>(x4[0]);
  if (x11.nrow() != N || x11.ncol() != N ||
      x22.nrow() != N || x22.ncol() != N ||
      x33.nrow() != N || x33.ncol() != N ||
      x44.nrow() != N || x44.ncol() != N)
    stop("Every covariate matrix must be N by N.");

  // Cache every transformed covariate exactly once. These matrices are shared
  // by all native likelihood and allocation paths below.
  NumericMatrix logx11 = log1p_matrix(x11);
  NumericMatrix logx22 = log1p_matrix(x22);
  NumericMatrix logx33 = log1p_matrix(x33);
  NumericMatrix logx44 = log1p_matrix(x44);

  // Standardised proposal coordinates for better conditioning. These values
  // never replace x_vars in the likelihood, and all returned coefficients are
  // back-transformed to the original log1p-covariate scale.
  arma::vec covariate_means(4), covariate_sds(4);
  cached_mean_sd(logx11, covariate_means[0], covariate_sds[0]);
  cached_mean_sd(logx22, covariate_means[1], covariate_sds[1]);
  cached_mean_sd(logx33, covariate_means[2], covariate_sds[2]);
  cached_mean_sd(logx44, covariate_means[3], covariate_sds[3]);

  // Thin QR whitening of [1, standardized log1p covariates]. Scaling R by
  // sqrt(n) keeps proposal coordinates comparable to standardized betas.
  // Only the proposal geometry changes; the target remains on beta.
  const int design_n = N * N;
  arma::mat design(design_n, 5, arma::fill::ones);
  int design_row = 0;
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++, design_row++) {
      design(design_row, 1) = (logx11(i, j) - covariate_means[0]) / covariate_sds[0];
      design(design_row, 2) = (logx22(i, j) - covariate_means[1]) / covariate_sds[1];
      design(design_row, 3) = (logx33(i, j) - covariate_means[2]) / covariate_sds[2];
      design(design_row, 4) = (logx44(i, j) - covariate_means[3]) / covariate_sds[3];
    }
  }
  arma::mat Q_design, R_scaled;
  bool qr_whitened = arma::qr_econ(Q_design, R_scaled, design) &&
    R_scaled.n_rows == 5 && R_scaled.n_cols == 5;
  if (qr_whitened) {
    R_scaled /= std::sqrt((double)design_n);
    for (int j = 0; j < 5; j++) {
      if (R_scaled(j, j) < 0.0) R_scaled.row(j) *= -1.0;
      if (std::fabs(R_scaled(j, j)) < 1e-10) qr_whitened = false;
    }
  }
  if (!qr_whitened) {
    // Collinearity affects proposal efficiency, not model validity. Identity
    // makes the helpers use centred/standardised coordinates without QR.
    R_scaled.eye(5, 5);
  }

  // epsilon is the ABC Gaussian-kernel bandwidth on the calibrated
  // join-count-excess summary statistic. distance_metric remains accepted
  // for source compatibility.
  double abc_user_epsilon = epsilon.isNotNull() ? as<double>(epsilon) : NA_REAL;
  if (epsilon.isNotNull() && (!(abc_user_epsilon > 0.0) || !std::isfinite(abc_user_epsilon)))
    stop("epsilon must be NULL or a finite positive ABC bandwidth.");
  (void) distance_metric;

  // Check if user_fixed_priors is needed
  if (!use_data_priors && user_fixed_priors.isNull()) {
    stop("Error: When use_data_priors is set to false, user_fixed_priors must be provided.");
  }

  // Initialize chains and other variables
  List chains = List::create(NumericMatrix(iterations + 1, 5),
                             NumericMatrix(iterations + 1, 5),
                             NumericMatrix(iterations + 1, 5));
  if (validation_freeze_beta && beta_start.isNull())
    stop("validation_freeze_beta requires a 3-by-5 beta_start matrix.");
  if (beta_start.isNotNull()) {
    NumericMatrix initial_beta(beta_start);
    if (initial_beta.nrow() != 3 || initial_beta.ncol() != 5)
      stop("beta_start must be a 3-by-5 matrix.");
    for (int component = 0; component < 3; ++component) {
      NumericMatrix component_chain = chains[component];
      for (int coefficient = 0; coefficient < 5; ++coefficient) {
        const double value = initial_beta(component, coefficient);
        if (!std::isfinite(value)) stop("beta_start must contain finite values.");
        component_chain(0, coefficient) = value;
      }
    }
  }
  NumericVector chain_gamma(iterations + 1, NA_REAL);
  chain_gamma[0] = gamma_prior; // Initialize with gamma_prior
  NumericVector theta(iterations + 1, NA_REAL);
  // Initialize theta with theta_start if provided and if dist is ZIP or ZINB
  if (theta_start.isNotNull() && (dist == "ZIP" || dist == "ZINB")) {
    theta[0] = as<double>(theta_start);
  } else {
    theta[0] = 0.5; // Default initialization for theta if not using ZIP or ZINB
  }
  if ((dist == "ZIP" || dist == "ZINB") &&
      (!(theta[0] > 0.0 && theta[0] < 1.0) || !std::isfinite(theta[0])))
    stop("theta_start must lie strictly between zero and one.");
  if (validation_freeze_theta &&
      (dist == "ZIP" || dist == "ZINB") && theta_start.isNull())
    stop("validation_freeze_theta requires theta_start for ZIP/ZINB.");


  // Initialize acceptance counts and QR-coordinate proposal standard
  // deviations for components 1, 2, and 3.  These proposal scales remain in
  // force for the complete run; they are adapted only before the retained
  // sample boundary and are then frozen.
  std::vector<int> acceptance_counts = {0, 0, 0};
  std::vector<long> beta_attempt_total = {0L, 0L, 0L};
  std::vector<long> beta_accept_total = {0L, 0L, 0L};
  std::vector<long> beta_warmup_attempt_total = {0L, 0L, 0L};
  std::vector<long> beta_warmup_accept_total = {0L, 0L, 0L};
  List sd_values = List::create(NumericVector::create(1.0, 0.5, 0.5, 0.5, 0.5),
                                NumericVector::create(1.0, 0.5, 0.5, 0.5, 0.5),
                                NumericVector::create(1.0, 0.5, 0.5, 0.5, 0.5));
  const int beta_dimension = 5;

  // Log-scale adaptive random-walk step for the overdispersion (size) parameter,
  // per component (NB / ZINB only). Proposing log(size) keeps size positive
  // and is scale-invariant.
  std::vector<double> size_logstep(3, 0.3);   // sd of the N(0, .) step on log(size)
  std::vector<int>    size_accept(3, 0);      // size acceptances in the current window
  std::vector<long>   size_attempt_total(3, 0L);
  std::vector<long>   size_accept_total(3, 0L);

  // Gamma uses a mixture of local, wide and prior-independent proposals on
  // (0, 1) inside the retained-state ABC transition.
  double gamma_logit_step = 0.5;
  int    gamma_local_accept_window = 0;
  long   gamma_local_proposals_window = 0;
  long   gamma_proposals_total = 0;
  long   gamma_local_proposals_total = 0;
  long   gamma_large_proposals_total = 0;
  long   gamma_independence_proposals_total = 0;
  long   gamma_accept_total = 0;
  long   gamma_local_accept_total = 0;
  long   gamma_large_accept_total = 0;
  long   gamma_independence_accept_total = 0;
  // Auxiliary Potts equilibration sweeps scale with lattice size (proportional
  // to the number of cells, baselined at 100 sweeps for the 20x20 lattice) so
  // that simulated fields reach the same clustering level regardless of N.
  // Under-equilibrated auxiliary fields are under-clustered, which biases gamma
  // upward; scaling removes that artifact at larger lattices.
  // Number of Potts sweeps per simulated field. 0 (the argument default) keeps
  // the N-dependent rule; a positive argument overrides it. This and
  // abc_sim_reps are the per-iteration ABC cost drivers (cost ~ abc_sim_reps *
  // abc_potts_sweeps * N^2). See gamma_update_interval for thinning the whole
  // gamma step.
  const int abc_potts_sweeps = (abc_potts_sweeps_arg > 0)
      ? abc_potts_sweeps_arg
      : std::max(100, (int) std::llround(100.0 * (double)(N * N) / 400.0));
  AbcPottsWorkspace abc_potts_workspace(N);
  const int abc_calibration_reps = 60;
  // Fixed RNG substream for tolerance calibration, so epsilon does not vary
  // with the chain's seed. Deliberately a constant: the ABC tolerance defines
  // which approximate posterior is being targeted, and that should not differ
  // between chains of the same analysis. Supply `epsilon` to pin it directly.
  const int abc_calibration_seed = 20260721;
  bool abc_calibrated = false;
  int abc_calibration_reps_executed = 0;
  double abc_scale = NA_REAL, abc_epsilon = NA_REAL, T_obs_fixed = NA_REAL;
  // Augmented ABC state: the auxiliary summary currently PAIRED with
  // chain_gamma[iter], and its distance to the observed summary. These persist
  // across iterations and are only replaced when a proposal is accepted.
  // Retained auxiliary state: the sufficient statistic of a Potts field drawn
  // at the current gamma. The simulator p(z | gamma) depends on gamma and N
  // alone, so this retained draw stays valid for the whole run.
  double abc_T_state = NA_REAL, abc_d_state = NA_REAL;
  NumericVector abc_distance_chain(iterations + 1, NA_REAL);
  NumericVector abc_s1_chain(iterations + 1, NA_REAL);
  NumericVector abc_s2_chain(iterations + 1, NA_REAL);

  // Crossing-barrier strength (kappa, barrier height) and width (w) for the
  // components 2/3 soft-separation penalty -- see crossing_barrier()'s
  // definition above. A width of 0.3 is narrower than the typical
  // signal-vs-false-signal separation on the log scale, so the barrier is
  // essentially inert once the chain has settled away from a tie in either
  // direction; a strength of 10 is moderate relative to typical per-iteration
  // log-posterior fluctuations, so it discourages lingering near a tie
  // without dominating the likelihood.
  // Strength (kappa) and width (w) now arrive as arguments; the defaults above
  // in the signature reproduce the previous hardcoded 10.0 / 0.3 exactly, so
  // existing results are unchanged unless a caller opts in to different
  // values. Both are recorded on the returned fit object.
  // Optional diagnostic toggles (environment-variable gated, off by default):
  // allow the crossing-barrier penalty and/or the non-local cluster move to
  // be disabled independently, for testing their individual contribution to
  // component 2/3 identifiability.
  const bool ablation_disable_blockmove  = std::getenv("HICPOTTS_DISABLE_BLOCKMOVE")  != nullptr;
  // The legacy HICPOTTS_DISABLE_BARRIER toggle is now expressed through the
  // public argument: setting the environment variable forces kappa to 0, which
  // is exactly what the old boolean did. Keeping the variable preserves any
  // existing ablation scripts.
  if (std::getenv("HICPOTTS_DISABLE_BARRIER") != nullptr) comp23_barrier_kappa = 0.0;
  const bool ablation_disable_barrier = !(comp23_barrier_kappa > 0.0);
  const bool branch_swap_enabled =
    use_noise_relationship_prior && branch_swap_interval > 0 &&
    !validation_freeze_z && !validation_freeze_beta &&
    !validation_freeze_size;
  int branch_swap_attempts = 0;
  int branch_swap_accepts = 0;
  int branch_swap_toward_intended_attempts = 0;
  int branch_swap_toward_intended_accepts = 0;
  int branch_swap_warmup_attempts = 0;
  int branch_swap_warmup_accepts = 0;
  int signal_block_move_attempts = 0;
  int signal_block_move_accepts = 0;
  int signal_block_warmup_attempts = 0;
  int signal_block_warmup_accepts = 0;
  NumericVector branch_gap_chain(iterations + 1, NA_REAL);

  // Initialize size_chain as a matrix for each component
  NumericMatrix size_chain(3, iterations + 1);  // 3 components, each with a chain of size `iterations + 1`
  if (dist == "NB" || dist == "ZINB") {
    if (size_start.isNotNull()) {
      NumericVector init_sizes = Rcpp::as<NumericVector>(size_start.get());
      if (init_sizes.size() != 3) stop("size_start should have 3 initial values for 3 components");
      for (int c = 0; c < 3; c++) {
        if (!(init_sizes[c] > 0.0) || !std::isfinite(init_sizes[c]))
          stop("Every size_start value must be positive and finite.");
        size_chain(c, 0) = init_sizes[c];
      }
    } else {
      stop("size_start must be provided for Negative Binomial or ZINB distribution.");
    }
  }
  if (validation_freeze_size && !(dist == "NB" || dist == "ZINB"))
    stop("validation_freeze_size is only meaningful for NB/ZINB.");

  // Adaptive tuning parameters for MCMC. The adaptation window scales with
  // the run length so that short runs still retune their proposal step sizes.
  int adaptation_interval = std::max(10, std::min(50, iterations / 10));
  double adaptation_scaling = 1.2;
  // Target acceptance band (around the ~0.22 optimum for these blocks). The step
  // size is only changed when the rate falls outside the band, avoiding the
  // oscillation caused by adjusting on every window.
  double accept_lo = 0.15;
  double accept_hi = 0.30;
  // The one-dimensional log-dispersion moves can use a wider optimal band
  // than the five-dimensional beta blocks.
  double size_accept_lo = 0.20;
  double size_accept_hi = 0.40;


  // Initialize z. All production likelihood and allocation calculations stay
  // inside the cached native engine; the public R helpers remain available for
  // backwards compatibility and independent equivalence tests.
  RObject fixed_prior_object = user_fixed_priors.isNotNull() ?
    RObject(user_fixed_priors.get()) : RObject(R_NilValue);

  // Initialise latent states for parameter inference. Sampling zeros are
  // possible under every component likelihood; observed zeros are forced to
  // component 1 only in the reported classification output below, not here.
  NumericMatrix z_current (N,N);
  // Z-warmup: when an informative z_start is supplied, the latent field is
  // held fixed at it for the first z_warmup updates while the regression,
  // dispersion and theta parameters equilibrate to that partition, before z
  // is released. This gives the parameters time to reach reasonable values
  // before the first latent-field update, rather than starting from
  // zero-initialised parameters that could erode a good initial partition
  // (particularly for sparse minority components). It changes only the
  // initial burn-in trajectory; once z is released, the sampler targets the
  // ordinary joint posterior. It is disabled automatically when no z_start
  // is supplied, since there is then no informative field to protect.
  const bool has_zstart = z_start.isNotNull();
  const int z_warmup = has_zstart ? std::min(600, iterations / 3) : 0;
  if (z_start.isNotNull()) {
    NumericMatrix z_initial(z_start);
    if (z_initial.nrow() != N || z_initial.ncol() != N)
      stop("z_start must be an N-by-N matrix.");
    for (int i = 0; i < N; i++) {
      for (int j = 0; j < N; j++) {
        double value = z_initial(i, j);
        if (!std::isfinite(value) || value < 1.0 || value > 3.0 ||
            std::fabs(value - std::round(value)) > 1e-8)
          stop("z_start values must be integer component labels 1, 2 or 3.");
        z_current(i, j) = std::round(value);
      }
    }
  } else {
    for (int i = 0; i < N; i++) {
      for (int j = 0; j < N; j++) {
        z_current(i, j) = std::floor(R::runif(1, 4));
      }
    }
  }

  // Classification is a posterior summary, not a new transition rule. Keep
  // online site-wise state counts after an automatic burn-in so downstream
  // calls can use the latent draws generated by this exact sampler instead of
  // reconstructing membership with a plug-in or mean-field approximation.
  // The default production run (20,000 iterations, MCSE checks after 10,000)
  // therefore starts accumulating at 5,000, matching fit_hicpotts_robust().
  const int automatic_z_probability_burnin =
    std::min(iterations, mcse_min_iterations) / 2;
  const int requested_z_probability_burnin =
    z_probability_burnin_arg < 0 ?
      automatic_z_probability_burnin : z_probability_burnin_arg;
  const int z_probability_burnin =
    std::max(z_warmup, requested_z_probability_burnin);

  // Soft-allocation empirical-Bayes state. Conditional membership
  // probabilities are refreshed during warm-up and converted to one fixed set
  // of Normal hyperparameters at the retained-sample boundary. The initial
  // partition supplies one-hot weights until the first released z sweep.
  NumericMatrix eb_weight1(N, N), eb_weight2(N, N), eb_weight3(N, N);
  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < N; ++j) {
      const int label = (int) z_current(i, j);
      eb_weight1(i, j) = label == 1 ? 1.0 : 0.0;
      eb_weight2(i, j) = label == 2 ? 1.0 : 0.0;
      eb_weight3(i, j) = label == 3 ? 1.0 : 0.0;
    }
  }
  RObject active_prior_object = fixed_prior_object;
  int eb_prior_updates = 0;
  if (use_data_priors) {
    active_prior_object = soft_eb_priors_impl(
      y, logx11, logx22, logx33, logx44,
      eb_weight1, eb_weight2, eb_weight3, true);
  }
  const int eb_freeze_iteration = z_probability_burnin;
  NumericMatrix z_probability_count1(N, N);
  NumericMatrix z_probability_count2(N, N);
  NumericMatrix z_probability_count3(N, N);
  int z_probability_draws = 0;

  // ---- batched latent-state frequencies -----------------------------------
  // Allocation draws are serially correlated, so the naive standard error
  // p(1-p)/sqrt(n) understates the Monte-Carlo error of a per-cell membership
  // probability, sometimes by an order of magnitude. Batch means fix that
  // without storing every draw: post-burn-in draws are split into contiguous
  // batches, each batch's mean membership is retained, and the spread BETWEEN
  // batch means gives an autocorrelation-aware MCSE and an effective sample
  // size per cell.
  //
  // Batch size follows the standard sqrt(n) rule against the number of
  // post-burn-in draws the run is expected to produce. The batch count is
  // capped so memory stays O(cap * 3 * N^2) rather than growing with the run
  // length: at the cap this is 32 * 3 * N^2 doubles.
  const int expected_z_draws = std::max(1, iterations - z_probability_burnin);
  const int z_batch_cap = 32;
  int z_n_batches = (int) std::floor(std::sqrt((double) expected_z_draws));
  if (z_n_batches < 2) z_n_batches = 2;
  if (z_n_batches > z_batch_cap) z_n_batches = z_batch_cap;
  int z_batch_size = expected_z_draws / z_n_batches;
  if (z_batch_size < 1) z_batch_size = 1;
  // counts[b * 3 * N * N + k * N * N + j * N + i]
  std::vector<double> z_batch_counts((size_t) z_n_batches * 3 * N * N, 0.0);
  std::vector<int> z_batch_draws((size_t) z_n_batches, 0);

  // ==========================================================================
  // MAIN MCMC LOOP
  // ==========================================================================

  // Diagnostic only: capture snapshots of z at a few checkpoint iterations so
  // the trajectory of latent-state assignment can be inspected against known
  // ground truth in simulation studies. Not used by the sampler itself; a
  // purely additive output.
  List z_checkpoints;
  std::vector<int> checkpoint_targets;
  {
    const char *dense_env = std::getenv("HICPOTTS_DENSE_CHECKPOINTS");
    if (dense_env != nullptr) {
      int dense_upto = std::atoi(dense_env);
      for (int t = 1; t <= dense_upto; t++) checkpoint_targets.push_back(t);
    } else {
      checkpoint_targets = {1, 5, 20, 50, 100, 250, 500, 1000, 2000};
    }
  }

  int executed_iterations = iterations;
  bool mcse_converged = false;
  double mcse_max_ratio = R_PosInf;
  std::vector<double> mcse_check_iterations;
  std::vector<double> mcse_check_ratios;
  typedef std::chrono::steady_clock sampler_clock;
  const sampler_clock::time_point sampler_started = sampler_clock::now();
  double seconds_branch = 0.0;
  double seconds_allocation = 0.0;
  double seconds_block_move = 0.0;
  double seconds_empirical_bayes = 0.0;
  double seconds_parameters = 0.0;
  double seconds_gamma = 0.0;
  double seconds_theta = 0.0;
  auto elapsed_seconds = [](sampler_clock::time_point start) {
    return std::chrono::duration<double>(sampler_clock::now() - start).count();
  };

  {
    NumericMatrix c1 = chains[0], c2 = chains[1], c3 = chains[2];
    branch_gap_chain[0] = noise_branch_gap(
      c1(0, _), c2(0, _), c3(0, _), covariate_sds,
      noise_link_sd, noise_order_strength, noise_order_width);
  }

  for (int iter = 0; iter < iterations; iter++) {
   if (verbose && (iter + 1) % progress_interval == 0) {
     Rcout << "Iteration: " << iter + 1 << std::endl;
    }

    const sampler_clock::time_point branch_started = sampler_clock::now();

    // ==== REVERSIBLE WHOLE COMPONENT-2/3 BRANCH TRANSITION =================
    // Propose T(x): exchange every latent 2/3 label together with the complete
    // beta and dispersion states of components 2 and 3. T is an involution, so
    // its proposal ratio is one. The Potts energy is invariant to a global
    // permutation and the component likelihood follows its parameters and
    // allocations, but we nevertheless evaluate the complete component-2/3
    // posterior ratio so component-specific user priors are handled exactly.
    // The labelled component-1/3 prior then favours the biologically intended
    // orientation. This is an in-sampler MH move, never post-hoc relabelling.
    //
    // The ordinary beta/size scan below immediately follows an accepted move,
    // providing the conditional parameter refresh without a second target or
    // a duplicated update implementation. Stored row `iter` is temporarily
    // used as the working state and restored before the iteration ends, so the
    // returned history remains the actual pre-transition sample.
    bool branch_swapped_this_iter = false;
    NumericVector saved_beta2(5), saved_beta3(5);
    double saved_size2 = size_chain(1, iter);
    double saved_size3 = size_chain(2, iter);
    const int branch_start = z_warmup;
    if (branch_swap_enabled && iter >= branch_start &&
        ((iter + 1) % branch_swap_interval == 0)) {
      branch_swap_attempts++;
      if (iter < eb_freeze_iteration) branch_swap_warmup_attempts++;
      NumericMatrix c1 = chains[0], c2 = chains[1], c3 = chains[2];
      NumericVector beta1 = c1(iter, _);
      NumericVector beta2 = c2(iter, _);
      NumericVector beta3 = c3(iter, _);
      for (int j = 0; j < 5; ++j) {
        saved_beta2[j] = beta2[j];
        saved_beta3[j] = beta3[j];
      }

      NumericMatrix z_swapped = clone(z_current);
      swap_signal_labels(z_swapped);
      // Use one adaptive-or-frozen EB object on both sides of this involution.
      // Re-estimating a different prior after swapping labels would put an
      // allocation-dependent hyperparameter change inside the MH ratio.
      const HicpottsNormalPrior prior2 = normal_prior_from_list(
        active_prior_object, 2);
      const HicpottsNormalPrior prior3 = normal_prior_from_list(
        active_prior_object, 3);
      const double current23 =
        component_logposterior_cached(
          beta2, z_current, y, 2, theta[iter], saved_size2, dist,
          logx11, logx22, logx33, logx44, prior2) +
        component_logposterior_cached(
          beta3, z_current, y, 3, theta[iter], saved_size3, dist,
          logx11, logx22, logx33, logx44, prior3);
      const double swapped23 =
        component_logposterior_cached(
          beta3, z_swapped, y, 2, theta[iter], saved_size3, dist,
          logx11, logx22, logx33, logx44, prior2) +
        component_logposterior_cached(
          beta2, z_swapped, y, 3, theta[iter], saved_size2, dist,
          logx11, logx22, logx33, logx44, prior3);
      const double relationship_current = noise_relationship_logprior(
        beta1, beta2, beta3, covariate_sds, noise_link_sd,
        noise_order_strength, noise_order_width);
      const double relationship_swapped = noise_relationship_logprior(
        beta1, beta3, beta2, covariate_sds, noise_link_sd,
        noise_order_strength, noise_order_width);
      const double current_gap = noise_branch_gap(
        beta1, beta2, beta3, covariate_sds, noise_link_sd,
        noise_order_strength, noise_order_width);
      const double proposed_gap = -current_gap;
      const bool toward_intended = proposed_gap > current_gap;
      if (toward_intended) branch_swap_toward_intended_attempts++;
      const double log_alpha_branch =
        (swapped23 - current23) +
        (relationship_swapped - relationship_current);
      if (std::isfinite(log_alpha_branch) &&
          std::log(R::runif(0.0, 1.0)) < log_alpha_branch) {
        z_current = z_swapped;
        for (int j = 0; j < 5; ++j) {
          c2(iter, j) = saved_beta3[j];
          c3(iter, j) = saved_beta2[j];
        }
        size_chain(1, iter) = saved_size3;
        size_chain(2, iter) = saved_size2;
        branch_swapped_this_iter = true;
        branch_swap_accepts++;
        if (iter < eb_freeze_iteration) branch_swap_warmup_accepts++;
        if (toward_intended) branch_swap_toward_intended_accepts++;
      }
    }

    seconds_branch += elapsed_seconds(branch_started);
    const sampler_clock::time_point allocation_started = sampler_clock::now();

    // ==== SEQUENTIAL CHECKERBOARD (RED-BLACK) GIBBS z-SWEEP ===================
    // The three exact component conditionals are now evaluated directly for
    // each cell from cached log1p covariates and current neighbour counts. The
    // previous route generated two random label permutations, constructed six
    // N-by-N temporary matrices per sweep and called back into R repeatedly;
    // it represented the same three-way categorical distribution but with
    // substantially more allocation and interpreter overhead.
    double beta_temp = 1.0;
    if (tempering_warmup > 0 && iter < tempering_warmup) {
      double phase = (double)(iter % tempering_cycle) / (double)tempering_cycle;
      beta_temp = tempering_beta_min + (1.0 - tempering_beta_min) * phase;
    }
    NumericMatrix z_next = clone(z_current);
    if (validation_freeze_z || iter < z_warmup) {
      // Z-WARMUP: leave z_next == z_current (frozen at z_start) so parameters
      // equilibrate before the field is released. No z-sweep this iteration.
    } else {
      checkerboard_z_sweep_cached(
        z_next, y, chains, size_chain, iter, theta[iter], chain_gamma[iter],
        dist, beta_temp, logx11, logx22, logx33, logx44,
        eb_weight1, eb_weight2, eb_weight3);
    }
    seconds_allocation += elapsed_seconds(allocation_started);

    // Non-local 2<->3 cluster move, layered on top of the crossing-barrier
    // penalty; it runs before the beta/size updates read the new z field.
    const sampler_clock::time_point block_move_started = sampler_clock::now();
    if (!validation_freeze_z && !ablation_disable_blockmove &&
        signal_block_move_interval > 0 &&
        ((iter + 1) % signal_block_move_interval == 0)) {
      signal_block_move_attempts++;
      if (iter < eb_freeze_iteration) signal_block_warmup_attempts++;
      bool blockmove_accepted = signal_cluster_swap_23(z_next, y,
                                                       logx11, logx22, logx33, logx44,
                                                       chains, size_chain, N, iter,
                                                       chain_gamma[iter], dist);
      if (blockmove_accepted) signal_block_move_accepts++;
      if (blockmove_accepted && iter < eb_freeze_iteration)
        signal_block_warmup_accepts++;
      const char *dbg_env = std::getenv("HICPOTTS_DEBUG_BLOCKMOVE");
      if (dbg_env != nullptr) {
        Rcout << "[blockmove] iter=" << (iter + 1)
              << " accepted=" << blockmove_accepted << std::endl;
      }
    }
    seconds_block_move += elapsed_seconds(block_move_started);

    // Update empirical-Bayes hyperparameters from soft membership weights only
    // during warm-up. At eb_freeze_iteration they become fixed and every
    // retained beta/size/branch transition uses the same posterior target.
    const sampler_clock::time_point empirical_bayes_started = sampler_clock::now();
    if (use_data_priors && iter < eb_freeze_iteration) {
      active_prior_object = soft_eb_priors_impl(
        y, logx11, logx22, logx33, logx44,
        eb_weight1, eb_weight2, eb_weight3, true);
      eb_prior_updates++;
    }
    seconds_empirical_bayes += elapsed_seconds(empirical_bayes_started);
    HicpottsNormalPrior iteration_priors[3] = {
      normal_prior_from_list(active_prior_object, 1),
      normal_prior_from_list(active_prior_object, 2),
      normal_prior_from_list(active_prior_object, 3)
    };
    const sampler_clock::time_point parameters_started = sampler_clock::now();

    // MCMC proposal update for betas, per component ---------------
    //
    // (Barrier conditioning is read directly from the chains at the point of
    // each component's update -- see the barrier block below. Iteration-start
    // snapshots are deliberately NOT used: component 3 must condition on the
    // component-2 value accepted earlier in this same sweep.)

    for (int comp = 1; comp <= 3; comp++) {
      NumericMatrix chain_matrix = chains[comp - 1];
      NumericVector proposal(5);
      NumericVector sd_component = as<NumericVector>(sd_values[comp - 1]);

      // Current state on the original coefficient scale, then mapped into
      // centred/scaled proposal coordinates.
      arma::vec beta_cur(beta_dimension);
      for (int j = 0; j < beta_dimension; j++)
        beta_cur[j] = chain_matrix(iter, j);
      arma::vec u_cur = beta_to_qr(beta_cur, covariate_means, covariate_sds, R_scaled);

      // Use one QR-whitened random-walk kernel throughout.  The earlier
      // implementation replaced this kernel at iterations / 10 with a
      // full-history adaptive covariance estimated while z and the
      // empirical-Bayes hyperparameters were still changing.  On realistic
      // matrices that covariance generated proposals with zero acceptance.
      // Worse, the rejection-rate adaptation below continued changing
      // sd_values even though the replacement kernel no longer used them, so
      // it could never recover.  Retaining this symmetric kernel preserves the
      // posterior target while keeping the tuned scale live until burn-in.
      arma::vec u_prop(beta_dimension);
      for (int j = 0; j < beta_dimension; j++)
        u_prop[j] = R::rnorm(u_cur[j], sd_component[j]);
      arma::vec beta_prop = qr_to_beta(
        u_prop, covariate_means, covariate_sds, R_scaled);
      for (int j = 0; j < beta_dimension; j++)
        proposal[j] = beta_prop[j];
      if (validation_freeze_beta) {
        for (int j = 0; j < beta_dimension; ++j)
          proposal[j] = chain_matrix(iter, j);
      } else {
        beta_attempt_total[comp - 1]++;
        if (iter < eb_freeze_iteration)
          beta_warmup_attempt_total[comp - 1]++;
      }

      double current_size = size_chain(comp - 1, iter);
      double posterior_current = component_logposterior_cached(
        chain_matrix(iter, _), z_next, y, comp, theta[iter], current_size,
        dist, logx11, logx22, logx33, logx44,
        iteration_priors[comp - 1]);
      double posterior_proposal = component_logposterior_cached(
        proposal, z_next, y, comp, theta[iter], current_size,
        dist, logx11, logx22, logx33, logx44,
        iteration_priors[comp - 1]);

      // Identifiability / label-switching handling for components 2/3:
      // no hard ordering constraint is imposed (order_ok is unconditionally
      // true); instead the symmetric crossing-barrier penalty (see
      // crossing_barrier()'s definition near the top of this file) is added to
      // component 2's and component 3's own acceptance ratio.
      //
      // Each component conditions on the OTHER component's CURRENT state at
      // the moment of its own update, which is what makes this a valid
      // Metropolis-within-Gibbs scan: every kernel targets the correct full
      // conditional given the rest of the state.
      //
      // Components are processed in order 2 then 3, so:
      //   comp 2 - component 3 has NOT yet moved this sweep, its current value
      //            is still row `iter`.
      //   comp 3 - component 2 HAS already moved, its current value is row
      //            `iter + 1`.
      // Previously BOTH read iteration-start snapshots, so component 3
      // conditioned on a stale component-2 intercept. That is neither a valid
      // sequential conditional nor a jointly accepted block proposal, and the
      // barrier term it produced could differ from the correct one by orders
      // of magnitude. The snapshots were introduced to keep the penalty
      // independent of processing order, but order dependence in a sequential
      // scan is standard and valid - staleness is not.
      //
      // The barrier's functional form, kappa and w are unchanged, so the
      // model and its implicit repulsive prior on a3 - a2 are exactly as
      // before; only the conditioning is corrected.
      bool order_ok = true;
      double barrier_diff = 0.0;
      if (!ablation_disable_barrier && comp == 2) {
        double a2_curr = chain_matrix(iter, 0);
        double a3_curr = as<NumericMatrix>(chains[2])(iter, 0);      // not yet updated
        barrier_diff = crossing_barrier(a3_curr - proposal[0], comp23_barrier_kappa, comp23_barrier_w)
                     - crossing_barrier(a3_curr - a2_curr,     comp23_barrier_kappa, comp23_barrier_w);
      } else if (!ablation_disable_barrier && comp == 3) {
        double a3_curr = chain_matrix(iter, 0);
        double a2_curr = as<NumericMatrix>(chains[1])(iter + 1, 0);  // already updated
        barrier_diff = crossing_barrier(proposal[0] - a2_curr, comp23_barrier_kappa, comp23_barrier_w)
                     - crossing_barrier(a3_curr     - a2_curr, comp23_barrier_kappa, comp23_barrier_w);
      }

      // Coupled biological prior. Build the live sequential state: components
      // already visited in this sweep use row iter+1 and later components use
      // row iter. Replacing only the component currently proposed yields the
      // exact joint-prior ratio for this Metropolis-within-Gibbs update.
      double relationship_diff = 0.0;
      if (use_noise_relationship_prior) {
        NumericVector live[3] = {
          NumericVector(5), NumericVector(5), NumericVector(5)
        };
        NumericVector proposed_live[3] = {
          NumericVector(5), NumericVector(5), NumericVector(5)
        };
        for (int c = 0; c < 3; c++) {
          NumericMatrix cm = chains[c];
          const int row = c < (comp - 1) ? iter + 1 : iter;
          for (int j = 0; j < 5; j++) {
            live[c][j] = cm(row, j);
            proposed_live[c][j] = cm(row, j);
          }
        }
        for (int j = 0; j < 5; j++) proposed_live[comp - 1][j] = proposal[j];
        const double prior_current = noise_relationship_logprior(
          live[0], live[1], live[2], covariate_sds, noise_link_sd,
          noise_order_strength, noise_order_width);
        const double prior_proposal = noise_relationship_logprior(
          proposed_live[0], proposed_live[1], proposed_live[2], covariate_sds,
          noise_link_sd, noise_order_strength, noise_order_width);
        relationship_diff = prior_proposal - prior_current;
      }

      // Call proposaldensity_combined to get log proposal densities for current and proposed states
       // Symmetric random-walk proposal: the proposal-density ratio cancels.
      double log_alpha = (posterior_proposal - posterior_current) +
                         barrier_diff + relationship_diff;
      if (order_ok && std::log(R::runif(0, 1)) < log_alpha) {
        for (int k = 0; k < 5; k++)
          chain_matrix(iter + 1, k) = proposal[k];
        acceptance_counts[comp - 1]++;
        if (!validation_freeze_beta) {
          beta_accept_total[comp - 1]++;
          if (iter < eb_freeze_iteration)
            beta_warmup_accept_total[comp - 1]++;
        }
      } else {
        for (int k = 0; k < 5; k++)
          chain_matrix(iter + 1, k) = chain_matrix(iter, k);
      }

      // Update size for NB or ZINB with a log-scale adaptive random walk.
      // Propose log(size') = log(size) + N(0, step)  =>  size' = size * exp(step).
      // This is always positive and scale-invariant. Because the target is on the
      // size scale, the change of variables contributes a Jacobian term
      // (+ log(size') - log(size)) to the Metropolis acceptance ratio.
      if ((dist == "NB" || dist == "ZINB") && validation_freeze_size) {
        size_chain(comp - 1, iter + 1) = current_size;
      } else if (dist == "NB" || dist == "ZINB") {
        size_attempt_total[comp - 1]++;
        double log_step      = R::rnorm(0.0, size_logstep[comp - 1]);
        double size_proposal = current_size * std::exp(log_step);
        if (size_proposal > 0.0 && std::isfinite(size_proposal)) {
          double posterior_current_size = component_logposterior_cached(
            chain_matrix(iter + 1, _), z_next, y, comp, theta[iter],
            current_size, dist, logx11, logx22, logx33, logx44,
            iteration_priors[comp - 1]);
          double posterior_proposal_size = component_logposterior_cached(
            chain_matrix(iter + 1, _), z_next, y, comp, theta[iter],
            size_proposal, dist, logx11, logx22, logx33, logx44,
            iteration_priors[comp - 1]);

          double log_alpha_s = posterior_proposal_size - posterior_current_size
                             + std::log(size_proposal) - std::log(current_size);  // Jacobian
          if (std::log(R::runif(0, 1)) < log_alpha_s) {
            size_chain(comp - 1, iter + 1) = size_proposal;
            size_accept[comp - 1]++;
            size_accept_total[comp - 1]++;
          } else {
            size_chain(comp - 1, iter + 1) = current_size;
          }
        } else {
          size_chain(comp - 1, iter + 1) = current_size;
        }
      }
    }

    // Adaptive tuning of sd_values
    if (((iter+1) % adaptation_interval == 0) &&
        ((iter + 1) <= eb_freeze_iteration)) {
      for (int comp = 1; comp <= 3; comp++) {
        double acceptance_rate = static_cast<double>(acceptance_counts[comp - 1]) / adaptation_interval;
        NumericVector sd_component = as<NumericVector>(sd_values[comp - 1]);
        if (acceptance_rate < accept_lo) {
          sd_component = sd_component / adaptation_scaling;   // too few accepts -> smaller steps
        } else if (acceptance_rate > accept_hi) {
          sd_component = sd_component * adaptation_scaling;   // too many accepts -> larger steps
        }
        sd_values[comp - 1] = sd_component;
        acceptance_counts[comp - 1] = 0;

        // Adapt the log-size step the same way (NB / ZINB only).
        if (dist == "NB" || dist == "ZINB") {
          double s_rate = static_cast<double>(size_accept[comp - 1]) / adaptation_interval;
          if (s_rate < size_accept_lo)
            size_logstep[comp - 1] /= adaptation_scaling;
          else if (s_rate > size_accept_hi)
            size_logstep[comp - 1] *= adaptation_scaling;
          size_accept[comp - 1] = 0;
        }
      }

      // Adapt gamma's logit-scale step the same way.
      // Denominator is the number of gamma proposals actually made in this
      // window, not adaptation_interval: with gamma_update_interval > 1 the
      // gamma step fires only every k-th iteration, so dividing by the full
      // interval would understate the acceptance rate and wrongly shrink the
      // step.
      if (gamma_local_proposals_window > 0) {
        double g_rate = static_cast<double>(gamma_local_accept_window) /
                        (double) gamma_local_proposals_window;
        if (g_rate < accept_lo)      gamma_logit_step /= adaptation_scaling;
        else if (g_rate > accept_hi) gamma_logit_step *= adaptation_scaling;
      }
      gamma_local_accept_window = 0;
      gamma_local_proposals_window = 0;

    }

    seconds_parameters += elapsed_seconds(parameters_started);
    const sampler_clock::time_point gamma_started = sampler_clock::now();

    // ==== RETAINED-STATE ABC UPDATE FOR THE POTTS COUPLING gamma =============
    // gamma enters the model only through p(z | gamma) = exp(gamma S(z))/Z(gamma),
    // where S(z) is the raw number of agreeing unordered edges. S(z) is the
    // sufficient statistic for gamma and Z(gamma) is intractable. A Gaussian
    // ABC kernel supplies the likelihood-free conditional update.
    //
    // Conditioning value: S(z_next), the sufficient statistic of the sampler's
    // own current latent field. This makes the step an update of the CONDITIONAL
    // gamma | z, the correct conditional in the surrounding Gibbs scheme, so it
    // is recomputed every sweep.
    //
    // Simulator: a Potts draw at the candidate gamma. It is a function of gamma
    // and N alone, so it cannot drift as beta/size/theta move -- the retained
    // auxiliary draw is always a draw from exactly the distribution used in the
    // acceptance ratio. Only its DISTANCE is re-scored against the moving S(z);
    // the draw itself is never re-simulated. See src/abc_s1s2_gamma.h.
    if (validation_freeze_gamma) {
      T_obs_fixed = potts_label_agreement(z_next);
      chain_gamma[iter + 1] = chain_gamma[iter];
      abc_distance_chain[iter + 1] = NA_REAL;
      abc_s1_chain[iter + 1] = T_obs_fixed;
      abc_s2_chain[iter + 1] = NA_REAL;
    } else {
    T_obs_fixed = potts_label_agreement(z_next);
    if (!abc_calibrated) {
      // A supplied epsilon defines the target directly, so no calibration
      // simulations are required. Repeated-chain wrappers exploit this by
      // calibrating once and passing the shared value to every worker.
      if (std::isfinite(abc_user_epsilon)) {
        abc_epsilon = abc_user_epsilon;
        abc_scale = NA_REAL;
      } else {
        // A standalone low-level fit retains the historical deterministic
        // calibration: run one prior-predictive pass on an isolated fixed RNG
        // substream, then restore the caller's chain stream exactly.
        Environment base_env = Environment::base_env();
        Environment global_env = Environment::global_env();
        Function set_seed = base_env["set.seed"];

        PutRNGstate();
        RObject saved_seed = R_NilValue;
        const bool had_seed = global_env.exists(".Random.seed");
        if (had_seed) saved_seed = global_env[".Random.seed"];

        set_seed(abc_calibration_seed);
        GetRNGstate();
        const AbcToleranceCalibration calibration = calibrate_abc_tolerance(
          abc_potts_workspace, gamma_prior_shape1, gamma_prior_shape2,
          abc_epsilon_quantile, abc_calibration_reps, abc_sim_reps,
          abc_potts_sweeps);
        PutRNGstate();
        if (had_seed) global_env.assign(".Random.seed", saved_seed);
        GetRNGstate();

        abc_epsilon = calibration.epsilon;
        abc_scale = calibration.scale;
        abc_calibration_reps_executed = abc_calibration_reps;
      }

      // Initialise the augmented state: draw the auxiliary summary paired with
      // the CURRENT gamma once, here. From now on it is carried forward and
      // only replaced on acceptance.
      {
        double acc = 0.0;
        for (int q = 0; q < abc_sim_reps; q++)
          acc += simulate_potts_abc_summary(
            abc_potts_workspace, chain_gamma[iter], abc_potts_sweeps);
        abc_T_state = acc / abc_sim_reps;
        abc_d_state = std::abs(abc_T_state - T_obs_fixed);
      }
      abc_calibrated = true;
    }

    // Thinned gamma update: the ABC step (the per-iteration cost driver) fires
    // only every gamma_update_interval-th iteration. On the other iterations
    // gamma and its auxiliary state are carried forward unchanged, so the
    // cheap, well-mixing beta/size/theta/z updates still run every iteration
    // while the expensive Potts simulations run 1/interval as often. iter == 0
    // always satisfies the test, so calibration above is never skipped.
    if ((iter % gamma_update_interval) == 0) {
      double gamma_cur = chain_gamma[iter];
      double logit_cur = std::log(gamma_cur / (1.0 - gamma_cur));
      // Proposal type: 0 = adaptive local logit RW, 1 = wide logit RW,
      // 2 = independent draw from the configured Beta prior. The complete
      // mixture density, rather than only the selected component density, is
      // evaluated below in both directions.
      (void) logit_cur;
      int gamma_proposal_type = 0;
      double gamma_prop = draw_gamma_mixture(
        gamma_cur, gamma_logit_step,
        gamma_large_jump_probability, gamma_large_jump_multiplier,
        gamma_independence_probability,
        gamma_prior_shape1, gamma_prior_shape2, gamma_proposal_type);
      double gamma_next = gamma_cur;

      // Simulated Potts sufficient statistic at a given gamma (averaged over reps).
      // Draws from p(z | gamma), which depends on gamma and N alone, so the
      // retained auxiliary state and this proposal come from exactly the same
      // simulator family that the acceptance ratio assumes.
      auto sim_summary = [&](double g) {
        double acc = 0.0;
        for (int q = 0; q < abc_sim_reps; q++)
          acc += simulate_potts_abc_summary(
            abc_potts_workspace, g, abc_potts_sweeps);
        return acc / abc_sim_reps;
      };

      // ---- augmented-state ABC-MCMC -------------------------------------
      // State is the PAIR (gamma, x) where x is the auxiliary dataset,
      // summarised by T. The target is
      //     pi(gamma, x | y)  ∝  p(gamma) p(x | gamma) K_eps(|T(x) - T(y)|)
      // whose gamma-marginal is the ABC posterior. A move proposes gamma' and
      // draws x' ~ p(. | gamma') afresh; on ACCEPTANCE the pair becomes
      // (gamma', x'), on REJECTION the pair (gamma, x) is retained IN FULL -
      // including the auxiliary draw.
      //
      // Previously the current gamma's auxiliary summary was re-simulated at
      // every iteration and discarded on rejection. That is
      // Monte-Carlo-within-Metropolis (noisy MCMC), not ABC-MCMC: the
      // acceptance ratio compared two fresh draws rather than a retained state
      // against a proposal, so the chain did not target the stated ABC
      // posterior. Averaging abc_sim_reps summaries before applying the
      // Gaussian kernel is also not a non-negative unbiased likelihood
      // estimator, so the pseudo-marginal route was not available either.
      // The retained auxiliary DRAW is kept (never re-simulated); only its
      // distance is re-scored against the current field's sufficient
      // statistic, because the conditioning value S(z) moves as z moves. This
      // is what keeps the step a valid update of the conditional gamma | z.
      double d_cur = std::abs(abc_T_state - T_obs_fixed);
      abc_d_state = d_cur;
      double sel_d = d_cur;
      double sel_T = abc_T_state;

      if (gamma_prop > 0.0 && gamma_prop < 1.0 && std::isfinite(gamma_prop)) {
        gamma_proposals_total++;
        if (gamma_proposal_type == 0) {
          gamma_local_proposals_total++;
          gamma_local_proposals_window++;
        } else if (gamma_proposal_type == 1) {
          gamma_large_proposals_total++;
        } else {
          gamma_independence_proposals_total++;
        }
        double T_prop = sim_summary(gamma_prop);     // only the PROPOSAL is drawn
        double d_prop = std::abs(T_prop - T_obs_fixed);
        double log_kernel_cur  = -0.5 * std::pow(d_cur  / abc_epsilon, 2.0);
        double log_kernel_prop = -0.5 * std::pow(d_prop / abc_epsilon, 2.0);
        // Configurable Beta prior (default flat Beta(1,1)). The previous
        // hardcoded Beta(2,2) had mean 0.5, above any realistic Potts coupling,
        // so with a loose tolerance it pulled gamma upward.
        double log_prior_ratio =
          R::dbeta(gamma_prop, gamma_prior_shape1, gamma_prior_shape2, true) -
          R::dbeta(gamma_cur,  gamma_prior_shape1, gamma_prior_shape2, true);
        const double log_proposal_ratio =
          gamma_mixture_logdensity(
            gamma_cur, gamma_prop, gamma_logit_step,
            gamma_large_jump_probability, gamma_large_jump_multiplier,
            gamma_independence_probability,
            gamma_prior_shape1, gamma_prior_shape2) -
          gamma_mixture_logdensity(
            gamma_prop, gamma_cur, gamma_logit_step,
            gamma_large_jump_probability, gamma_large_jump_multiplier,
            gamma_independence_probability,
            gamma_prior_shape1, gamma_prior_shape2);
        double log_alpha = (log_kernel_prop - log_kernel_cur) +
                           log_prior_ratio + log_proposal_ratio;
        if (std::log(R::runif(0.0, 1.0)) < log_alpha) {
          gamma_next = gamma_prop;
          sel_d = d_prop;
          sel_T = T_prop;
          gamma_accept_total++;
          if (gamma_proposal_type == 0) {
            gamma_local_accept_total++;
            gamma_local_accept_window++;
          } else if (gamma_proposal_type == 1) {
            gamma_large_accept_total++;
          } else {
            gamma_independence_accept_total++;
          }
        }
      }
      // Carry the accepted (or retained) auxiliary forward as the new state.
      abc_T_state = sel_T;
      abc_d_state = sel_d;
      chain_gamma[iter + 1] = gamma_next;
      // Both diagnostics now describe the SAME retained auxiliary state.
      // Previously T_sim held a freshly simulated current-gamma summary while
      // abc_distance held the accepted proposal's distance, so the two could
      // describe different auxiliary draws.
      abc_distance_chain[iter + 1] = sel_d;
      abc_s1_chain[iter + 1] = T_obs_fixed;
      abc_s2_chain[iter + 1] = sel_T;
    } else {
      // Off-cycle: gamma and the auxiliary draw are retained. The conditioning
      // field z may nevertheless have moved, so re-score the retained draw
      // against the current T_obs rather than carrying a stale distance.
      chain_gamma[iter + 1] = chain_gamma[iter];
      abc_d_state = std::abs(abc_T_state - T_obs_fixed);
      abc_distance_chain[iter + 1] = abc_d_state;
      abc_s1_chain[iter + 1] = T_obs_fixed;
      abc_s2_chain[iter + 1] = abc_T_state;
    }
    }

    seconds_gamma += elapsed_seconds(gamma_started);
    const sampler_clock::time_point theta_started = sampler_clock::now();

     // Update theta with the Beta(1,1)-Bernoulli conjugate Gibbs step.
    if ((dist == "ZIP" || dist == "ZINB") && validation_freeze_theta) {
      theta[iter + 1] = theta[iter];
    } else if (dist == "ZIP" || dist == "ZINB") {
      theta[iter + 1] = update_theta_gibbs(
        z_next, y, chains, iter + 1, size_chain,
        logx11, logx22, logx33, logx44, theta[iter], dist);
    } else {
      theta[iter + 1] = theta[iter];
    }
    seconds_theta += elapsed_seconds(theta_started);

    {
      NumericMatrix c1 = chains[0], c2 = chains[1], c3 = chains[2];
      branch_gap_chain[iter + 1] = noise_branch_gap(
        c1(iter + 1, _), c2(iter + 1, _), c3(iter + 1, _),
        covariate_sds, noise_link_sd, noise_order_strength,
        noise_order_width);
      if (branch_swapped_this_iter) {
        for (int j = 0; j < 5; ++j) {
          c2(iter, j) = saved_beta2[j];
          c3(iter, j) = saved_beta3[j];
        }
        size_chain(1, iter) = saved_size2;
        size_chain(2, iter) = saved_size3;
      }
    }
     z_current = z_next;

    const int completed_for_z = iter + 1;
    if (completed_for_z > z_probability_burnin) {
      // Contiguous batch index; the final batch absorbs any remainder so no
      // draw is discarded.
      int b = z_probability_draws / z_batch_size;
      if (b >= z_n_batches) b = z_n_batches - 1;
      const size_t batch_base = (size_t) b * 3 * N * N;
      for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
          const int label = (int)z_current(i, j);
          if (label == 1) z_probability_count1(i, j) += 1.0;
          else if (label == 2) z_probability_count2(i, j) += 1.0;
          else if (label == 3) z_probability_count3(i, j) += 1.0;
          if (label >= 1 && label <= 3)
            z_batch_counts[batch_base + (size_t)(label - 1) * N * N +
                           (size_t) j * N + i] += 1.0;
        }
      }
      z_batch_draws[(size_t) b]++;
      z_probability_draws++;
    }

    for (size_t ck = 0; ck < checkpoint_targets.size(); ck++) {
      if (iter + 1 == checkpoint_targets[ck]) {
        NumericMatrix z_reported_checkpoint = clone(z_current);
        z_checkpoints.push_back(z_reported_checkpoint, "iter_" + std::to_string(iter + 1));
        break;
      }
    }

    const int completed = iter + 1;
    executed_iterations = completed;
    // Early stopping must not cut the run short of the latent-state
    // accumulation window. The parameter chains can satisfy their MCSE
    // criterion well before the z-probability burn-in elapses; stopping there
    // leaves ZERO post-burn-in allocation draws, so the fit converges on its
    // parameters yet cannot be classified at all. Requiring a minimum number
    // of retained allocation draws before honouring the stop keeps the two
    // criteria consistent.
    const int min_alloc_draws_before_stop = 100;
    if (mcse_stop && completed >= mcse_min_iterations &&
        completed > tempering_warmup + 200 &&
        z_probability_draws >= min_alloc_draws_before_stop &&
        completed % mcse_check_interval == 0) {
      int start = std::max(tempering_warmup, completed / 2);
      double worst = 0.0;
      bool valid = true;
      for (int comp = 0; comp < 3; comp++) {
        NumericMatrix cm = chains[comp];
        for (int col = 0; col < 5; col++) {
          std::vector<double> values;
          values.reserve(completed - start + 1);
          for (int row = start; row <= completed; row++)
            values.push_back(cm(row, col));
          double ratio = relative_mcse(values);
          if (!std::isfinite(ratio)) valid = false;
          else worst = std::max(worst, ratio);
        }
        if (dist == "NB" || dist == "ZINB") {
          std::vector<double> values;
          values.reserve(completed - start + 1);
          for (int row = start; row <= completed; row++)
            values.push_back(size_chain(comp, row));
          double ratio = relative_mcse(values);
          if (!std::isfinite(ratio)) valid = false;
          else worst = std::max(worst, ratio);
        }
      }
      std::vector<double> gamma_values;
      gamma_values.reserve(completed - start + 1);
      for (int row = start; row <= completed; row++)
        gamma_values.push_back(chain_gamma[row]);
      double gamma_ratio = relative_mcse(gamma_values);
      if (!std::isfinite(gamma_ratio)) valid = false;
      else worst = std::max(worst, gamma_ratio);
      if (dist == "ZIP" || dist == "ZINB") {
        std::vector<double> theta_values;
        theta_values.reserve(completed - start + 1);
        for (int row = start; row <= completed; row++) theta_values.push_back(theta[row]);
        double theta_ratio = relative_mcse(theta_values);
        if (!std::isfinite(theta_ratio)) valid = false;
        else worst = std::max(worst, theta_ratio);
      }
      mcse_max_ratio = valid ? worst : R_PosInf;
      mcse_check_iterations.push_back((double)completed);
      mcse_check_ratios.push_back(mcse_max_ratio);
      if (valid && worst <= mcse_relative_threshold) {
        mcse_converged = true;
        break;
      }
    }

  }

  const double seconds_total = elapsed_seconds(sampler_started);
  const double seconds_profiled = seconds_branch + seconds_allocation +
    seconds_block_move + seconds_empirical_bayes + seconds_parameters +
    seconds_gamma + seconds_theta;

  List chains_out(3);
  for (int comp = 0; comp < 3; comp++) {
    NumericMatrix cm = chains[comp];
    chains_out[comp] = cm(Range(0, executed_iterations), _);
  }
  NumericVector gamma_out = chain_gamma[Range(0, executed_iterations)];
  NumericVector theta_out = theta[Range(0, executed_iterations)];
  NumericMatrix size_out = size_chain(_, Range(0, executed_iterations));
  NumericVector abc_distance_out = abc_distance_chain[Range(0, executed_iterations)];
  NumericVector abc_s1_out = abc_s1_chain[Range(0, executed_iterations)];
  NumericVector abc_s2_out = abc_s2_chain[Range(0, executed_iterations)];
  NumericVector branch_gap_out = branch_gap_chain[Range(0, executed_iterations)];

  gamma_out.attr("gamma_update") = "ABC_potts_label_agreement";
  gamma_out.attr("gamma_method") = gamma_method;
  gamma_out.attr("abc_distance") = abc_distance_out;    // |T_sim - T_obs| per iter
  gamma_out.attr("T_obs") = abc_s1_out;                 // current-label adjacency excess
  gamma_out.attr("T_sim") = abc_s2_out;                 // auxiliary-label adjacency excess
  gamma_out.attr("abc_epsilon") = abc_epsilon;
  gamma_out.attr("abc_scale") = abc_scale;
  // Record both the selected gamma method and its auxiliary-simulation
  // configuration so fits remain auditable and sensitivity runs comparable.
  gamma_out.attr("abc_epsilon_quantile") = abc_epsilon_quantile;
  gamma_out.attr("abc_epsilon_user_supplied") = std::isfinite(abc_user_epsilon);
  gamma_out.attr("abc_epsilon_shared_calibration") = false;
  gamma_out.attr("abc_calibration_reps_executed") =
    abc_calibration_reps_executed;
  gamma_out.attr("abc_simulator_engine") = "compact_lookup_v1";
  gamma_out.attr("gamma_prior_shape1") = gamma_prior_shape1;
  gamma_out.attr("gamma_prior_shape2") = gamma_prior_shape2;
  gamma_out.attr("gamma_prior_mean") =
    gamma_prior_shape1 / (gamma_prior_shape1 + gamma_prior_shape2);
  // Auxiliary-Potts cost/thinning settings actually used.
  gamma_out.attr("gamma_update_interval") = gamma_update_interval;
  gamma_out.attr("abc_potts_sweeps") = abc_potts_sweeps;
  gamma_out.attr("abc_sim_reps") = abc_sim_reps;
  gamma_out.attr("gamma_acceptance_rate") = gamma_proposals_total > 0 ?
    (double)gamma_accept_total / (double)gamma_proposals_total : NA_REAL;
  gamma_out.attr("gamma_proposals") = (double)gamma_proposals_total;
  gamma_out.attr("gamma_local_proposals") = (double)gamma_local_proposals_total;
  gamma_out.attr("gamma_large_proposals") = (double)gamma_large_proposals_total;
  gamma_out.attr("gamma_independence_proposals") =
    (double)gamma_independence_proposals_total;
  gamma_out.attr("gamma_local_acceptance_rate") =
    gamma_local_proposals_total > 0 ?
      (double)gamma_local_accept_total / (double)gamma_local_proposals_total : NA_REAL;
  gamma_out.attr("gamma_large_acceptance_rate") =
    gamma_large_proposals_total > 0 ?
      (double)gamma_large_accept_total / (double)gamma_large_proposals_total : NA_REAL;
  gamma_out.attr("gamma_independence_acceptance_rate") =
    gamma_independence_proposals_total > 0 ?
      (double)gamma_independence_accept_total /
        (double)gamma_independence_proposals_total : NA_REAL;
  gamma_out.attr("gamma_final_local_logit_step") = gamma_logit_step;
  gamma_out.attr("gamma_large_jump_probability") =
    gamma_large_jump_probability;
  gamma_out.attr("gamma_large_jump_multiplier") = gamma_large_jump_multiplier;
  gamma_out.attr("gamma_independence_probability") =
    gamma_independence_probability;
  gamma_out.attr("iterations_requested") = iterations;
  gamma_out.attr("iterations_completed") = executed_iterations;
  gamma_out.attr("mcse_stopping_enabled") = mcse_stop;
  gamma_out.attr("mcse_converged") = mcse_converged;
  gamma_out.attr("mcse_relative_threshold") = mcse_relative_threshold;
  gamma_out.attr("mcse_max_relative") = mcse_max_ratio;
  gamma_out.attr("mcse_check_iterations") = wrap(mcse_check_iterations);
  gamma_out.attr("mcse_check_max_relative") = wrap(mcse_check_ratios);
  gamma_out.attr("tempering_warmup") = tempering_warmup;
  gamma_out.attr("tempering_beta_min") = tempering_beta_min;
  theta_out.attr("theta_update") = "Beta(1,1)-Bernoulli conjugate Gibbs";
  theta_out.attr("theta_prior_shape1") = 1.0;
  theta_out.attr("theta_prior_shape2") = 1.0;
  // H7: the PRIMARY returned allocation is now the unforced posterior
  // allocation z_current -- the labels the sampler actually drew.
  //
  // Previously the primary object forced every observed zero to component 1.
  // That is not a posterior classification: zeros have positive sampling
  // probability under all three components, and under plain Poisson/NB there
  // is no structural-zero mechanism at all, so the override applied even where
  // the model contains nothing to justify it. It biased reported labels toward
  // "noise" and could contradict both the fitted latent field and the
  // probability output.
  //
  // The forced mask is still available, under an explicit name
  // ("zero_forced_classification"), for anyone who wants the domain
  // convention "every observed zero is noise". It is a mask, not a posterior,
  // and is documented as such.
  NumericMatrix z_zero_forced = clone(z_current);
  IntegerVector parameter_counts(3, 0), classification_counts(3, 0);
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      int internal_label = (int)z_current(i, j);
      if (internal_label >= 1 && internal_label <= 3)
        parameter_counts[internal_label - 1]++;
      if (y(i, j) == 0.0) z_zero_forced(i, j) = 1.0;
      int forced_label = (int)z_zero_forced(i, j);
      if (forced_label >= 1 && forced_label <= 3)
        classification_counts[forced_label - 1]++;
    }
  }
  NumericMatrix z_reported = clone(z_current);   // unforced: the posterior draw
  parameter_counts.attr("names") = CharacterVector::create("component1", "component2", "component3");
  classification_counts.attr("names") = CharacterVector::create("component1", "component2", "component3");
  z_reported.attr("parameter_component_counts") = parameter_counts;
  z_reported.attr("classification_component_counts") = classification_counts;
  z_reported.attr("zero_allocation") = "z_final is the UNFORCED posterior allocation drawn by the sampler; the attribute 'zero_forced_classification' holds the alternative mask in which every observed zero is set to component 1 (a domain convention, not a posterior classification)";
  // The forced mask, explicitly named and explicitly not a posterior.
  z_reported.attr("zero_forced_classification") = z_zero_forced;
  // Retained for backwards compatibility: callers that read z_parameter_final
  // to recover the unforced field still work, and now get the same matrix as
  // the primary return value.
  z_reported.attr("z_parameter_final") = z_current;
  z_reported.attr("minimum_recommended_parameter_cells") = 100;
  // Per-component occupancy is reported through the
  // 'parameter_component_counts' attribute and is surfaced by the occupancy
  // gate in diagnose_hicpotts_fit(); the sampler itself does not warn.

  NumericVector proposal_means = wrap(covariate_means);
  NumericVector proposal_sds = wrap(covariate_sds);
  proposal_means.attr("names") = CharacterVector::create("distance", "GC", "TES", "ACC");
  proposal_sds.attr("names") = CharacterVector::create("distance", "GC", "TES", "ACC");
  for (int comp = 0; comp < 3; comp++) {
    NumericMatrix chain_matrix = chains_out[comp];
    chain_matrix.attr("coefficient_scale") = "original log1p-covariate scale";
    chain_matrix.attr("proposal_preconditioning") =
      qr_whitened ? "QR-whitened centred and standardised log1p covariates" :
      "centred and standardised log1p covariates (QR rank-deficiency fallback)";
    chain_matrix.attr("proposal_covariate_means") = proposal_means;
    chain_matrix.attr("proposal_covariate_sds") = proposal_sds;
    chain_matrix.attr("proposal_qr_R_scaled") = wrap(R_scaled);
    chains_out[comp] = chain_matrix;
  }

  // A run that stops at or before the automatic burn-in still returns a valid
  // (one-draw) summary. This matters for smoke tests and very short exploratory
  // runs; production classification should use converged multi-chain fits.
  // A run that accumulates no post-burn-in latent-state draws cannot be
  // classified. The previous behaviour substituted the single final allocation
  // and reported it as one "draw", which yields degenerate 0/1 probabilities
  // that look like a posterior summary but carry no uncertainty at all: the
  // margin and entropy diagnostics computed from them are meaningless, and
  // downstream code could not distinguish this from a genuine result. It is
  // now an error, so a too-short run fails loudly instead of silently
  // producing an over-confident classification.
  if (z_probability_draws == 0) {
    stop("No post-burn-in latent-state draws were accumulated (iterations = %d, "
         "z-probability burn-in = %d). Run more iterations or lower "
         "z_probability_burnin_arg; classification from the single final "
         "latent state is not supported.",
         iterations, z_probability_burnin);
  }
  const double z_draw_denominator = (double)z_probability_draws;
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      z_probability_count1(i, j) /= z_draw_denominator;
      z_probability_count2(i, j) /= z_draw_denominator;
      z_probability_count3(i, j) /= z_draw_denominator;
    }
  }
  List z_probabilities = List::create(
    Named("component1") = z_probability_count1,
    Named("component2") = z_probability_count2,
    Named("component3") = z_probability_count3);
  z_probabilities.attr("burnin") = z_probability_burnin;
  z_probabilities.attr("draws") = z_probability_draws;
  // Batch MEANS as an N x N x 3 x n_used array, plus the draw count backing
  // each batch. Empty batches (possible only on very short runs) are dropped
  // so downstream variance calculations never see a zero-denominator batch.
  int z_batches_used = 0;
  for (int b = 0; b < z_n_batches; b++) if (z_batch_draws[(size_t) b] > 0) z_batches_used++;
  NumericVector z_batch_means((size_t) N * N * 3 * std::max(1, z_batches_used));
  IntegerVector z_batch_sizes(std::max(1, z_batches_used));
  {
    int out_b = 0;
    for (int b = 0; b < z_n_batches; b++) {
      const int nb = z_batch_draws[(size_t) b];
      if (nb <= 0) continue;
      const size_t in_base  = (size_t) b     * 3 * N * N;
      const size_t out_base = (size_t) out_b * 3 * N * N;
      for (size_t t = 0; t < (size_t) 3 * N * N; t++)
        z_batch_means[out_base + t] = z_batch_counts[in_base + t] / (double) nb;
      z_batch_sizes[out_b] = nb;
      out_b++;
    }
  }
  z_batch_means.attr("dim") =
    IntegerVector::create(N, N, 3, std::max(1, z_batches_used));
  z_probabilities.attr("source") =
    "post-burn-in unrestricted latent-state draws from the fitted MCMC";

  NumericVector size_final_logstep = wrap(size_logstep);
  NumericVector size_attempts = wrap(size_attempt_total);
  NumericVector size_acceptances = wrap(size_accept_total);
  NumericVector size_acceptance_rate(3, NA_REAL);
  for (int comp = 0; comp < 3; ++comp) {
    if (size_attempt_total[comp] > 0)
      size_acceptance_rate[comp] =
        (double) size_accept_total[comp] / (double) size_attempt_total[comp];
  }
  CharacterVector component_labels = CharacterVector::create(
    "component1", "component2", "component3");
  size_final_logstep.attr("names") = component_labels;
  size_attempts.attr("names") = component_labels;
  size_acceptances.attr("names") = component_labels;
  size_acceptance_rate.attr("names") = component_labels;

  NumericVector beta_attempts = wrap(beta_attempt_total);
  NumericVector beta_acceptances = wrap(beta_accept_total);
  NumericVector beta_warmup_attempts = wrap(beta_warmup_attempt_total);
  NumericVector beta_warmup_acceptances = wrap(beta_warmup_accept_total);
  NumericVector beta_retained_attempts(3), beta_retained_acceptances(3);
  NumericVector beta_acceptance_rate(3, NA_REAL);
  NumericVector beta_retained_acceptance_rate(3, NA_REAL);
  NumericMatrix beta_final_step(3, beta_dimension);
  for (int comp = 0; comp < 3; ++comp) {
    beta_retained_attempts[comp] =
      (double)(beta_attempt_total[comp] - beta_warmup_attempt_total[comp]);
    beta_retained_acceptances[comp] =
      (double)(beta_accept_total[comp] - beta_warmup_accept_total[comp]);
    if (beta_attempt_total[comp] > 0)
      beta_acceptance_rate[comp] =
        (double)beta_accept_total[comp] / (double)beta_attempt_total[comp];
    if (beta_retained_attempts[comp] > 0)
      beta_retained_acceptance_rate[comp] =
        beta_retained_acceptances[comp] / beta_retained_attempts[comp];
    NumericVector step = as<NumericVector>(sd_values[comp]);
    for (int coefficient = 0; coefficient < beta_dimension; ++coefficient)
      beta_final_step(comp, coefficient) = step[coefficient];
  }
  beta_attempts.attr("names") = component_labels;
  beta_acceptances.attr("names") = component_labels;
  beta_warmup_attempts.attr("names") = component_labels;
  beta_warmup_acceptances.attr("names") = component_labels;
  beta_retained_attempts.attr("names") = component_labels;
  beta_retained_acceptances.attr("names") = component_labels;
  beta_acceptance_rate.attr("names") = component_labels;
  beta_retained_acceptance_rate.attr("names") = component_labels;
  beta_final_step.attr("dimnames") = List::create(
    component_labels,
    CharacterVector::create("intercept", "distance", "GC", "TES", "ACC"));

  return List::create(Named("chains") = chains_out,
                      Named("gamma") = gamma_out,
                      Named("theta") = theta_out,
                      Named("size") = size_out,
                      Named("z_final") = z_reported,
                      Named("z_checkpoints") = z_checkpoints,
                      Named("z_probabilities") = z_probabilities,
                      Named("z_probability_draws") = z_probability_draws,
                      Named("z_probability_burnin") = z_probability_burnin,
                      Named("beta_mixing") = List::create(
                        Named("proposal") =
                          "adaptive-scale QR-whitened random walk",
                        Named("adaptation_end") = eb_freeze_iteration,
                        Named("target_acceptance_lower") = accept_lo,
                        Named("target_acceptance_upper") = accept_hi,
                        Named("attempts") = beta_attempts,
                        Named("acceptances") = beta_acceptances,
                        Named("acceptance_rate") = beta_acceptance_rate,
                        Named("warmup_attempts") = beta_warmup_attempts,
                        Named("warmup_acceptances") =
                          beta_warmup_acceptances,
                        Named("retained_attempts") = beta_retained_attempts,
                        Named("retained_acceptances") =
                          beta_retained_acceptances,
                        Named("retained_acceptance_rate") =
                          beta_retained_acceptance_rate,
                        Named("final_qr_step") = beta_final_step,
                        Named("separate_from_branch_swaps") = true),
                      Named("branch_mixing") = List::create(
                        Named("enabled") = branch_swap_enabled,
                        Named("interval") = branch_swap_interval,
                        Named("attempts") = branch_swap_attempts,
                        Named("acceptances") = branch_swap_accepts,
                        Named("acceptance_rate") = branch_swap_attempts > 0 ?
                          (double)branch_swap_accepts /
                            (double)branch_swap_attempts : NA_REAL,
                        Named("toward_intended_attempts") =
                          branch_swap_toward_intended_attempts,
                        Named("toward_intended_acceptances") =
                          branch_swap_toward_intended_accepts,
                        Named("warmup_attempts") = branch_swap_warmup_attempts,
                        Named("warmup_acceptances") = branch_swap_warmup_accepts,
                        Named("retained_attempts") =
                          branch_swap_attempts - branch_swap_warmup_attempts,
                        Named("retained_acceptances") =
                          branch_swap_accepts - branch_swap_warmup_accepts,
                        Named("joint_state") =
                          "all component-2/3 labels, betas and dispersions",
                        Named("joint_proposal") = true,
                        Named("immediate_parameter_refresh") = true,
                        Named("branch_gap") = branch_gap_out),
                      Named("signal_block_mixing") = List::create(
                        Named("enabled") =
                          !ablation_disable_blockmove &&
                            signal_block_move_interval > 0,
                        Named("interval") = signal_block_move_interval,
                        Named("attempts") = signal_block_move_attempts,
                        Named("acceptances") = signal_block_move_accepts,
                        Named("acceptance_rate") =
                          signal_block_move_attempts > 0 ?
                            (double)signal_block_move_accepts /
                              (double)signal_block_move_attempts : NA_REAL,
                        Named("warmup_attempts") = signal_block_warmup_attempts,
                        Named("warmup_acceptances") = signal_block_warmup_accepts,
                        Named("retained_attempts") =
                          signal_block_move_attempts - signal_block_warmup_attempts,
                        Named("retained_acceptances") =
                          signal_block_move_accepts - signal_block_warmup_accepts,
                        Named("immediate_parameter_refresh") = true),
                      Named("dispersion_mixing") = List::create(
                        Named("proposal") = "adaptive log-scale random walk",
                        Named("adaptation_end") = eb_freeze_iteration,
                        Named("target_acceptance_lower") = size_accept_lo,
                        Named("target_acceptance_upper") = size_accept_hi,
                        Named("attempts") = size_attempts,
                        Named("acceptances") = size_acceptances,
                        Named("acceptance_rate") = size_acceptance_rate,
                        Named("final_log_step") = size_final_logstep),
                      // Recorded so every fit carries the identifiability
                      // settings it was produced under, and so sensitivity
                      // runs can be told apart without consulting the call.
                      Named("comp23_barrier") = List::create(
                        Named("kappa") = comp23_barrier_kappa,
                        Named("width") = comp23_barrier_w,
                        Named("enabled") = !ablation_disable_barrier),
                      Named("noise_relationship") = List::create(
                        Named("enabled") = use_noise_relationship_prior,
                        Named("link_sd") = noise_link_sd,
                        Named("order_strength") = noise_order_strength,
                        Named("order_width") = noise_order_width,
                        Named("slope_scale") =
                          "standardised log1p-covariate effects",
                        Named("component3_coupled") = true,
                        Named("component2_unrestricted") = true,
                        Named("internal_labels_exchangeable") = false),
                      Named("regression_prior") = List::create(
                        Named("method") = use_data_priors ?
                          "soft_allocation_empirical_Bayes_warmup_frozen" :
                          "user_fixed_Normal",
                        Named("reestimated_each_iteration") = false,
                        Named("soft_allocation_updates") = use_data_priors,
                        Named("updates_completed") = eb_prior_updates,
                        Named("freeze_iteration") = eb_freeze_iteration,
                        Named("fixed_for_retained_draws") = true,
                        Named("frozen_hyperparameters") = active_prior_object,
                        Named("centre") = use_data_priors ?
                          "component-wise weighted ridge regression of log1p(count)" :
                          "user supplied",
                        Named("small_component_fallback_cells") = 8),
                      Named("dispersion_prior") = List::create(
                        Named("family") = "Gamma",
                        Named("shape") = NumericVector::create(3.0, 2.0, 2.0),
                        Named("rate") = NumericVector::create(1.0, 0.2, 0.2)),
                      // Records the full-matrix emission convention. Symmetry
                      // remains available for classification reflection only.
                      Named("pair_weighting") = List::create(
                        Named("symmetric_input") = hicpotts_symmetric_pairs(y),
                        Named("offdiagonal_weight") = 1.0,
                        Named("likelihood_domain") = "complete_N_by_N_matrix"),
                      // Batch means of the post-burn-in latent-state
                      // frequencies, for autocorrelation-aware per-cell MCSE
                      // and allocation ESS. Summary only: the three-way MAP
                      // label is computed from z_probabilities and is not
                      // affected by anything derived from these.
                      Named("z_probability_batches") = List::create(
                        Named("batch_means") = z_batch_means,
                        Named("batch_draws") = z_batch_sizes,
                        Named("n_batches") = z_batches_used,
                        Named("nominal_batch_size") = z_batch_size),
                      Named("performance") = List::create(
                        Named("engine") = "cached_native_v1",
                        Named("cached_log1p_covariates") = true,
                        Named("direct_checkerboard_allocation") = true,
                        Named("native_parameter_posterior") = true,
                        Named("preallocated_chain_storage") = true,
                        Named("elapsed_seconds") = seconds_total,
                        Named("seconds_per_iteration") =
                          seconds_total / std::max(1, executed_iterations),
                        Named("branch_seconds") = seconds_branch,
                        Named("allocation_seconds") = seconds_allocation,
                        Named("block_move_seconds") = seconds_block_move,
                        Named("empirical_bayes_seconds") = seconds_empirical_bayes,
                        Named("parameter_seconds") = seconds_parameters,
                        Named("gamma_seconds") = seconds_gamma,
                        Named("theta_seconds") = seconds_theta,
                        Named("other_seconds") =
                          std::max(0.0, seconds_total - seconds_profiled)),
                      Named("sampler_settings") = List::create(
                        Named("N") = N,
                        Named("distribution") = dist,
                        Named("iterations_requested") = iterations,
                        Named("iterations_completed") = executed_iterations,
                        Named("verbose") = verbose,
                        Named("progress_interval") = progress_interval,
                        Named("z_probability_burnin") = z_probability_burnin,
                        Named("gamma_method") = gamma_method,
                        Named("gamma_update_interval") = gamma_update_interval,
                        Named("abc_potts_sweeps") = abc_potts_sweeps,
                        Named("abc_sim_reps") = abc_sim_reps,
                        Named("abc_epsilon") = abc_epsilon,
                        Named("abc_epsilon_user_supplied") =
                          std::isfinite(abc_user_epsilon),
                        Named("abc_epsilon_shared_calibration") = false,
                        Named("abc_calibration_reps_executed") =
                          abc_calibration_reps_executed,
                        Named("abc_simulator_engine") =
                          "compact_lookup_v1",
                        Named("validation_freeze_z") = validation_freeze_z,
                        Named("validation_freeze_beta") = validation_freeze_beta,
                        Named("validation_freeze_gamma") = validation_freeze_gamma,
                        Named("validation_freeze_size") = validation_freeze_size,
                        Named("validation_freeze_theta") = validation_freeze_theta));
}
