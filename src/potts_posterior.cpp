#include <RcppArmadillo.h>
#include "hicpotts_pair_weight.h"
using namespace Rcpp;

// Helper function to get subset indices where z == comp.
//
// Traverses in column-major order (j outer, i inner) so that the resulting
// (rows, cols) ordering exactly matches R's column-major linear storage
// order, which is what pred_combined()'s R-side logical-matrix subsetting
// (`x_vars[[k]][[1]][mask]`) produces. This keeps each cell's observed count
// and Potts neighbour-agreement count paired with the correct fitted mean
// from pred_combined() in pz_123() below.
List get_indices(const NumericMatrix &z, int comp) {
  std::vector<int> rows, cols;
  for (int j = 0; j < z.ncol(); j++)
    for (int i = 0; i < z.nrow(); i++)
      if (z(i, j) == comp) {
        rows.push_back(i);
        cols.push_back(j);
      }
  return List::create(Named("rows") = wrap(rows), Named("cols") = wrap(cols));
}

// See the matching note above Neighbours_combined()'s definition
// (src/Neighbours_function.cpp): this function is called via a nested
// `Function` R-callback from within run_metropolis_MCMC_betas()'s own
// RNGScope-wrapped body, every iteration. pz_123() itself uses no
// random-number generation internally (purely deterministic density and
// probability computation), so `rng = false` is correct here and avoids the
// nested RNGScope.
//' @title Posterior Function for the Potts Model
//'
//' @description
//' This function computes the posterior of the Potts model combined with an interaction count model 
//' (e.g., from genomic interactions data), incorporating mixture components, neighbor relationships, 
//' and distributional assumptions (Poisson, Negative Binomial, Zero-Inflated Poisson, Zero-Inflated Negative Binomial).
//' It returns the log-probabilities for each site on the lattice under the specified model configuration.
//'
//' @usage
//' pz_123(z, sum_neighbours, y, pred_combined, chains, chain_gamma,
//'        x_vars, theta, size_chain, N, iter, dist)
//'
//' @param z A numeric \eqn{N \times N} matrix representing the component assignments for each site on the lattice.
//'   Each element is a component label (e.g., 1, 2, or 3), indicating which mixture component that interaction currently belongs to.
//'
//' @param sum_neighbours An \eqn{N \times N} numeric matrix representing the number of neighboring sites 
//'   (up, down, left, right) that agree in state with each site. Typically computed by the function 
//'   \code{Neighbours_combined()}.
//'
//' @param y An \eqn{N \times N} numeric matrix of observed interaction counts (e.g., genomic interaction data).
//'   Each element corresponds to the count observed at that site.
//'
//' @param pred_combined An R function (passed from R to C++) that computes predicted values (linear predictors) 
//'   for the specified component. This is often a reference to the \code{pred_combined()} function in R.
//'
//' @param chains A \code{list} of parameter vectors (MCMC chains) for the model. Each element corresponds to 
//'   a component and contains the parameters associated with that component.
//'
//' @param chain_gamma A numeric vector (using ABC algorithm) for the interaction parameter \eqn{\gamma} that governs 
//'   neighbor effects in the Potts model. The \code{iter}-th element of this vector is used in the calculation.
//'
//' @param x_vars A list of covariates (sources of biases), each a nested structure as used by 
//'   \code{pred_combined()}. These covariates are used to compute component-specific predicted values (\eqn{\lambda}).
//'
//' @param theta A numeric scalar representing the zero-inflation parameter. Used when the distribution is 
//'   Zero-Inflated Poisson (ZIP) or Zero-Inflated Negative Binomial (ZINB), typically for the first component.
//'
//' @param size_chain An \eqn{3 \times M} matrix representing the MCMC draws of the size (overdispersion) parameter 
//'   for Negative Binomial or Zero-Inflated Negative Binomial distributions. One row per component, and columns 
//'   correspond to iterations. The value at \code{size_chain[comp - 1, iter]} is used for the specified iteration.
//'
//' @param N An integer specifying the dimension of the lattice (\eqn{N \times N}).
//'
//' @param iter An integer indicating which iteration of the MCMC chains (parameters, \eqn{\gamma}, size) to use 
//'   in the calculation.
//'
//' @param dist A character string specifying the distribution:
//'   \itemize{
//'     \item \code{"Poisson"}: Poisson distribution
//'     \item \code{"NB"}: Negative Binomial distribution
//'     \item \code{"ZIP"}: Zero-Inflated Poisson distribution
//'     \item \code{"ZINB"}: Zero-Inflated Negative Binomial distribution
//'   }
//'
//' @details
//' This function calculates the posterior probabilities (on a log-scale) for each site in an \eqn{N \times N} lattice, 
//' given:
//' \enumerate{
//'   \item The current component assignment matrix \code{z}, indicating which component each site belongs to.
//'   \item The observed interaction counts \code{y}.
//'   \item The neighbor structure encoded by \code{sum_neighbours}.
//'   \item Model parameters drawn from MCMC chains (\code{chains}, \code{chain_gamma}, \code{size_chain}).
//'   \item The specified distribution (\code{dist}).
//' }
//'
//' Each component can follow a Poisson, NB, ZIP, or ZINB distribution. For ZIP and ZINB, zero-inflation is applied 
//' only to the first component. The negative binomial distributions also incorporate the \code{size_value} parameter, 
//' drawn from \code{size_chain}.
//'
//' The \code{pred_combined} function is called to obtain linear predictors \eqn{\lambda} for each component. 
//' Probabilities are computed for each observation under the chosen distribution and adjusted by the Potts model 
//' neighbor effect using \eqn{\exp(\gamma \times \text{sum\_neighbours})}.
//'
//' Finally, the probabilities are logged and returned as an \eqn{N \times N} matrix of log-posteriors.
//'
//' @return
//' A numeric \eqn{N \times N} matrix of log-posterior values for each site under the given model configuration 
//' at the specified iteration.
//'
//' @examples
//' N <- 3
//' z <- matrix(c(1, 1, 2, 3, 1, 3, 2, 3, 1), nrow = N)
//' y <- matrix(rpois(N * N, 3), N, N)
//' x_vars <- list(
//'   list(matrix(runif(N * N, 1, 10), nrow = N)),
//'   list(matrix(runif(N * N, 1, 10), nrow = N)),
//'   list(matrix(runif(N * N, 1, 10), nrow = N)),
//'   list(matrix(runif(N * N, 1, 10), nrow = N))
//' )
//' sum_neighbours <- Neighbours_combined(z, N)
//' chains <- list(matrix(c(1, 0, 0, 0, 0), nrow = 1),
//'                matrix(c(2, 0, 0, 0, 0), nrow = 1),
//'                matrix(c(3, 0, 0, 0, 0), nrow = 1))
//' chain_gamma <- 0.5
//' size_chain <- matrix(c(3, 8, 20), nrow = 3, ncol = 1)
//' pz_vals <- pz_123(z, sum_neighbours, y, pred_combined, chains, chain_gamma,
//'   x_vars, theta = 0, size_chain = size_chain, N = N, iter = 0, dist = "Poisson")
//' head(pz_vals)
//'
//' @noRd
// [[Rcpp::export(rng = false)]]
NumericMatrix pz_123(NumericMatrix z, NumericMatrix sum_neighbours, NumericMatrix y,
                     Function pred_combined, List chains, NumericVector chain_gamma,
                     List x_vars, double theta, NumericMatrix size_chain, int N, int iter,
                     std::string dist) {
  NumericMatrix prob_sum(N, N); // Matrix to store probabilities
  const double epsilon = 1e-300; // near machine-floor; avoids log(0) without biasing

  // Detected once per call, not per cell: the scan is O(N^2) against O(N^2)
  // density evaluations that follow. Mirrors .hicpotts_symmetric_pairs().
  const bool y_pairs_symmetric = hicpotts_symmetric_pairs(y);

  // Loop over components
  for (int comp = 1; comp <= 3; comp++) {
    // Get subset indices where z == comp
    List indices = get_indices(z, comp);
    IntegerVector rows = indices["rows"];
    IntegerVector cols = indices["cols"];
    int num_indices = rows.size();
    if (num_indices == 0) continue;

    // Extract y_comp and sum_neighbours_comp based on indices
    NumericVector y_comp(num_indices), sum_neighbours_comp(num_indices);
    for (int k = 0; k < num_indices; k++) {
      int i = rows[k];
      int j = cols[k];
      y_comp[k] = y(i, j);
      sum_neighbours_comp[k] = sum_neighbours(i, j);
    }

    // Get parameters for current component and iteration
    double gamma_iter = chain_gamma[iter];
    NumericMatrix chain_mat = chains[comp - 1];
    NumericVector chain1 = chain_mat(iter, _);
    NumericVector pred_values = as<NumericVector>(pred_combined(chain1, z, x_vars, comp, N));
    NumericVector lambda = exp(pred_values); // Ensure lambda is a vector and lambda = mu = exp(eta)

    // Retrieve the specific size for the current component and iteration
    double size_value = size_chain(comp - 1, iter);

    // Compute probabilities based on component and distribution
    NumericVector prob_comp(num_indices);
    for (int k = 0; k < num_indices; k++) {
      double pk = 0.0;
      if (dist == "ZIP" && comp == 1) {
        // Zero-Inflated Poisson (ZIP) for component 1. Use lambda[k]
        pk = (y_comp[k] == 0)
               ? theta + (1 - theta) * std::exp(-lambda[k])
               : (1 - theta) * R::dpois(y_comp[k], lambda[k], /*give_log=*/false);
      } else if (dist == "ZINB" && comp == 1) {
        // Zero-Inflated Negative Binomial (ZINB) for component 1
        pk = (y_comp[k] == 0)
               ? theta + (1 - theta) * R::dnbinom_mu(0, size_value, lambda[k], false)
               : (1 - theta) * R::dnbinom_mu(y_comp[k], size_value, lambda[k], false);
      } else if (dist == "Poisson" || (dist == "ZIP" && comp != 1)) {
        // Poisson for all components, and for ZIP if component != 1
        pk = R::dpois(y_comp[k], lambda[k], false);
      } else if (dist == "NB" || (dist == "ZINB" && comp != 1)) {
        // Negative Binomial for all components, and for ZINB if component != 1
        pk = R::dnbinom_mu(y_comp[k], size_value, lambda[k], false);
      } else {
        stop("Invalid distribution specified.");
      }
      // Full-matrix weighting: every stored cell contributes once to the z
      // conditional, matching beta/size/theta and the complete N x N Potts
      // field. See src/hicpotts_pair_weight.h.
      pk = hicpotts_weight_density(
             pk, hicpotts_pair_weight(rows[k], cols[k], y_pairs_symmetric));
      prob_comp[k] = pk;
    }

    // Apply neighbor interaction effect
    // Multiply by the Potts neighbour factor (natural scale), then assign
        for (int k = 0; k < num_indices; k++) {
      double v = std::exp(gamma_iter * sum_neighbours_comp[k]) * prob_comp[k];
      if (!std::isfinite(v) || v < epsilon) v = epsilon;
      prob_sum(rows[k], cols[k]) = v;
    }
  }

  // Clip values in prob_sum to avoid log(0). Log-transform the N x N natural-scale matrix
  for (int i = 0; i < prob_sum.nrow(); i++) {
    for (int j = 0; j < prob_sum.ncol(); j++) {
      double v = prob_sum(i, j);
      prob_sum(i, j) = (v < epsilon) ? std::log (epsilon): std::log(v);
    }
  }
  return prob_sum; // Return log of probabilities
}
