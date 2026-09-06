#ifndef HICPOTTS_ABC_S1S2_GAMMA_H
#define HICPOTTS_ABC_S1S2_GAMMA_H

// -----------------------------------------------------------------------------
// ABC support for the Potts coupling gamma.
//
// gamma enters the model ONLY through the Potts prior on the latent field:
//
//     p(z | gamma) = exp( gamma * S(z) ) / Z(gamma),
//     S(z) = number of same-label four-neighbour adjacent pairs.
//
// S(z) is the sufficient statistic for gamma, and Z(gamma) is the intractable
// normalising constant that forces a likelihood-free update in the first place.
//
// Two design points follow, and together they fix both defects of the previous
// implementation:
//
// 1. The simulator is p(z | gamma) -- a Potts draw at the candidate gamma.
//    It depends on NOTHING except gamma and the lattice size. It therefore
//    cannot drift as beta/size/theta move, so a retained auxiliary state is
//    always a valid draw from exactly the distribution the acceptance ratio
//    uses. (The previous simulator drew COUNTS using the running
//    beta/size/theta, so p(x | gamma) changed every sweep and the retained
//    state and the proposal came from different distributions.)
//
// 2. The summary statistic is S(z) rather than a join-count excess of
//    simulated COUNTS. Comparing simulated counts to observed counts routed
//    gamma's signal through the entire emission model, which diluted it to the
//    point where gamma was barely identified: on simulated data with true
//    gamma of 0.30 and 0.60 the old scheme returned posterior means of ~0.50
//    and ~0.55, i.e. essentially the prior. S(z) is the exact sufficient
//    statistic, so the comparison is made in the one quantity that actually
//    carries information about gamma.
//
// The conditioning value is S(z_current), the sufficient statistic of the
// sampler's own current latent field. This makes the step a likelihood-free
// update of the CONDITIONAL gamma | z, which is the correct conditional in the
// surrounding Gibbs scheme. S(z_current) changes as z moves, so the retained
// auxiliary statistic is re-scored against the current field each sweep; the
// retained draw itself is untouched, which is what validity requires.
// -----------------------------------------------------------------------------

// Reusable compact workspace for the auxiliary Potts simulator. The previous
// implementation allocated an R numeric matrix for every simulated field and
// recomputed boundaries, checkerboard parity and three exponentials at every
// site update. At N=50 that meant billions of avoidable operations per chain.
//
// This engine preserves the exact transition: uniform random initialization,
// the same sweep/colour/i/j site order, the same four-neighbour conditional
// probabilities and exactly one R RNG draw per initialized or updated site.
// Only the representation and invariant calculations are cached.
class AbcPottsWorkspace {
 private:
  int N_;
  int cells_;
  std::vector<unsigned char> labels_;
  std::vector<std::array<int, 4> > neighbours_;
  std::array<std::vector<int>, 2> checkerboard_;
  std::array<double, 125> probability1_;
  std::array<double, 125> probability12_;
  std::array<double, 125> probability_total_;
  double prepared_gamma_;

  static inline int pattern_index(int n1, int n2, int n3) {
    return (n1 * 5 + n2) * 5 + n3;
  }

  void prepare_probabilities(double gamma) {
    if (gamma == prepared_gamma_) return;
    for (int n1 = 0; n1 <= 4; ++n1) {
      for (int n2 = 0; n2 <= 4; ++n2) {
        for (int n3 = 0; n3 <= 4; ++n3) {
          const double lp1 = gamma * n1;
          const double lp2 = gamma * n2;
          const double lp3 = gamma * n3;
          const double m = std::max(lp1, std::max(lp2, lp3));
          const double p1 = std::exp(lp1 - m);
          const double p2 = std::exp(lp2 - m);
          const double p3 = std::exp(lp3 - m);
          const int pattern = pattern_index(n1, n2, n3);
          probability1_[pattern] = p1;
          probability12_[pattern] = p1 + p2;
          probability_total_[pattern] = p1 + p2 + p3;
        }
      }
    }
    prepared_gamma_ = gamma;
  }

  inline unsigned char draw_label(int n1, int n2, int n3) const {
    const int pattern = pattern_index(n1, n2, n3);
    const double u = R::runif(0.0, probability_total_[pattern]);
    if (u < probability1_[pattern]) return 1;
    if (u < probability12_[pattern]) return 2;
    return 3;
  }

 public:
  explicit AbcPottsWorkspace(int N)
      : N_(N), cells_(N * N), labels_(N * N), neighbours_(N * N),
        prepared_gamma_(NA_REAL) {
    checkerboard_[0].reserve((cells_ + 1) / 2);
    checkerboard_[1].reserve(cells_ / 2);
    for (int i = 0; i < N_; ++i) {
      for (int j = 0; j < N_; ++j) {
        const int index = i * N_ + j;
        checkerboard_[(i + j) & 1].push_back(index);
        neighbours_[index] = {{
          i > 0 ? index - N_ : -1,
          i + 1 < N_ ? index + N_ : -1,
          j > 0 ? index - 1 : -1,
          j + 1 < N_ ? index + 1 : -1
        }};
      }
    }
  }

  void simulate(double gamma, int sweeps) {
    prepare_probabilities(gamma);

    // Preserve the historical i-then-j initialization order exactly.
    for (int i = 0; i < N_; ++i)
      for (int j = 0; j < N_; ++j)
        labels_[i * N_ + j] =
          static_cast<unsigned char>(std::floor(R::runif(1.0, 4.0)));

    for (int sweep = 0; sweep < sweeps; ++sweep) {
      for (int colour = 0; colour <= 1; ++colour) {
        const std::vector<int> &sites = checkerboard_[colour];
        for (size_t position = 0; position < sites.size(); ++position) {
          const int index = sites[position];
          int n1 = 0, n2 = 0, n3 = 0;
          const std::array<int, 4> &site_neighbours = neighbours_[index];
          for (int direction = 0; direction < 4; ++direction) {
            const int neighbour = site_neighbours[direction];
            if (neighbour < 0) continue;
            const unsigned char label = labels_[neighbour];
            if (label == 1) ++n1;
            else if (label == 2) ++n2;
            else ++n3;
          }
          labels_[index] = draw_label(n1, n2, n3);
        }
      }
    }
  }

  double agreement() const {
    double same = 0.0;
    int pairs = 0;
    int c1 = 0, c2 = 0, c3 = 0;
    for (int i = 0; i < N_; ++i) {
      for (int j = 0; j < N_; ++j) {
        const int index = i * N_ + j;
        const unsigned char label = labels_[index];
        if (label == 1) ++c1;
        else if (label == 2) ++c2;
        else ++c3;
        if (i + 1 < N_) {
          ++pairs;
          if (labels_[index + N_] == label) same += 1.0;
        }
        if (j + 1 < N_) {
          ++pairs;
          if (labels_[index + 1] == label) same += 1.0;
        }
      }
    }
    if (pairs == 0) return 0.0;
    const double n = static_cast<double>(cells_);
    const double p1 = c1 / n, p2 = c2 / n, p3 = c3 / n;
    const double expected = static_cast<double>(pairs) *
      (p1 * p1 + p2 * p2 + p3 * p3);
    return (same - expected) / static_cast<double>(pairs);
  }

  NumericMatrix as_matrix() const {
    NumericMatrix result(N_, N_);
    for (int i = 0; i < N_; ++i)
      for (int j = 0; j < N_; ++j)
        result(i, j) = labels_[i * N_ + j];
    return result;
  }
};

static inline double simulate_potts_abc_summary(
    AbcPottsWorkspace &workspace, double gamma, int sweeps) {
  workspace.simulate(gamma, sweeps);
  return workspace.agreement();
}

// Matrix-returning compatibility path used by the public truth simulator.
static NumericMatrix simulate_potts_abc(int N, double gamma, int sweeps) {
  AbcPottsWorkspace workspace(N);
  workspace.simulate(gamma, sweeps);
  return workspace.as_matrix();
}

// Spatial-clustering statistic for gamma: the EXCESS same-label four-neighbour
// adjacency of a labelling, over and above what its own marginal label
// proportions already imply.
//
//   excess = ( same_pairs - pairs * sum_k p_k^2 ) / pairs
//
// The correction is essential. The raw same-label proportion is dominated by
// composition, not by coupling: a field that is 85% noise shows ~0.73
// agreement even with gamma = 0, whereas a Potts draw from the simulator has
// roughly equal label proportions and so shows ~0.33 at the same gamma.
// Matching the two on the raw proportion drives gamma to its upper bound
// regardless of the truth (empirically ~0.99 for every true gamma from 0.2 to
// 0.8). Subtracting the composition-implied expectation removes that
// confound and leaves a statistic that responds to spatial coupling alone,
// which is what gamma controls. Each unordered pair is visited once (right and
// down neighbours only), and normalising by the pair count makes the statistic
// lattice-size free so the ABC bandwidth means the same thing across N.
static double potts_label_agreement(const NumericMatrix &z) {
  const int N = z.nrow();
  if (N < 2) return 0.0;
  double same = 0.0;
  int pairs = 0;
  int c1 = 0, c2 = 0, c3 = 0;
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      const int lij = (int) z(i, j);
      if (lij == 1) c1++; else if (lij == 2) c2++; else c3++;
      if (i + 1 < N) { pairs++; if ((int) z(i + 1, j) == lij) same += 1.0; }
      if (j + 1 < N) { pairs++; if ((int) z(i, j + 1) == lij) same += 1.0; }
    }
  }
  if (pairs == 0) return 0.0;
  const double n = (double) (N * N);
  const double p1 = c1 / n, p2 = c2 / n, p3 = c3 / n;
  const double expected = (double) pairs * (p1 * p1 + p2 * p2 + p3 * p3);
  return (same - expected) / (double) pairs;
}

// Linear interpolation over the empirical quantile of a numeric sample,
// used to calibrate the ABC kernel bandwidth from prior-predictive draws.
static double abc_numeric_quantile(std::vector<double> x, double p) {
  if (x.empty()) return NA_REAL;
  std::sort(x.begin(), x.end());
  double h = (x.size() - 1) * p;
  size_t lo = (size_t) std::floor(h), hi = (size_t) std::ceil(h);
  if (lo == hi) return x[lo];
  return x[lo] + (h - lo) * (x[hi] - x[lo]);
}

struct AbcToleranceCalibration {
  double epsilon;
  double scale;
};

// Prior-predictive tolerance calibration shared by the standalone sampler and
// repeated-chain wrappers. The caller owns RNG isolation: the sampler uses its
// historical fixed substream and R wrappers use withr::with_seed().
static AbcToleranceCalibration calibrate_abc_tolerance(
    AbcPottsWorkspace &workspace,
    double gamma_prior_shape1, double gamma_prior_shape2,
    double epsilon_quantile, int simulations, int sim_reps, int sweeps) {
  std::vector<double> simulated_statistics;
  simulated_statistics.reserve(simulations);
  for (int repetition = 0; repetition < simulations; ++repetition) {
    const double gamma = R::rbeta(gamma_prior_shape1, gamma_prior_shape2);
    double statistic = 0.0;
    for (int replicate = 0; replicate < sim_reps; ++replicate)
      statistic += simulate_potts_abc_summary(workspace, gamma, sweeps);
    simulated_statistics.push_back(statistic / sim_reps);
  }

  const double centre = abc_numeric_quantile(simulated_statistics, 0.5);
  std::vector<double> distances;
  distances.reserve(simulated_statistics.size());
  for (size_t index = 0; index < simulated_statistics.size(); ++index)
    distances.push_back(std::abs(simulated_statistics[index] - centre));

  double scale = abc_numeric_quantile(distances, 0.5);
  if (!(scale > 1e-8) || !std::isfinite(scale)) scale = 0.1;
  double epsilon = abc_numeric_quantile(distances, epsilon_quantile);
  if (!(epsilon > 1e-8) || !std::isfinite(epsilon)) epsilon = 0.5 * scale;
  AbcToleranceCalibration result = {epsilon, scale};
  return result;
}

#endif
