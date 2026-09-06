#ifndef HICPOTTS_PAIR_WEIGHT_H
#define HICPOTTS_PAIR_WEIGHT_H

#include <Rcpp.h>
#include <cmath>

// -----------------------------------------------------------------------------
// Full-matrix emission weighting (shared by every MCMC transition).
//
// HiCPotts analyses the complete N x N lattice returned by process_data().
// Consequently every stored cell contributes once to the emission likelihood,
// including both (i,j) and (j,i) when the processed matrix is symmetric.  The
// Potts field and the emission likelihood therefore use the same full lattice.
// The symmetry detector is retained only as output metadata and for optional
// reflection averaging of posterior classification probabilities.
// -----------------------------------------------------------------------------

// Mirror of .hicpotts_symmetric_pairs() in R/likelihood.R. Exact equality is
// the correct test: mirrored cells are literal copies of one measured value.
// If floating-point drift ever made a genuinely mirrored matrix compare
// unequal the effect is to skip the halving -- the conservative direction.
static inline bool hicpotts_symmetric_pairs(const Rcpp::NumericMatrix &y) {
  const int nr = y.nrow(), nc = y.ncol();
  if (nr != nc || nr < 2) return false;
  for (int j = 0; j < nc; j++) {
    for (int i = 0; i < nr; i++) {
      const double a = y(i, j);
      if (Rcpp::NumericMatrix::is_na(a)) return false;
      if (i < j && a != y(j, i)) return false;
    }
  }
  return true;
}

// Every cell in the structured N x N input contributes exactly once.
static inline double hicpotts_pair_weight(int i, int j, bool symmetric) {
  (void) i;
  (void) j;
  (void) symmetric;
  return 1.0;
}

// Kept as a shared helper so all transitions continue to use one convention.
static inline double hicpotts_weight_density(double dens, double w) {
  (void) w;
  return dens;
}

#endif // HICPOTTS_PAIR_WEIGHT_H
