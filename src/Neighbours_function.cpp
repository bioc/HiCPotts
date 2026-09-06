#include <RcppArmadillo.h>
using namespace Rcpp;

// Helper function to create shifted matrices
// Helper: shift_matrix(data, dir, N)
//
// Returns an N x N matrix whose (i,j) entry is:
//   "down"  : data(i-1, j)     (i.e. the value ABOVE (i,j) in the original)
//   "up"    : data(i+1, j)     (the value BELOW  (i,j) in the original)
//   "right" : data(i,  j-1)    (the value to the LEFT  of (i,j))
//   "left"  : data(i,  j+1)    (the value to the RIGHT of (i,j))
// Border cells (where the neighbour falls off the lattice) are left at 0.
// -----------------------------------------------------------------------------
NumericMatrix shift_matrix(const NumericMatrix &data,
                           std::string shift_direction,
                           int N) {
  NumericMatrix shifted(N, N);
 
  if (shift_direction == "down") {
    for (int i = 1; i < N; i++)
      for (int j = 0; j < N; j++)
        shifted(i, j) = data(i - 1, j);
 
  } else if (shift_direction == "up") {
    for (int i = 0; i < N - 1; i++)
      for (int j = 0; j < N; j++)
        shifted(i, j) = data(i + 1, j);
 
  } else if (shift_direction == "right") {
    for (int i = 0; i < N; i++)
      for (int j = 1; j < N; j++)
        shifted(i, j) = data(i, j - 1);
 
  } else if (shift_direction == "left") {
    for (int i = 0; i < N; i++)
      for (int j = 0; j < N - 1; j++)
        shifted(i, j) = data(i, j + 1);
 
  } else {
    stop("Unknown shift_direction: %s", shift_direction);
  }
 
  return shifted;
}


// -----------------------------------------------------------------------------
// Helper: compute_neighbours(a, b, N)
//   returns an N x N 0/1 matrix with 1 wherever a(i,j) == b(i,j).
// -----------------------------------------------------------------------------
NumericMatrix compute_neighbours(const NumericMatrix &a,
                                 const NumericMatrix &b,
                                 int N) {
  NumericMatrix r(N, N);
  for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++)
      r(i, j) = (a(i, j) == b(i, j)) ? 1.0 : 0.0;
  return r;
}

// -----------------------------------------------------------------------------
//Add Note:
//
// For each site (i,j) we want the COUNT of its four (up/down/left/right)
// neighbours whose state equals the state at (i,j) (or at proposed_value(i,j)
// when a proposed configuration is supplied). 
//
// The reference configuration (the "centre" of each site) is:
//   - proposed_value, if supplied
//   - potts_data, otherwise
// The four neighbour configurations always come from potts_data (the *current*
// Potts state). This is the usual single-site update convention.
//
// Border cells naturally get a 0 contribution from the missing side: the
// shifted matrix is 0 at the border, and the reference takes values in
// {1,2,3}, so the comparison is false.
// -----------------------------------------------------------------------------

// Rcpp automatically wraps every `[[Rcpp::export]]` function's generated
// wrapper in an `Rcpp::RNGScope` (a GetRNGstate()/PutRNGstate() bracket)
// unless told otherwise. `run_metropolis_MCMC_betas()` (src/MCMC_chain.cpp)
// calls this function via a `Function` R-callback object from within its own
// already-RNGScope-wrapped body, every MCMC iteration; a nested RNGScope from
// an inner exported call would otherwise disrupt the outer function's random
// draws for the remainder of that top-level call. `Neighbours_combined` and
// `pz_123` (below) use no random-number generation internally (both are
// purely deterministic), so marking them `rng = false` suppresses Rcpp's
// auto-inserted RNGScope in their wrappers, avoiding the nested scope and
// preserving independently-advancing RNG behaviour in the caller.
//' @name Neighbours_combined
//' @title Neighbours Function for the Potts Model
//'
//' @description
//' Computes the number of agreeing neighbors for each site in a Potts model lattice.
//' The Potts model is a generalization of the Ising model, where each site on a lattice can 
//' take on multiple states. By counting how many of a site's four neighbors share the same 
//' state, this function provides a measure of clustering within the lattice.
//'
//' @usage
//' Neighbours_combined(potts_data, N, proposed_value = NULL)
//'
//' @param potts_data A numeric \eqn{N \times N} matrix representing the current configuration 
//'   of the Potts model. Each element corresponds to the state of a particular site.
//'
//' @param N An integer specifying the dimension of the \code{potts_data} lattice (i.e., the 
//'   lattice has \code{N} rows and \code{N} columns).
//'
//' @param proposed_value (Optional) A numeric \eqn{N \times N} matrix representing a proposed 
//'   configuration of the Potts model. If provided, the function computes neighbor relationships 
//'   relative to this proposed configuration; otherwise, it uses a shifted version of \code{potts_data} 
//'   to determine neighbors.
//'
//' @details
//' The function checks each site and compares it with its four directional neighbors 
//' (up, down, left, right). For each neighbor that has the same state, a value of 1 is assigned. 
//' Summing these values for each site results in a measure of how "similar" the immediate neighborhood 
//' is to that site.
//'
//' Internally, the function:
//' \enumerate{
//'   \item Uses \code{shift_matrix} to create versions of the lattice shifted in each of the four directions.
//'   \item Uses \code{compute_neighbours} to determine if neighboring positions match in state.
//'   \item Sums the contributions from all four directions to form the \code{Neighbours_total} matrix.
//' }
//'
//' By providing an optional \code{proposed_value} matrix, this function can be integrated into an 
//' MCMC algorithm where new configurations are proposed and evaluated 
//' based on their local coherence.
//'
//' @return
//' A numeric \eqn{N \times N} matrix, \code{Neighbours_total}, where each element is the count of 
//' how many of the four neighbors match the state of that site.
//'
//' @examples
//' 
//' # Suppose we have a 5x5 Potts model lattice:
//' potts_data <- matrix(sample(1:3, 25, replace = TRUE), ncol = 5)
//' N <- 5
//'
//' # Compute the neighbors without a proposed configuration:
//' neigh <- Neighbours_combined(potts_data, N)
//' neigh
//'
//' # Suppose we propose a new configuration:
//' proposed <- matrix(sample(1:3, 25, replace = TRUE), ncol = 5)
//'
//' # Evaluate neighbors for the proposed configuration:
//' neigh_proposed <- Neighbours_combined(
//'   potts_data, N,
//'   proposed_value = proposed
//' )
//' neigh_proposed
//' 
//' @noRd
//'
// [[Rcpp::export(rng = false)]]
NumericMatrix Neighbours_combined(NumericMatrix potts_data,
                                  int N,
                                  Nullable<NumericMatrix> proposed_value = R_NilValue) {
  NumericMatrix ref = proposed_value.isNotNull()
                       ? as<NumericMatrix>(proposed_value)
                       : potts_data;
 
  // Build the four shifted copies of potts_data (the current state)
  NumericMatrix up    = shift_matrix(potts_data, "down",  N); // neighbour above each cell
  NumericMatrix down  = shift_matrix(potts_data, "up",    N); // neighbour below
  NumericMatrix left  = shift_matrix(potts_data, "right", N); // neighbour to the left
  NumericMatrix right = shift_matrix(potts_data, "left",  N); // neighbour to the right
 
  // Compare the reference configuration to each shifted copy
  NumericMatrix m_up    = compute_neighbours(ref, up,    N);
  NumericMatrix m_down  = compute_neighbours(ref, down,  N);
  NumericMatrix m_left  = compute_neighbours(ref, left,  N);
  NumericMatrix m_right = compute_neighbours(ref, right, N);
 
  // Sum agreements (count in 0..4)
  NumericMatrix total(N, N);
  for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++)
      total(i, j) = m_up(i, j) + m_down(i, j) + m_left(i, j) + m_right(i, j);
 
  return total;
}
 
