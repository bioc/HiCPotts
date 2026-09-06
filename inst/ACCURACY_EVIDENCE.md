# Evidence protocol for HiCPotts estimation and allocation

Engineering equivalence and statistical recovery answer different questions.
The following tests must be reported separately and in this order.

1. **Native/reference equivalence.** Compare the cached C++ component
   posterior and cell-wise allocation probabilities with
   `posterior_combined()` and `pz_123()` for every emission family. Failure is
   an implementation defect; success only proves that the faster engine uses
   the same numerical target.

2. **Oracle-label parameter recovery.** Hold simulated labels at their truth
   and estimate beta, dispersion and zero inflation. Failure localises the
   problem to the emission likelihood, covariate design, priors or parameter
   updates; changing allocation moves cannot repair it.

3. **Oracle-parameter allocation recovery.** Hold beta, dispersion, theta and
   gamma at their generating values and update or score labels. Failure means
   the three components are not distinguishable under the stated emission and
   Potts model for that scenario; longer chains cannot create missing
   information.

4. **Joint-fit ablations.** Repeat the fit with gamma, dispersion and theta
   fixed at truth, one block at a time. The first fixed block that restores
   recovery identifies the coupled parameter responsible for the failure.

5. **Mixing evidence.** Run at least four dispersed chains and inspect
   rank-normalised split-Rhat, bulk/tail ESS, MCSE, component occupancy,
   component-2/3 branch acceptance and agreement of the full allocation maps.
   Poor diagnostics call for more effective transitions or computation, not an
   accuracy claim.

6. **Known-truth recovery.** Use `validate_hicpotts_simulation()` across at
   least 30 replicates per scenario. Report coefficient and dispersion bias,
   RMSE, 95% interval coverage, class-specific sensitivity and specificity,
   calibration, and between-chain allocation disagreement—not overall accuracy
   alone.

7. **Identifiability grid.** Include noise-dominant 80/12/8 and 90/7/3
   proportions, exact and approximate component-1/3 slope sharing, several
   dispersions, N = 20/30/50, and multiple gamma and zero-inflation values.
   Record where component-2/3 separation becomes possible or impossible.

8. **Version-controlled comparison.** Run identical generated datasets and
   seeds under the published baseline and candidate release. An improvement is
   supported only when recovery metrics improve across held-out scenarios and
   Monte Carlo uncertainty, not merely on the simulations used to choose a
   setting.

The diagnostic interpretation is decisive: oracle-label failure implicates
parameter estimation; oracle-parameter failure implicates allocation
identifiability; both oracle tests passing while the joint fit fails implicates
posterior coupling or mixing; persistent bias with good multi-chain mixing
implicates model/prior mismatch and requires additional scientific information
or a methodological change.
