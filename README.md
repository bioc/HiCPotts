# HiCPotts <img src="man/figures/logo.png" height="110" align="right" />

*A Bayesian framework for detection of enriched Hi‑C interactions and experimental biases in Hi-C data*  

[![R‑CMD‑check](https://github.com/igosungithub/HiCPotts/actions/workflows/R‑CMD‑check.yml/badge.svg)](https://github.com/igosungithub/HiCPotts/actions) &nbsp;
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

---

## Development version

A newer development version of HiCPotts is available through Bioconductor devel.

Users who require the latest features and updates can install the development version by first switching to Bioconductor devel:

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(version = "devel")
BiocManager::install("HiCPotts")

---

## 1  Why HiCPotts?

Hi‑C interaction counts are  

* **sparse** (most bin pairs are zero),  
* **over‑dispersed** and often **zero‑inflated**, and  
* **spatially correlated** along the genome.

**HiCPotts** deals with these challenges by combining

| Layer | Details |
|-------|---------------------------------------|
| **Spatial layer** | **Potts model** (Wu 1982) with interaction parameter \(\gamma\) to capture neighbourhood dependence. |
| **Count layer** | Mixture of three components (“noise”, “signal”, “false‑positive”) modelled with Poisson | NB | ZIP | ZINB. |
| **Bias regression** | Genomic distance, GC content, TE density and chromatin accessibility are covariates in a log‑linear model for the mean. |
| **Inference** | Metropolis‑within‑Gibbs MCMC updates betas, \(\gamma\), zero‑inflation \(\theta\) and dispersion *size* (if NB/ZINB). |

---

## 2  Installation

```r
## released version (once on Bioconductor)
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("HiCPotts")

Requirements: R ≥ 4.2, C++17 compiler, plus
Rcpp, RcppArmadillo, parallel.
```
The optional genome support uses BSgenome.Dmelanogaster.UCSC.dm6 (only needed to compute GC content from raw coordinates, hence moved to Suggests):

```r
BiocManager::install("BSgenome.Dmelanogaster.UCSC.dm6")
```

## 3  Quick start

```r
library(HiCPotts)


## 1  Load long‑format interaction table
df <- read.csv("hic_interactions.csv")  # start, end, interactions, GC, ACC, TES

## 2  Convert to N×N matrices
prep   <- process_data(df, N = 40, standardization_y = FALSE)
x_vars <- prep$x_vars
y_list <- prep$y                      # list of 40×40 matrices

## 3  Run the three‑component MCMC
res <- run_chain_betas(
  N            = 40,
  iterations   = 5000,
  x_vars       = x_vars,
  y_sim_list   = y_list,
  use_data_prior = TRUE,
  dist         = "ZINB",
  size_initial = c(1, 1, 1), # initial NB size for the 3 comps
  thetap       = 0.5,        # initial θ
  mc_cores     = 1
)

## 4  Posterior component probabilities
probs <- compute_HMRFHiC_probabilities(
  data        = df,
  chain_betas = results,                    # full list, NOT results[["chains"]]
  iterations  = 2000,
  dist        = "ZINB"
)
head(probs)
```

```mermaid
graph TD
  A[get_data] --> B[process_data]
  B --> C[run_chain_betas]
  C --> D[Posterior summaries + probability maps]
```


Key exported functions

| Function                                                            | What it does                                                                                                 |
| ------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------ |
| `process_data()`                                                    | Converts long‑format counts + covariates to lists of $N \times N$ matrices (counts optionally scaled 1–500). |
| `run_chain_betas()`                                                 | Parallel front‑end to `run_metropolis_MCMC_betas()` for multiple matrices; **always uses three components**. |
| `run_metropolis_MCMC_betas()`                                       | Core Metropolis‑within‑Gibbs sampler.                                                                        |
| `compute_HMRFHiC_probabilities()`                                   | Converts posterior mean betas into per‑pixel probabilities for components 1, 2, 3.                           |
| `prior_combined()`, `likelihood_combined()`, `posterior_combined()` | Compute log‑priors, log‑likelihoods and log‑posteriors.                                                      |
| `pred_combined()`                                                   | Calculates log‑linear predictors $\lambda$.                                                                  |
| `Neighbours_combined()` (C++)                                       | Counts matching neighbours for the Potts model.                                                              |
| `pz_123()` (C++)                                                    | Site‑wise component log‑posterior calculation.                                                               |

```r
mcmc1 <- res[[1]]

## trace of γ
plot(mcmc1$gamma, type = "l", col = "#1f77b4",
     main = "Potts interaction γ", ylab = "γ")

## posterior means of βs for component 1
colMeans(mcmc1$chains[[1]][-(1:2500), ])
```
Advanced options
Distribution choice – dist = "Poisson", "NB" or "ZIP".

Fixed vs. data‑driven priors – set use_data_priors = FALSE in the sampler and supply user_fixed_priors.

Proposal tuning – edit the sd_values vectors inside run_metropolis_MCMC_betas().

```r
set.seed(1)
N <- 10
fake <- data.frame(
  start        = rep(1:N, each = N),
  end          = rep(1:N, N),
  interactions = rpois(N*N, 5),
  GC           = runif(N*N),
  ACC          = runif(N*N),
  TES          = runif(N*N)
)

prep   <- process_data(fake, N)
res    <- run_chain_betas(prep$y, N = N, iterations = 100,
                          x_vars = prep$x_vars,
                          thetap = 0.5, size_initial = c(1,1,1),
                          dist   = "Poisson", mc_cores = 1)
```

Another example is using the test_data2.csv file inside the folder; inst/extdata

```r
mydata=read.csv("~/inst/extdata/test_data2.csv")

N=20
thetap=0.6

iterations=2000


# Assuming mydata is a data frame
colnames(mydata)[colnames(mydata) %in% c("start.i.", "start.i", "start")] <- "start"
colnames(mydata)[colnames(mydata) %in% c("end.j.", "end.j", "end")] <- "end"
colnames(mydata)[colnames(mydata) %in% c("Acc", "ACC")] <- "ACC"
colnames(mydata)[colnames(mydata) %in% c("GC", "Gc")] <- "GC"
colnames(mydata)[colnames(mydata) %in% c("Tes", "TEs", "TES")] <- "TES"
colnames(mydata)[colnames(mydata) %in% c("interactions", "interaction")] <- "interactions"



scaled_data<-process_data(mydata, N, scale_max = 500, standardization_y = TRUE)


chains=run_chain_betas(
  N,
  gamma_prior=0.3,
  iterations=iterations,
  x_vars=scaled_data[["x_vars"]],
  y=scaled_data[["y"]],
  use_data_priors=TRUE,
  dist="Poisson",
  distance_metric="manhattan",
  mc_cores = 22
)
```
## 4  Important changes from 1.0.0 -> 1.0.1
The 1.0.1 patch release fixes a number of statistical and plumbing bugs surfaced by external review. Most user code will continue to run unchanged. Highlights:
 
* **Component 2 and 3 likelihoods** now use every data point assigned to the component, and does not force stop if a component was not assigned any data during iteration.
* **ZIP/ZINB zero probability** written for more clarify to $\theta + (1-\theta)\exp(-\mu)$.
* **Neighbour counting** in `Neighbours_combined()` now returns counts in $\{0,1,2,3,4\}$ as intended.
* **ABC update for $\gamma$** simulates synthetic data with the proposed $\gamma$ with a more refined distribution.
* **MH acceptance** for regression coefficients no longer includes a "proposal density" term; the symmetric Gaussian random-walk proposal cancels correctly.
* **Data-driven priors** are more deterministic given the same (parameters, data, assignments) and refined distribution.
* `compute_HMRFHiC_probabilities()` no longer silently caps interactions at 500; pass `max_interactions = 500` to restore the old behaviour.
* Windows: `mc_cores > 1` emits a warning and falls back to 1 (previously silently).
* `.hic` input to `get_data()` now raises an informative error (previously failed inside `h5ls`).
See `NEWS` for the full list.
 
---
 
## 5. Feedback and bug reports
 
Please file issues at <https://github.com/igosungithub/HiCPotts/issues>.
 
License
GPL-3 © 2025 Itunu Osuntoki
