# TraitScan
A trait subset method for multi-trait GWAS
TraitScan is a powerful trait subset scan algorithm for mult-trait GWAS. It can also test the association significance between multiple traits and a single genetic variant (SNP). By taking the trait correlations into account, TraitScan is more powerful than PheWAS under scenarios with sparse and moderate true signals. TraitScan also has linear computational time over the number of traits and thus can handle over 1,000 traits. To run the TraitScan algorithm, first a Monte Carlo simulation is required using \code{TraitScan.MC}. With the summary-level GWAS data available, \code{TraitScan.summary.level} can run the test and identify a subset of traits most likely associated with the SNP.
The paper draft of TraitScan is under prepration.

# Required Packages

R >= 4.1.0, "whitening" >= 1.3.0

# Installation
devtools::install_github("RuiCao34/TraitScan")

# Example for summary level data
set.seed(1)
p = 10    # number of traits
alpha_grid = sort(exp(seq(log(0.05), log(0.05/p), length.out=p*2)))    # alpha grid in HC statistic
nsim = 1e4    # number of MC iterations
null_dist <- TraitScan.MC(nsim = nsim, p = p, alpha_grid = alpha_grid)

beta_hat = c(rnorm(8,0,0.1), 0.3, 0.3)
se_hat = rep(0.1, p)
null_cov = matrix(0.2, nrow = p, ncol = p)
diag(null_cov) = 1
