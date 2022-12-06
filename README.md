# TraitScan
A trait subset method for multi-trait GWAS
TraitScan is a powerful trait subset scan algorithm for mult-trait GWAS. It can also test the association significance between multiple traits and a single genetic variant (SNP). By taking the trait correlations into account, TraitScan is more powerful than PheWAS under scenarios with sparse and moderate true signals. TraitScan also has linear computational time over the number of traits and thus can handle over 1,000 traits. To run the TraitScan algorithm, first a Monte Carlo simulation is required using *TraitScan.MC*. With the summary-level GWAS data available, *TraitScan.summary.level* can run the test and identify a subset of traits most likely associated with the SNP.<br />
The paper draft of TraitScan is under prepration.

# Required Packages

R $\ge$ 4.1.0, "whitening" $\ge$ 1.3.0

# Installation
devtools::install_github("RuiCao34/TraitScan")

# Example for summary level data
set.seed(1) <br />
p = 10    # number of traits <br />
alpha_grid = sort(exp(seq(log(0.05), log(0.05/p), length.out=p*2)))    # alpha grid in HC statistic <br />
nsim = 1e4    # number of MC iterations <br />
null_dist <- TraitScan.MC(nsim = nsim, p = p, alpha_grid = alpha_grid) <br />

beta_hat = c(rnorm(8,0,0.1), 0.3, 0.3) <br />
se_hat = rep(0.1, p) <br />
null_cov = matrix(0.2, nrow = p, ncol = p) <br />
diag(null_cov) = 1 <br />
example <- TraitScan.summary.level(beta_hat, se_hat, null_cov, alpha_grid, null_dist)
