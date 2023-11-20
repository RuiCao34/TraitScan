#' TraitScan Test Using Summary-level GWAS Data
#'
#' Run HC and TC tests on summary-level GWAS data and combine tests by minp strategy.
#'
#' @param beta_hat a numeric vector, estimated regression coefficients of traits.
#' @param se_hat a numeric vector, estimated standard errors of beta_hat. When only z scores are available, it is also acceptable to input z scores into \code{beta_hat} and \code{rep(1,length(beta_hat))} into \code{se_hat}.
#' @param null_cov a \code{p} by \code{p} positive definite matrix, the covariance matrix of z scores under the null. Can be estimated from null SNPs or LD score regression.
#' @param alpha_grid the alpha grid for p value thresholds. Should be the same as \code{alpha_grid} when calling \code{TraitScan.MC}
#' @param null_dist the output of function \code{TraitScan.MC}. Note the null distribution only depends on \code{alpha_grid} and \code{p}. Also downloadable from https://drive.google.com/drive/folders/1qZsJLrrvpMERxvqeCByDpNhYI7PccbUw?usp=sharing. All downloadable null distributions were simulated under 1e4 iterations, and the alpha_grid = sort(exp(seq(log(0.05), log(0.05/p), length.out=200))).
#' @param tol When \code{null_cov} is rank-deficient, it is recommended that the eigenvalues smaller than \code{tol} are filtered out. Additional sensitivity analysis may also be carried on \code{null_cov} to ensure the output is replicable.
#'
#' @return A list with the selected subset of traits
#' \item{trait.pval.decor}{The decorrelated GWAS p-values for each trait.}
#' \item{select.traits.minp}{The selected traits by minp strategy of HC and TC. The value represented the order of selected traits in \code{beta_hat}.}
#' \item{select.traits.HC}{The selected traits by HC.}
#' \item{select.traits.TC}{The selected traits by TC.}
#' \item{tests.minp}{Which test has a smaller p value, either HC or TC.}
#' \item{p.value}{The empirical p-value of combined test from Monte Carlo}
#' @author Rui Cao, Tianzhong Yang
#' @examples
#'
#' set.seed(1)
#' p = 10    # number of traits
#' alpha_grid = sort(exp(seq(log(0.05), log(0.05/p), length.out=200)))    # alpha grid in HC statistic
#' nsim = 1e4    # number of MC iterations
#' null_dist <- TraitScan.MC(nsim = nsim, p = p, alpha_grid = alpha_grid)
#'
#' beta_hat = c(rnorm(8,0,0.1), 0.3, 0.3)
#' se_hat = rep(0.1, p)
#' null_cov = matrix(0.2, nrow = p, ncol = p)
#' diag(null_cov) = 1
#'
#' example <- TraitScan.main(beta_hat, se_hat, null_cov, alpha_grid, null_dist)
#' @references
#' Cao, Rui, et al. "Subset scanning for multi-trait analysis using GWAS summary statistics." medRxiv (2023): 2023-07.
#' @rdname summary_test
#' @export
#'

TraitScan.main <- function(beta_hat, se_hat, null_cov, alpha_grid, null_dist, tol = 1e-2){

  rank_null_cov = Matrix::rankMatrix(null_cov)
  if(rank_null_cov < nrow(null_cov)){
    warning("The covariance matrix is not full-rank. Filtering small eigenvalues...")
    P = cov2cor(null_cov)
    V = diag(null_cov)
    temp = eigen(P)
    temp$values[temp$values < tol] = 0
    inverse_square = 1/sqrt(temp$values)
    inverse_square[inverse_square == Inf] = 0
    w = temp$vectors %*% diag(inverse_square) %*% solve(temp$vectors) %*% sqrt(diag(1/V))
  }else{
    w <- whitening::whiteningMatrix(null_cov, method='ZCA-cor')
  }
  z.decor <- as.vector(beta_hat/se_hat %*% t(w))
  pval.de <- (1-pnorm(abs(z.decor)))*2

  interest.data <- list("summary" = cbind(z.decor, z.decor, z.decor))
  test.output <- combined.test(interest.data, diag(1,length(z.decor)), null_dist$null_dist_TC_prior, alpha_grid, null_dist$null_dist_TC_stat)
  p.value <- sum(null_dist$null_dist_combined_test <= test.output$test.statistics) / length(null_dist$null_dist_combined_test)
  if(p.value == 0) p.value <- paste0("<", 1/length(null_dist$null_dist_combined_test))
  return( list(trait.pval.decor = pval.de,
               select.traits.minp = test.output$select,
               select.traits.HC = test.output$select.HC,
               select.traits.TC = test.output$select.TC,
               tests.minp = test.output$minp.test,
               p.value = p.value
               ) )

}
