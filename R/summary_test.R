#' Monte Carlo Simulation
#'
#' Monte Carlo Simulation for the null distributions of TC and combined tests
#'
#' @param beta_hat a numeric vector, estimated regression coefficients of traits.
#' @param se_hat a numeric vector, estimated standard errors of beta_hat.
#' @param null_cov a \code{p} by \code{p} positive definite matrix, the covariance matrix of z scores under the null. Can be estimated from null SNPs or LD score regression.
#' @param alpha_grid the alpha grid for p value thresholds. Should be the same as \code{alpha_grid} when calling \code{TraitScan.MC}
#' @param null_dist the output of function \code{TraitScan.MC}.
#'
#' @return A list with the selected subset of traits
#' \item{select.traits.minp}{The selected traits by minp strategy of HC and TC. The value represented the order of selected traits in \code{beta_hat}.}
#' \item{select.traits.HC}{The selected traits by HC.}
#' \item{select.traits.TC}{The selected traits by TC.}
#' \item{tests.minp}{Which test has a smaller p value, either HC or TC.}
#' \item{p.value}{The empirical p-value of combined test from Monte Carlo}
#' @author Rui Cao
#' @examples
#'
#' set.seed(1)
#' p = 10
#' alpha_grid = sort(exp(seq(log(0.05), log(0.05/p), length.out=p*2)))
#' nsim = 1e4
#' null_dist <- TraitScan.MC(nsim = nsim, p = p, alpha_grid = alpha_grid)
#'
#' beta_hat = c(rnorm(8,0,0.1), 0.3, 0.3)
#' se_hat = rep(0.1, p)
#' null_cov = matrix(0.2, nrow = p, ncol = p)
#' diag(null_cov) = 1
#'
#' example <- TraitScan.summary.level(beta_hat, se_hat, null_cov, alpha_grid, null_dist)
#' @rdname summary_test
#' @export
#'

TraitScan.summary.level <- function(beta_hat, se_hat, null_cov, alpha_grid, null_dist){
  interest.data <- list("summary" = cbind(beta_hat, se_hat, beta_hat/se_hat))
  test.output <- combined.test(interest.data, null_cov, null_dist$null_dist_TC_prior, alpha_grid, null_dist$null_dist_TC_stat)
  p.value <- sum(null_dist$null_dist_combined_test <= test.output$test.statistics) / length(null_dist$null_dist_combined_test)
  if(p.value == 0) p.value <- paste0("<", 1/length(null_dist$null_dist_combined_test))
  return( list(select.traits.minp = test.output$select,
               select.traits.HC = test.output$select.HC,
               select.traits.TC = test.output$select.TC,
               tests.minp = test.output$minp.test,
               p.value = p.value) )

}
