#' Monte Carlo Simulation
#'
#' Monte Carlo Simulation for the null distributions of TC and combined tests
#'
#' @param nsim int, number of MC iterations.
#' @param p int, number of traits.
#' @param alpha_grid numeric vector less than 1: the alpha grid for p value thresholds. Default as a sequence from the Bonferroni-significant threshold to 0.05.
#'
#' @return A list with the null distributions of TC priority function, TC statistic, and combimed test statistic
#' @note Use this function to simulate the null distributions before run tests \code{TraitScan.summary.level} or \code{TraitScan.individual.level}
#' @author Rui Cao, Tianzhong Yang
#' @examples
#'
#' p = 10    # number of traits
#' alpha_grid = sort(exp(seq(log(0.05), log(0.05/p), length.out=p*2)))    # alpha grid in HC statistic
#' nsim = 1e4    # number of MC iterations
#' null_dist <- TraitScan.MC(nsim = nsim, p = p, alpha_grid = alpha_grid)
#' @references
#' Cao, Rui, et al. "Subset scanning for multi-trait analysis using GWAS summary statistics." medRxiv (2023): 2023-07.
#' @rdname MC
#' @export
#'
TraitScan.MC <- function(nsim = 1e4, p, alpha_grid = sort(exp(seq(log(0.05), log(0.05/p), length.out=p*2)))){
  return( null.dist.Monte.Carlo(N=nsim, p, alpha_grid,para=F) )
}

