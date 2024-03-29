#' TraitScan extension to multiple SNPs
#'
#' In TraitScan, single SNP association with multiple traits can be extended to multiple SNPs (i.e. a genetic score) with multiple traits. By Pattee and Pan (2020), the multi-SNP association test statistic can be approximated by marginal SNP-wise statistics and an external reference panel.
#'
#' @param weight a numeric vector, the weights of SNPs in a genetic score. Can be obtained from public polygenic risk score or transcript imputation data.
#' @param beta_hat a numeric vector, estimated regression coefficients of traits.
#' @param se_hat a numeric vector, estimated standard errors of beta_hat.
#' @param ref_panel a genotype matrix with \code{length(weight)} columns. Rows are the samples from the reference panel and columns are the SNPs corresponding to \code{weight}.
#' @param N the sample size of the GWAS data (NOT the reference panel).
#' @param cont the trait is continuous (linear regression) or binary (logitic regression)
#' @param ratio_control2case the control-case ratio if the trait is binary.
#'
#' @return A z-statistic of the association between the genetic score and each trait.
#' @author Rui Cao, Tianzhong Yang
#' @examples
#'
#' set.seed(1)
#' p = 10    # number of traits
#' rho = 0.2
#' weight = 1:10/10
#' N = 1000
#' sigma = matrix(rho,nrow=p,ncol=p)
#' diag(sigma) = 1
#' X = MASS::mvrnorm(N,rep(0,p),sigma)
#' Y = rnorm(N)
#' beta_hat = se = rep(0,p)
#' X_s = 0
#' for(i in 1:p){
#'   temp = summary(lm(Y~X[,i]))
#'   beta_hat[i] = temp$coefficients[2,1]
#'   se[i] = temp$coefficients[2,2]
#'   X_s = X_s + weight[i] * X[,i]
#' }
#' model = lm(Y~X_s)
#' summary(model)$coefficients
#' TraitScan.multiSNPs(weight,beta_hat,se,ref_panel=X,N,cont = T)
#'
#' @references
#' Cao, Rui, et al. "Subset scanning for multi-trait analysis using GWAS summary statistics." Bioinformatics (2024): btad777.
#' Pattee, Jack, and Wei Pan. "Penalized regression and model selection methods for polygenic scores on summary statistics." PLoS computational biology 16.10 (2020): e1008271.
#' @rdname TraitScan_multiSNPs
#' @export
#'

TraitScan.multiSNPs <- function(weight,beta_hat,se,ref_panel,N,cont = T,ratio_control2case = NULL){
  z <- z_multi_snps(weight,beta_hat,se,ref_panel,N,cont,ratio_control2case)
  return(z)
}
