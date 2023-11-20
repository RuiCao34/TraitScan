z_multi_snps_cont <- function(c,beta_hat,se,ref_panel,N){
  beta_hat = as.matrix(beta_hat)
  se = as.matrix(se)
  ref_panel = as.matrix(ref_panel)
  n = nrow(ref_panel)
  q = length(c)
  s2 = apply(ref_panel, MARGIN = 2, var)
  # 1/N XsTXs
  XsTXs = 0
  for(i in 1:q){
    for(j in 1:q){
      XsTXs = XsTXs + 1/n * c[i] * c[j] * sum(ref_panel[,i] * ref_panel[,j])
    }
  }
  
  # 1/N XsTy
  XsTy = sum(c * s2 * beta_hat)
  
  # 1/N yTy
  yTy = N * s2 * se^2 + s2 * beta_hat^2
  yTy = median(yTy)
  
  sigma_s_square_hat = N/(N-q) * (yTy - XsTy/XsTXs*XsTy)

  beta_s_hat = XsTy/XsTXs
  se_beta_s = sqrt(sigma_s_square_hat/(XsTXs*N))
  z = beta_s_hat/se_beta_s
  return(z)
}

z_multi_snps <- function(c,beta_hat,se,ref_panel,N,cont = T,ratio_control2case = NULL){
  if(cont == T){
    z = z_multi_snps_cont(c,beta_hat,se,ref_panel,N)
  }else{
    exp_neg_b0h = ratio_control2case
    beta_hat_linear = exp_neg_b0h / (1 + exp_neg_b0h)^2 * beta_hat
    se_linear = exp_neg_b0h / (1 + exp_neg_b0h)^2 * se
    z = z_multi_snps_cont(c,beta_hat_linear,se_linear,ref_panel,N)
  }
  return(z)
}