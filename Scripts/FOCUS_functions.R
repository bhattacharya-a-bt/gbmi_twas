estimate_cor = function(wmat,ldmat,intercept=F){
  wcov = t(wmat) %*% ldmat %*% wmat
  scale = diag(1/sqrt(diag(wcov)))
  wcor = scale %*% wcov %*% scale
  if (intercept){
    inter = scale %*% t(wmat) %*% ldmat
    return(list(wcor,
                inter))
  } else {
    return(list(wcor,NA))
  }
}

bayes_factor = function(zscores,
                        idx_set,
                        wcor,
                        prior_chisq = 40,
                        prb = 1e-3,
                        use_log = T){
  m = length(zscores)
  nc = length(idx_set)
  cur_chi2 = prior_chisq/nc
  cur_wcor = wcor[idx_set,idx_set]
  cur_zscors = zscores[idx_set]
  if (nc > 1){
    sss = svd(cur_wcor)
    cur_U = sss$u
    cur_EIG = sss$d
    rm(sss)
  } else {
    cur_U = 1
    cur_EIG = 1
  }
  scaled_chisq = cur_zscors^2
  
  cur_bf = .5 * -1*sum(log(1 + cur_chi2 %*% cur_EIG)) +
    .5 * sum((cur_chi2 / (1 + cur_chi2 %*% cur_EIG)) * scaled_chisq) +
    nc * log(prb) + (m-nc) * log(1-prb)
  if (use_log){
    return(cur_bf)
  } else {
    return(exp(cur_bf))
  }
  
}


get_resid = function(zscores, swld, wcor){
  m = nrow(wcor)
  p = ncol(swld)
  intercept = swld %*% rep(1,p)
  wcor_inv = MASS::ginv(as.matrix(wcor))
  rank = Matrix::rankMatrix(wcor_inv)[1]
  numer = t(intercept) %*% wcor_inv %*% zscores
  denom = t(intercept) %*% wcor_inv %*% intercept
  alpha = numer/denom
  resid = zscores -
    (t(intercept) * alpha)
  s2 = (resid %*% wcor_inv %*%
          t(resid))/(rank-1)
  inter_se = sqrt(s2/denom)
  inter_z = alpha/inter_se
  return(list(resid,
              inter_z))
}