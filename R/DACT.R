
DACT = function(p_a,p_b,correction=NULL){
  Z_a = stats::qnorm(p_a,lower.tail = F)
  Z_b = stats::qnorm(p_b,lower.tail = F)
  pi0a = nonnullPropEst(Z_a,0,1)
  pi0b = nonnullPropEst(Z_b,0,1)
  #pi0a = locfdr::locfdr(Z_a,nulltype = 0)$fp0[5,3]
  #pi0b = locfdr::locfdr(Z_b,nulltype = 0)$fp0[5,3]
  if(pi0a > 1){
    pi0a = 1
  }
  if(pi0b >1){
    pi0b = 1
  }
  p.mat = cbind(p_a,p_b)
  p3 = (apply(p.mat,1,max))^2
  wg1 = pi0a*(1-pi0b)
  wg2 = (1-pi0a)*pi0b
  wg3 = pi0a*pi0b
  wg.sum = wg1 + wg2 + wg3
  wg.std = c(wg1,wg2,wg3)/wg.sum
  p_dact = wg.std[1]*p_a + wg.std[2]*p_b + wg.std[3]*p3
  if(correction == "Efron"){
  p_dact = EfronCorrect(p_dact)
  }
  if(correction == "JC"){
    p_dact = JCCorrect(p_dact)
  }
  return(p_dact)
}
