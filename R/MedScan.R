library("DACT")
library("HDMT")
library("locfdr")
library("qqman")


##### methods function #####
Sobel.est = function(z.alpha,z.beta){
  T.sobel = z.beta/sqrt(1+(z.beta/z.alpha)^2)
  p.sobel = 2*pnorm(-abs(T.sobel),mean = 0, sd=1, lower.tail = TRUE)
  return(p.sobel)
}

maxP.est = function(p.alpha,p.beta){
  return(pmax(p.alpha,p.beta))
}

JTcomp.est = function(z.alpha,z.beta){
  
  # Running time: 17.61s
  cat("The range of test statistics, from 0 to (default=10): ")
  # int<-as.numeric(readLines(con=stdin(),1))
  int <- max(c(abs(z.alpha*z.beta)/sd(z.alpha),abs(z.alpha*z.beta)/sd(z.beta)))+2
  
  B<-20000; pdf<-NULL
  for (i in 1:B){
    pdfi<-besselK(x=int*i/B, nu=0)
    pdf<-c(pdf, pdfi)
    flush.console()
  }
  
  myp<-function(cut){
    select<-(int*1:B/B)>cut
    pdf.sub<-pdf[select]
    pval<-sum(pdf.sub)/sum(pdf)
    return(pval)
  }
  
  MT_Comp<-function(a, b){
    ab<-a*b
    pp0<-sapply(abs(ab)/sqrt(1), myp)
    pp1<-sapply(abs(ab)/sd(a), myp)
    pp2<-sapply(abs(ab)/sd(b), myp)
    pp.comp<-pp1+pp2-pp0
    return(pp.comp)
  }
  
  p.JTcomp <- MT_Comp(z.alpha, z.beta)
  p.JTcomp[which(p.JTcomp<=0)] = min(p.JTcomp[which(p.JTcomp>0)])
  p.JTcomp[which(p.JTcomp>=1)] = 1-1e-8
  
  return(p.JTcomp)
  
}

p_value_underH1 = function(input_pvalues, alpha1, alpha2) {
  
  # input_pvalues contains two columns. 
  # The first column is p_1j (p for alpha=0). The second column is p_2j (p for beta=0).
  # alpha1: proportion of null alpha=0
  # alpha2: proportion of null beta=0
  
  pmax <- apply(input_pvalues,1,max)
  nmed <- length(pmax)
  efdr1 <- rep(0,nmed)
  
  nmed  <- nrow(input_pvalues)  
  cdf12 <- input_pvalues
  
  xx1 <- c(0,input_pvalues[order(input_pvalues[,1]),1])
  yy1 <- c(0,seq(1,nmed,by=1)/nmed)
  gfit1<- gcmlcm(xx1,yy1,type="lcm")
  xknots1 <- gfit1$x.knots[-1]
  Fknots1 <- cumsum(diff(gfit1$x.knots)*gfit1$slope.knots)
  
  xx2 <- c(0,input_pvalues[order(input_pvalues[,2]),2])
  yy2 <- c(0,seq(1,nmed,by=1)/nmed)
  gfit2<- gcmlcm(xx2,yy2,type="lcm")
  xknots2 <- gfit2$x.knots[-1]
  Fknots2 <- cumsum(diff(gfit2$x.knots)*gfit2$slope.knots)
  
  if (alpha1!=1) Fknots1 <- (Fknots1 - alpha1*xknots1)/(1-alpha1) else Fknots1 <- rep(0,length(xknots1))
  if (alpha2!=1) Fknots2 <- (Fknots2 - alpha2*xknots2)/(1-alpha2) else Fknots2 <- rep(0,length(xknots2))
  
  
  orderq1 <- pmax
  orderq2 <- pmax
  
  gcdf1 <- pmax
  gcdf2 <- pmax
  for (i in 1:length(xknots1)) {
    if (i==1) {
      gcdf1[orderq1<=xknots1[i]] <- (Fknots1[i]/xknots1[i])*orderq1[orderq1<=xknots1[i]] 
    } else {   
      if (sum(orderq1>xknots1[i-1] & orderq1<=xknots1[i])>0){
        print(i)
        temp <- orderq1[orderq1>xknots1[i-1] & orderq1<=xknots1[i]] 
        gcdf1[orderq1>xknots1[i-1] & orderq1<=xknots1[i]] <- Fknots1[i-1] + (Fknots1[i]-Fknots1[i-1])/(xknots1[i]-xknots1[i-1])*(temp-xknots1[i-1])
      }
    }
  }
  
  for (i in 1:length(xknots2)) {
    if (i==1) {
      gcdf2[orderq2<=xknots2[i]] <- (Fknots2[i]/xknots2[i])*orderq2[orderq2<=xknots2[i]] 
    } else {   
      if (sum(orderq2>xknots2[i-1] & orderq2<=xknots2[i])>0){
        temp <- orderq2[orderq2>xknots2[i-1] & orderq2<=xknots2[i]] 
        gcdf2[orderq2>xknots2[i-1] & orderq2<=xknots2[i]] <- Fknots2[i-1] + (Fknots2[i]-Fknots2[i-1])/(xknots2[i]-xknots2[i-1])*(temp-xknots2[i-1])
      } 
    }
  }
  
  
  gcdf1 <- ifelse(gcdf1>1,1,gcdf1)
  gcdf2 <- ifelse(gcdf2>1,1,gcdf2)
  
  cdf12[,1] <- gcdf1 # p_1j under H10
  cdf12[,2] <- gcdf2 # p_2j under H01
  
  return(cdf12 = cdf12)
  
}

HDMT.est = function(p.alpha,p.beta){
  input_pvalues = cbind(p.alpha,p.beta)
  pmax <- apply(input_pvalues,1,max)
  nullprop<- null_estimation(input_pvalues)
  c = nullprop$alpha10+nullprop$alpha01+nullprop$alpha00
  pi.01.est = nullprop$alpha01/c
  pi.10.est = nullprop$alpha10/c
  pi.00.est = nullprop$alpha00/c
  
  x <- tryCatch(p_value_underH1(input_pvalues = input_pvalues, alpha1 = nullprop$alpha1, alpha2 = nullprop$alpha2), error = function(e) e)
  if(!inherits(x, "error")){
    correction = p_value_underH1(input_pvalues = input_pvalues, alpha1 = nullprop$alpha1, alpha2 = nullprop$alpha2)
    p_1j_H10 = correction[,1]
    p_2j_H01 = correction[,2]
    p.HDMT = pi.10.est*pmax*p_1j_H10+pi.01.est*pmax*p_2j_H01+pi.00.est*pmax^2
  }else if(inherits(x, "error")){
    p.HDMT = pi.10.est*pmax+pi.01.est*pmax+pi.00.est*pmax^2
  }
  
  p.HDMT[which(p.HDMT<=0)] = min(p.HDMT[which(p.HDMT>0)])
  p.HDMT[which(p.HDMT>=1)] = 1-1e-8
  
  return(p.HDMT)
}

Sobelcomp.est = function(z.alpha,z.beta,p.alpha,p.beta){
  input_pvalues = cbind(p.alpha,p.beta)
  pmax <- apply(input_pvalues,1,max)
  nullprop<- null_estimation(input_pvalues)
  c = nullprop$alpha10+nullprop$alpha01+nullprop$alpha00
  pi.01.est = nullprop$alpha01/c
  pi.10.est = nullprop$alpha10/c
  pi.00.est = nullprop$alpha00/c
  
  T.sobel = z.beta/sqrt(1+(z.beta/z.alpha)^2)
  
  p.sobel.comp.01 = p.sobel.comp.10 = 2*pnorm(-abs(T.sobel),mean = 0, sd=1, lower.tail = TRUE)
  p.sobel.comp.00 = 2*pnorm(-abs(T.sobel),mean = 0, sd=0.5, lower.tail = TRUE)
  p.sobel.comp = (pi.01.est+pi.10.est)*p.sobel.comp.01+pi.00.est*p.sobel.comp.00
  
  return(p.sobel.comp)
}

DACT.est = function(p.alpha,p.beta){
  x <- tryCatch(DACT(p_a=p.alpha,p_b=p.beta,correction='JC'), error = function(e) e)
  if(!inherits(x, "error")){
    p.DACT = DACT(p_a=p.alpha,p_b=p.beta,correction='JC')
  }else if(inherits(x, "error")){
    warning("DACT fails")
    p.DACT = rep(NA,length(p.alpha))
  }
  
  return(p.DACT)
}

##### main function #####
medScan = function(z.alpha, z.beta, method){
  # z.alpha: the z-test statistic for alpha (exposure effect on the mediator)
  # z.beta:  the z-test statistic for beta (mediator effect on the outcome)
  # method: choose one from Sobel, MaxP, JT_comp, HDMT, DACT, Sobel_comp
  
  methods.all = c("Sobel","MaxP", "JT_comp", "HDMT", "DACT", "Sobel_comp")
  if(!(method %in% methods.all)){
    warning("invalid method choice")
  }
  
  p.alpha = pnorm(-abs(z.alpha),lower.tail = TRUE)*2
  p.beta = pnorm(-abs(z.beta),lower.tail = TRUE)*2
  
  # Sobel
  if(method == "Sobel"){ pvalues = Sobel.est(z.alpha,z.beta)}
  
  # MaxP
  if(method == "MaxP"){ pvalues = maxP.est(p.alpha,p.beta)}
  
  # JT-comp
  if(method == "JT_comp"){ pvalues = JTcomp.est(z.alpha,z.beta)}
  
  # HDMT
  if(method == "HDMT"){ pvalues = HDMT.est(p.alpha,p.beta)}
  
  # DACT
  if(method == "DACT"){ pvalues = DACT.est(p.alpha,p.beta)}
  
  # Sobel_comp
  if(method == "Sobel_comp"){ pvalues = Sobelcomp.est(z.alpha,z.beta,p.alpha,p.beta)}
  
  # proportion estimation
  input_pvalues = cbind(p.alpha,p.beta)
  pmax <- apply(input_pvalues,1,max)
  nullprop<- null_estimation(input_pvalues)
  c = nullprop$alpha10+nullprop$alpha01+nullprop$alpha00
  pi.01.est = nullprop$alpha01/c
  pi.10.est = nullprop$alpha10/c
  pi.00.est = nullprop$alpha00/c
  pi = data.frame("pi00"=pi.00.est, "pi01"=pi.01.est,"pi10"=pi.10.est)
  
  return(list(pvalues = pvalues,
              pi = pi))
  
}


###########################
######## example ##########
###########################

# simulate data under the mixture null
n=10000 
u = runif(n,0,1)
z.alpha = z.beta = rep(NA,0)
pi00 = 0.98
pi10 = 0.01
pi01 = 0.01
for(i in 1:n){
  if(u[i]<=pi00){
    z.alpha[i] = rnorm(1, 0, 1)
    z.beta[i] = rnorm(1, 0, 1)
  } else if (u[i]<= pi00+pi10){
    z.alpha[i] = rnorm(1, 1, 1)
    z.beta[i] = rnorm(1, 0, 1)
  } else {
    z.alpha[i] = rnorm(1, 0, 1)
    z.beta[i] = rnorm(1, 1, 1)
  }
}

# obtain p-values

obj = medScan(z.alpha = z.alpha, z.beta = z.beta, method = "Sobel")
qqman::qq(obj$pvalues,  xlim = c(0,4), ylim = c(0,4), main = "Sobel")

obj = medScan(z.alpha = z.alpha, z.beta = z.beta, method = "MaxP")
qqman::qq(obj$pvalues,  xlim = c(0,4), ylim = c(0,4), main = "MaxP")

obj = medScan(z.alpha = z.alpha, z.beta = z.beta, method = "HDMT")
qqman::qq(obj$pvalues, xlim = c(0,4), ylim = c(0,4), main="HDMT")

obj = medScan(z.alpha = z.alpha, z.beta = z.beta, method = "Sobel_comp")
qqman::qq(obj$pvalues,  xlim = c(0,4), ylim = c(0,4), main = "Sobel-comp")

obj = medScan(z.alpha = z.alpha, z.beta = z.beta, method = "JT_comp")
qqman::qq(obj$pvalues, xlim = c(0,4), ylim = c(0,4), main="JT-comp")

obj = medScan(z.alpha = z.alpha, z.beta = z.beta, method = "DACT")
qqman::qq(obj$pvalues,  xlim = c(0,4), ylim = c(0,4), main="DACT")
