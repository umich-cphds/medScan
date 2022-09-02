
`medScan` ....

### How to run

```r
remotes::install_github("umich-cphds/medScan")
```

Example

```r
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
```
