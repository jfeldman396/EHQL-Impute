
<!-- README.md is generated from README.Rmd. Please edit that file -->

# EHQLImpute

<!-- badges: start -->
<!-- badges: end -->

The EHQL Gaussian copula is a powerful tool for imputing nonignorable
missing data with arbitrary marginal distributions. In what follows we
will walk you through a toy example

## Installation

You can install the development version of EHQLImpute from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("jfeldman396/EHQL-Impute")
```

``` r
library(EHQLImpute)
```

## Simulating toy data

We will start by simulating data with non-ignorable missingness from a
Gaussian copula with arbitrary marginal distribution functions. We begin
by simulating $(\boldsymbol Y,\boldsymbol R)$, and removing values
$\boldsymbol Y_{ij}$ for $R_{ij}>0$.

``` r
  library(LaplacesDemon)
  # simulate data generating copula correlation
  C_prior = diag(1,10)
  C_0= cov2cor(LaplacesDemon::rinvwishart(15, C_prior))
  
  # create intercept vector to encode marginal missingness: we'll use 50\% marginal missingness in this example. Recall that \alpha will be 0 for indices corresponding to Y. We'll have non-ignorable missingness for each study variable (p = 5)
  
  alpha  = c(rep(0,5),rep(qnorm(0.5),5))
  
  # Generate latent variables: we'll simulate 1000 observations
  library(MASS)
  Z <- MASS::mvrnorm(1000,mu = alpha, Sigma = C_0)
  
  # Create R
  
  R <- t(apply(Z[,6:10],1,function(z) ifelse(z>0,1,0)))
  
  #look at marginal incidence of missingness
  print("Missingness in Each Y_j:")
#> [1] "Missingness in Each Y_j:"
  print(colMeans(R))
#> [1] 0.524 0.494 0.535 0.517 0.507

  # Now we'll create the study variables. Here we use gamma, t, and beta marginals 
  
      Y_raw =NULL
      ind = 1
      
      for(j in 1:5){
        # the sequence of marginals goes gamma, t, beta
        if(ind == 1){
          Y_j= qgamma(pnorm(Z[,j]),1,1)
        
          
        }else if(ind == 2){
          Y_j = qt(pnorm(Z[,j]),5,2)
       
          
        }else if(ind == 3){
          Y_j<-qbeta(pnorm(Z[,j]),1,2)
        
        }
        if(ind == 3){
          ind = 1
        }else{
          ind = ind+1
        }
        Y_raw = cbind(Y_raw,Y_j)
      }
      
      # remove values based on R = 1
      Y_raw[which(R ==1, arr.ind = T)] = NA
      
      #join Y and R to get observed data
      Y_obs=cbind(Y_raw,R)
      print("Summary of (Y^obs, R)")
#> [1] "Summary of (Y^obs, R)"
      print(summary(Y_obs))
#>       Y_j              Y_j               Y_j              Y_j        
#>  Min.   :0.0001   Min.   :-0.7307   Min.   :0.0014   Min.   :0.0051  
#>  1st Qu.:0.1904   1st Qu.: 1.4349   1st Qu.:0.2534   1st Qu.:0.3093  
#>  Median :0.4680   Median : 2.1258   Median :0.4149   Median :0.7714  
#>  Mean   :0.7183   Mean   : 2.4168   Mean   :0.4312   Mean   :1.0853  
#>  3rd Qu.:1.0598   3rd Qu.: 3.1627   3rd Qu.:0.6005   3rd Qu.:1.5769  
#>  Max.   :5.6446   Max.   : 8.4221   Max.   :0.9458   Max.   :5.6819  
#>  NA's   :524      NA's   :494       NA's   :535      NA's   :517     
#>       Y_j               V6              V7              V8       
#>  Min.   :-1.243   Min.   :0.000   Min.   :0.000   Min.   :0.000  
#>  1st Qu.: 1.048   1st Qu.:0.000   1st Qu.:0.000   1st Qu.:0.000  
#>  Median : 1.747   Median :1.000   Median :0.000   Median :1.000  
#>  Mean   : 1.889   Mean   :0.524   Mean   :0.494   Mean   :0.535  
#>  3rd Qu.: 2.560   3rd Qu.:1.000   3rd Qu.:1.000   3rd Qu.:1.000  
#>  Max.   : 7.084   Max.   :1.000   Max.   :1.000   Max.   :1.000  
#>  NA's   :507                                                     
#>        V9             V10       
#>  Min.   :0.000   Min.   :0.000  
#>  1st Qu.:0.000   1st Qu.:0.000  
#>  Median :1.000   Median :1.000  
#>  Mean   :0.517   Mean   :0.507  
#>  3rd Qu.:1.000   3rd Qu.:1.000  
#>  Max.   :1.000   Max.   :1.000  
#> 
```

## Estimating the EHQL copula

Now we’ll estimate the EHQL copula. To do so requires several pieces of
information specified by the user:

- $\texttt{YR}$: This is the combined observed study variables and
  missingness indicators. It is given by $\texttt{Yobs}$ above.
- $\texttt{ncolY}$: The number of study variables. In this example we
  have $\texttt{ncolY} = 5$:
- $\texttt{ncolR}$: The number of study variables modeled as
  non-ignorable. In this example we have $\texttt{ncolR} = 5$
- aux_quantiles: auxiliary quantiles assumed known for each study
  variable. This is a list of length $\texttt{ncolY}$, where the
  auxquantiles\[\[j\]\] = $(0,\{\tau_{k_{j}}\}_{k=1}^{q},...,1\}$,
  i.e. there are $q_{j}+2$ auxiliary quantiles assumed known for each
  variable. If $\texttt{aux_quantiles[[j]]}$ is null, empirical deciles
  will be used. Note the quantiles for each variable do not need to be
  unique.
- aux_infos: This is a list of length $\texttt{ncolY}$ where
  aux_infos\[\[j\]\] =
  $(F_{j}^{-1}(0),\{F_{j}^{-1}(\tau_{k_{j}})\}_{k=1}^{q},...,F_{j}^{-1}(1)\}$
- $\texttt{MA}$: vector of length ncolY indicating whether or not to
  compute the margin adjustment. This is strongly recommended for all
  levels of auxiliary information, and especially when auxiliary
  information is extremely sparse, as it propagates uncertainty about
  each $F_{j}$ which is beneficial for prediction and imputation
- $\texttt{nImps}$: The number of completed data sets to create
- $\texttt{nsamp}$: Number of iterations for the MCMC
- $\texttt{burn}$: Burn-in iterations for the MCMC

``` r
YR = Y_obs
ncolY = 5; ncolR = 5
MA = rep(T,5)
aux_quantiles = vector('list',5)
aux_infos = vector('list',5)

for(j in 1:5){
  aux_quantiles[[j]] = c(0,0.5,1) #assume access to the median
  if(j ==1 | j == 4){
    aux_infos[[j]] = qgamma(aux_quantiles[[j]],1,1)
  }else if(j ==2 |j == 5){
    aux_infos[[j]] = qt(aux_quantiles[[j]],5,2)
  }else{
    aux_infos[[j]] = qbeta(aux_quantiles[[j]],1,2)
  }
  
}

imps<- EHQLImpute(YR,
                  ncolY,
                  ncolR,
                  aux_quantiles = aux_quantiles,
                  aux_infos = aux_infos,
                  MA = MA,
                  nImps = 20,
                  nsamp = 5000,
                  burn =4000)
#> Sampling is: 20 percent doneNULL
#> Sampling is: 40 percent doneNULL
#> Sampling is: 60 percent doneNULL
#> Sampling is: 80 percent doneNULL
#> Sampling is: 100 percent doneNULL
```

## Analyze results

Now, we can look at the imputations

``` r
summary(imps$YImpute[[2]])
#>       Y_j                Y_j               Y_j                Y_j         
#>  Min.   :0.000064   Min.   :-0.7307   Min.   :-0.01336   Min.   :0.00189  
#>  1st Qu.:0.254413   1st Qu.: 1.3588   1st Qu.: 0.10563   1st Qu.:0.25395  
#>  Median :0.686923   Median : 2.1461   Median : 0.29950   Median :0.68230  
#>  Mean   :1.029893   Mean   : 2.3823   Mean   : 0.32601   Mean   :1.03111  
#>  3rd Qu.:1.438245   3rd Qu.: 3.1300   3rd Qu.: 0.50142   3rd Qu.:1.40615  
#>  Max.   :6.194027   Max.   : 8.4221   Max.   : 0.94584   Max.   :6.14363  
#>       Y_j        
#>  Min.   :-1.243  
#>  1st Qu.: 1.232  
#>  Median : 2.065  
#>  Mean   : 2.343  
#>  3rd Qu.: 3.169  
#>  Max.   : 7.702

plot(density(imps$YImpute[[2]][,3], bw = .1),col = "blue", lwd = 3, main = "Imputed vs. Truth vs. Observed")
lines(density(rbeta(1000,1,2),bw = .1),lwd = 3, lty = 2)
lines(density(Y_raw[,3],na.rm = T, bw = .1), lwd = 3, col = "gray")
legend("topright", c("Imputed", "Truth", "Observed"), col= c("blue", "black", "gray"), lwd = 3)
```

<img src="man/figures/README-unnamed-chunk-4-1.png" width="100%" />

We see that the imputations correct bias in the observed data.

### Analyze margin Adjustment

``` r
# 


  for( j in 1:3){

    if(j == 2){
    plot(imps$x_ma[[j]],imps$MAs[[j]][1,], type = 'n', xlab = expression(Y[2]),
         ylab = expression(paste(P(Y[2]<=y[2]))),
         main =paste("Margin Adjustment Posterior Samples for Y_2"),
         ylim = c(0,1),
         lwd = 3, 
         cex = 2)
    }else if(j==1){
      plot(imps$x_ma[[j]],imps$MAs[[j]][1,], type = 'n', xlab = expression(Y[2]),
         ylab = expression(paste(P(Y[1]<=y[1]))),
         main =paste("Margin Adjustment Posterior Samples for Y_1"),
         ylim = c(0,1),
         lwd = 3, 
         cex = 2)
    }else{
      plot(imps$x_ma[[j]],imps$MAs[[j]][1,], type = 'n', xlab = expression(Y[2]),
         ylab = expression(paste(P(Y[3]<=y[3]))),
         main =paste("Margin Adjustment Posterior Samples for Y_3"),
         ylim = c(0,1),
         lwd = 3, 
         cex = 2)
    }
    sapply(1:1000, function(i)
      lines(imps$x_ma[[j]], imps$MAs[[j]][i,], type = "b", col = "gray", lwd = 2))
    if(j == 1){
      points(imps$x_ma[[j]],pgamma(imps$x_ma[[j]],1,1), pch = 16, cex =4)
    }else if(j == 2){
      points(imps$x_ma[[j]],pt(imps$x_ma[[j]],5,2), pch = 16, cex = 4)
    }else{
      points(imps$x_ma[[j]],pbeta(imps$x_ma[[j]],1,2), pch = 16, cex = 3)
    }

  }
```

<img src="man/figures/README-unnamed-chunk-5-1.png" width="100%" /><img src="man/figures/README-unnamed-chunk-5-2.png" width="100%" /><img src="man/figures/README-unnamed-chunk-5-3.png" width="100%" />
