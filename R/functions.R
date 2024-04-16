#' Title
#'
#' @param YR Observed data matrix. It should be partititioned into the study variables Y followed by missingness indicators R.
#' @param ncolY Number of study variables
#' @param ncolR Number of study variables modeled as noningorable
#' @param aux_quantiles List of length ncolY containing auxiliary quantiles assumed known for each study variable. aux_quantiles[[j]] should be a vector containg 0, intermediate quantiles, and 1.
#' @param aux_infos List of length ncolY containing the known auxiliary values of the study variables. aux_infos[[j]] should be a vector containing \eqn{(F_j^{-1}(0),\dots,F_{j}^{-1}(\tau),\dots,F_{j}^{-1}(1))}. In addition aux_infos[[j]] should be the same length as aux_quantiles[[j]]
#' @param MA Vector logicals of length ncolY indicating whether or not to estimate the margin adjustment at intermediate quantiles using the margin adjustment. Default is True for each study variable.
#' @param nImps Number of completed data sets to produce
#' @param namp Number of iterations for the MCMC
#' @param burn Number of Burn-in samples for the MCMC
#'
#' @return C.post,alpha.post, and MAs, which are the nsamp-burn posterior samples of the copula correlation, latent intercept vector, and margin adjustment (list or length ncolY, MAs[[j]] computes the margin adjustment at x_ma[[j]]), respectively. In addition, YImpute is a list of length nImps containing imputed data sets.
#' @export

EHQLImpute<- function(YR,
                      ncolY,
                      ncolR,
                      aux_quantiles,
                      aux_infos,
                      MA = rep(T, ncolY),
                      nImps,
                      nsamp,
                      burn){
  
  
  F_invs = vector('list',ncolY) #interpolated cdf
  MAs = vector('list', ncolY)
  aux_bins = vector('list',ncolY) #get aux bin membership
  Completed_Data = vector('list',nImps)# store completed data
  
  imp_every = floor(nsamp/nImps)
  
  ind =1
  Y_binned = NULL
  for(j in 1:ncolY){
    Y_j = YR[,j]
    if(is.null(aux_quantiles[[j]])){
      
      # if there is no auxiliary information, use empirical deciles quantile for EQL
      auxq = seq(0,1,by = 0.10)
      
      if(MA[j]){
        aux_info=sort(unique(c(seq(min(Y_j,na.rm =T),max(Y_j,na.rm = T), length.out = 15),
                               quantile(YR[,j],
                                        seq(0,1,by = 0.10),
                                        type = 1))))
        #Bin data for EQL/EHQL
        Y_j_binned = aux_info[findInterval(Y_j,aux_info, left.open = T)+1]
      }else{
        aux_info =   quantile(YR[,j],auxq, type = 1)
        #Bin data for EQL/EHQL
        Y_j_binned = aux_info[findInterval(Y_j,aux_info, left.open = T)+1]
      }
      #get empirical quantiles
      
      
      aux_bins[[j]] = findInterval(Y_j,quantile(YR[,j],auxq, type = 1), left.open = T)+1 #get auxiliary bin membership
      aux_quantiles[[j]] = auxq # save aux quantiles
      
      
      #save interpolated cdf
      F_invs[[j]] = splinefun(seq(0,1,by = 0.5),auxq) #save interpolated ecdf
      
    }
    else{
      auxq = aux_quantiles[[j]]
      a_info = aux_infos[[j]]
      if(MA[j]){
        aux_info=sort(unique(c(seq(min(Y_j,na.rm = T),max(Y_j,na.rm = T), length.out = 15),a_info)))
        #Bin data for EQL/EHQL
        Y_j_binned = aux_info[findInterval(Y_j,aux_info, left.open = T)+1]
      }else{
        aux_info =   quantile(YR[,j],auxq, type = 1)
        #Bin data for EQL/EHQL
        Y_j_binned = aux_info[findInterval(Y_j,aux_info, left.open = T)+1]
      }
      #get empirical quantiles
      
      
      aux_bins[[j]] = findInterval(Y_j,a_info, left.open = T)+1 #get auxiliary bin membership
    }
    Y_binned = cbind(Y_binned, Y_j_binned)
    
  }
  
  #create binned data
  Y <- as.matrix(cbind(Y_binned,YR[,((ncolY+1):ncol(YR))]))
  n = nobs = dim(Y)[1]
  p <- dim(Y)[2]
  plugin.marginal = rep(F,p)
  
  
  
  #initiate ranks and latent data
  R <- NULL
  for (j in 1:p) {
    
    R <- cbind(R, match(Y[, j], sort(unique(Y[, j]))))
    
  }
  Rlevels <- apply(R, 2, max, na.rm = TRUE)
  
  Ranks <- apply(Y, 2, rank, ties.method = "max", na.last = "keep")
  N <- apply(!is.na(Ranks), 2, sum)
  U <- t(t(Ranks)/(N + 1))
  Z = NULL
  nobs = n
  for(j in 1:p){
    if(j <= dim){
      Z = cbind(Z, U[,j])
    }else{
      Zj = array(0,n)
      Zj[which(Y[,j] == 1)] = rtruncnorm(sum(Y[,j] == 1),a = 0)
      Zj[which(Y[,j] == 0)] = rtruncnorm(sum(Y[,j] == 0),b = 0)
      Z = cbind(Z,Zj)
    }
    
  }
  
  Zfill <- matrix(rnorm(n * p), n, p)
  Z[is.na(Y)] <- Zfill[is.na(Y)]
  Z[is.na(Z)] =  Zfill[is.na(Z)]
  
  
  # initialize covariance matrix
  
  k.star = floor(p/2)
  #
  # Initialize Factor Model Parameters: these are
  
  # Local df:
  nu = 3;
  
  # Global shrinkage:
  a1 = 2; a2 = 3;
  
  # Error variance:
  a.sigma = 1; b.sigma = 0.3
  
  
  # Use SVD for simplicity
  svd0 = svd(Z); # plot(1-cumsum(svd0$d^2)/sum(svd0$d^2)); #plot(svd0$d[-1]/svd0$d[-min(p,n)])
  
  # Factor loadings (p x k.star):
  Lambda = svd0$v[1:p,1:k.star]
  
  # Factors (n x k.star):
  eta = svd0$u[,1:k.star] #y%*%Lambda
  
  # Diagonal error variance (p-dimensional vector):
  Sigma.diag = rep(1,p)
  #
  
  #Local and global precision parameters:
  phi.jh = 1/Lambda^2; tau.h = delta.h = rep(1, k.star)
  
  # intercept
  alpha = rep(0,p)
  C.post= array(dim = c(p, p, floor(nsamp/odens)- burn))
  
  alpha.post<- array(dim = c(p, floor(nsamp/odens)- burn))
  have_aux = rep(T,ncolY)
  
  #for eql
  is_cat_bin = c(rep(0,ncolY), rep(1,ncolR))
  
  
  #for imputations
  
  m = 1
  for (ns in 1:nsamp){
    
    
    
    
    #Step 1: Sample the factor loadings
    
    cp.eta = crossprod(eta)
    for(j in 1:p){
      chQj = chol(diag(phi.jh[j,]*tau.h, k.star) + cp.eta/Sigma.diag[j])
      lj = crossprod(eta, Z[,j] - alpha[j])/Sigma.diag[j]
      Lambda[j,] = backsolve(chQj,forwardsolve(t(chQj), lj) + rnorm(k.star))
    }
    
    
    # Step 2: sample the error variances
    eps = Z - sweep(tcrossprod(eta, Lambda),2,alpha,'+')
    
    Sigma.diag = apply(eps, 2, function(x) 1/rgamma(n = 1, shape = a.sigma + n/2,
                                                    rate = b.sigma + 1/2*sum(x^2)))
    
    # Step 3: sample the factors
    chQeta = chol(diag(k.star) + crossprod(Lambda, diag(1/Sigma.diag))%*%Lambda)
    leta = tcrossprod(crossprod(Lambda, diag(1/Sigma.diag)), sweep(Z,2,alpha,'-'))
    eta = t(backsolve(chQeta,forwardsolve(t(chQeta), leta) + rnorm(n*k.star))) #for(i in 1:n) eta[i,]= backsolve(chQeta,forwardsolve(t(chQeta), leta[,i]) + rnorm(k.star))
    
    # Step 4: sample phi.jh
    phi.jh = matrix(rgamma(n = p*k.star, shape = (nu + 1)/2,
                           rate = (nu + Lambda^2*matrix(rep(tau.h, each = p), nr = p))/2), nr = p) #for(h in 1:k.star){for(j in 1:p) phi.jh[j,h] = rgamma(n = 1, shape = (nu + 1)/2, rate = (nu + Lambda[j,h]^2*tau.h[h])/2)
    
    # Step 5: sample tau.h via delta.h
    delta.h = sampleMGP(theta.jh = sqrt(phi.jh)*Lambda, delta.h = delta.h, a1 = a1, a2 = a2)
    tau.h = cumprod(delta.h)
    
    # Step 6: Sample alpha
    
    eps = Z - tcrossprod(eta, Lambda)
    for(j in (ncolY +1):p){
      lalpha = sum(eps[,j])/Sigma.diag[j]
      Qalpha = (n/Sigma.diag[j] + 1)^-1
      alpha[j] = rnorm(1,Qalpha*lalpha,sqrt(Qalpha))
    }
    
    
    if(ns>= 1){
      Z_past = Z
      
      for(j in 1:p){
        
        Zj = resample_Z_MA(j,Z,Y,R,Rlevels,alpha, Lambda, eta, Sigma.diag,
                           plugin.marginal,is_cat_bin,have_aux,aux_quantiles,aux_bins,
                           nobs = n,
                           Z_past)
        Z[,j] = Zj
        
      }
      
    }
    
    if(any(is.na(Z))){
      break
    }
    if( (ns%%100) == 0) {
      
      print(cat(round(100 * ns/nsamp), "percent done"))
      
    }
    if (ns%%odens == 0 &ns > burn ) {
      S<-tcrossprod(Lambda) + diag(Sigma.diag)
      
      
      C <- S/(sqrt(diag(S)) %*% t(sqrt(diag(S))))
      
      C.post[, , ns - burn ] = C
      alpha.post[,ns-burn] = alpha/sqrt(diag(S))
      
      if((ns-burn)%%imp_every == 0){ # imputation
        
        
        Y_copy = YR[,1:ncolY]
        for(j in 1:ncolY){
          x_ma =apply(Y[,1:ncolY],2, function(x) sort(unique(x)))
          
          
          #cut points
          cuts_MA<- lapply(1:ncolY, function(x)
            sapply(x_ma[[x]], function(y)
              max(Z[which(Y[,x] == y),x])))
          
          # compute the margin adjusmtent
          MA_j<- sapply(1:ncolY,function(j) pnorm(cuts_MA[[j]],sd = sqrt(S[j,j])))
          for(j in 1:ncolY){
            MAs[[j]] = rbind(MAs[[j]], MA_j[[j]])# save margin adjustment samples
          }
          lb_j = aux_infos[[j]][1]
          ub_j = aux_infos[[j]][length(aux_infos[[j]])]
          Fn_inv = sapply(1:ncolY, function(j)
            splinefun(c(0,MA_j[[j]],1),c(lb_j,x_ma[[j]],ub_j)))
          
          na_inds = which(is.na(YR[,j]))
          
          Y_mis<- Fn_inv[[j]](pnorm(Z[na_inds,j], sd = sqrt(S[j,j])))
          
          
          
          Y_copy[na_inds,j] = Y_mis
          
          
          
          
        }
        
        Completed_Data[[m]] = Y_copy
        m = m+1
        
      }
      
      
    }
  }
  
  return(list(
    C.post = C.post,
    alpha.post = alpha.post,
    MAs = MAs,
    x_ma,
    YImpute = Completed_Data
  ))
}


resample_Z_MA<-function(j,Z,Y,R,Rlevels,alpha, Lambda, eta, Sigma.diag,
                        plugin.marginal,is_cat_bin,have_aux,aux_quantiles,aux_bins,
                        nobs = n,
                        Z_past){
  
  
  muj = alpha[j] +Lambda[j,]%*%t(eta)
  sdj = sqrt(Sigma.diag[j])
  margsdj = sqrt((tcrossprod(Lambda) + diag(Sigma.diag))[j,j])
  Zj = Z[,j]
  Z_pastj = Z_past[,j]
  Yj = Y[,j]
  bounds = NULL
  if(is_cat_bin[j] == 0){ #if j is numeric
    if(have_aux[j]){
      aux_quantile = aux_quantiles[[j]]
      aux_bin = aux_bins[[j]]
      if(!plugin.marginal[j]){
        for(r in 1:Rlevels[j]){
          ir_obs <- (1:nobs)[R[(1:nobs), j] == r & !is.na(R[(1:nobs), j])]
          
          #find aux interval
          aux_lb =suppressWarnings(qnorm(aux_quantile[aux_bin[ir_obs]-1], 0, margsdj))
          aux_ub = suppressWarnings(qnorm(aux_quantile[aux_bin[ir_obs]],0, margsdj))    #Y_ij obs has more information about the binning
          
          
          #fix upper and lower bounds based on auxiliary information and adjacent interval
          lb <- pmax(suppressWarnings(max(Zj[R[, j] == r -
                                               1], na.rm = TRUE)),aux_lb)
          ub <- pmin(suppressWarnings(min(Zj[R[, j] == r+1], na.rm = TRUE)),aux_ub)
          
          
          
          if(length(lb) == 0 | length(ub) == 0){
            
            Zj[ir_obs] = Z_pastj[ir_obs]
          }
          else if(any(lb>ub)){
            
            lb <- suppressWarnings(max(Zj[R[, j] == r -
                                            1], na.rm = TRUE))
            ub <- suppressWarnings(min(Zj[R[, j] == r+1], na.rm = TRUE))
            
            Zj[ir_obs] =truncnorm::rtruncnorm(length(ir_obs), a = lb,b = ub,muj[ir_obs], sdj)
            
          }else{
            Zj[ir_obs] =truncnorm::rtruncnorm(length(ir_obs), a = lb,b = ub,muj[ir_obs], sdj)
            
            if(any(is.na(Zj))){
              Zj[ir_obs] = Z_pastj[ir_obs]
            }
            
          }
          
          
          
        }
      }
    }
    else{# if there's no auxiliary information on Y_j, just do rank likelihood sampling
      for(r in 1:Rlevels[j]){
        
        ir <- (1:n)[R[, j] == r & !is.na(R[, j])]
        lb <- suppressWarnings(max(Zj[R[, j] == r -
                                        1], na.rm = TRUE))
        ub <- suppressWarnings(min(Zj[R[, j] == r +
                                        1], na.rm = TRUE))
        
        Zj[ir] = truncnorm::rtruncnorm(length(ir), a = lb,b = ub,muj[ir],sdj)
      }
      
    }
  }
  else{ #if j is a categorical sub-level
    
    Zj[which(Y[,j] == 1)] =   rtruncnorm(n = sum(Y[,j] == 1, na.rm = T),a =0 ,mean = muj[which(Y[,j] == 1)],sd = sdj)
    Zj[which(Y[,j] == 0)] = rtruncnorm(n = sum(Y[,j] == 0, na.rm = T),b = 0,mean = muj[which(Y[,j] == 0)],sd = sdj)
    
  }
  
  #impute
  ir <- (1:n)[is.na(R[, j])]
  Zj[ir]<- rnorm(length(ir), muj[ir], sdj)
  
  return(Zj)
  
  
}

sampleMGP<-function(theta.jh, delta.h, a1 = 3, a2 = 5){
  # Store the dimensions locally
  p = nrow(theta.jh); K = ncol(theta.jh)
  
  # Sum over the (squared) replicates:
  sum.theta.l = colSums(theta.jh^2)
  
  # h = 1 case is separate:
  tau.not.1 = cumprod(delta.h)/delta.h[1]
  delta.h[1] = rgamma(n = 1, shape = a1 + p*K/2,
                      rate = 1 + 1/2*sum(tau.not.1*sum.theta.l))
  # h > 1:
  if(K > 1){for(h in 2:K){
    tau.not.h = cumprod(delta.h)/delta.h[h]
    delta.h[h] = rgamma(n = 1, shape = a2 + p/2*(K - h + 1),
                        rate = 1 + 1/2*sum(tau.not.h[h:K]*sum.theta.l[h:K]))
  }}
  delta.h #list(tau.h = cumprod(delta.h), delta.h = delta.h)
}

