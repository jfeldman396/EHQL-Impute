
#' Title
#'
#' @param j s
#' @param Z s
#' @param Y s
#' @param R s
#' @param Rlevels s
#' @param alpha s
#' @param Lambda s
#' @param eta s
#' @param Sigma.diag s
#' @param plugin.marginal s
#' @param is_cat_bin s
#' @param have_aux s
#' @param aux_quantiles s
#' @param aux_bins s
#' @param Z_past s
#'
#' @import truncnorm
#' @return s
#' @export

resample_Z_MA<-function(j,Z,Y,R,Rlevels,alpha, Lambda, eta, Sigma.diag,
                        plugin.marginal,is_cat_bin,have_aux,aux_quantiles,aux_bins,
                        Z_past){


  muj = alpha[j] +Lambda[j,]%*%t(eta)
  sdj = sqrt(Sigma.diag[j])
  margsdj = sqrt((tcrossprod(Lambda) + diag(Sigma.diag))[j,j])
  Zj = Z[,j]
  Z_pastj = Z_past[,j]
  Yj = Y[,j]
  nobs = nrow(Y)
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

        ir <- (1:nobs)[R[, j] == r & !is.na(R[, j])]
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
  ir <- (1:nobs)[is.na(R[, j])]
  Zj[ir]<- rnorm(length(ir), muj[ir], sdj)

  return(Zj)


}
