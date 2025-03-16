


# SS: stratification   +  SoF    (stratification+SoS  see package DRMR)

# parametric fitting using specific change-points (usually determined by SSS)


SS<-function(       Z ,
                    X ,
                    Y ,
                    pos , # the user specific change-point position vector for parametric inference
                    Ns_used=100, # the number of strata
                    x_baseline_used = NA ,  # the basline exposure level; default is mean(X)
                    precision =100    # the precision number used in the exposure quantile region; default is 100
){


  RES<-list()


  tpoints<- quantile(X, seq(0,1, length=precision+1)    ) # tpoints are used for weight function and its integration

  ### DR stratification with strata number = Ns_used ---------------------------------------------------------------------
  dat_<-data.frame(  Z=Z, X=X, Y=Y )   # X:= alcohol   with outcome SBP
  rdat<-Stratify(dat_, Ns = Ns_used) #DR stratification


  ### Scalar-on-Function regression-----------------------------------------------------------------------------------------

  ### stratum-specific results
  ZXest<-c()  #stratum-specific IV-exposure association
  ZYest<-c()  #stratum-specific IV-outcome association
  SEs_ZY<-c() #the s.e. of IV-outcome association  (so, the SoF model residual term uses first-order error)
  SEs_ZX<-c() #the s.e. of IV-exposure association; mainly for Q test for linearity
  Xmean<-c()  #stratum-specific exposure mean
  Weight_values_over_tpoints<-c()  #caution: non-negative for all strata?
  for(s in 1:max(rdat$DRstratum)){
    selected_dat<-rdat[rdat$DRstratum==s,]
    Xmean<-c( Xmean , mean( selected_dat$X )  )
    # ZX
    fitZX<-lm(    selected_dat$X ~  selected_dat$Z  );bZX<-summary(fitZX)$coef[-1,1] ;seZX<- summary(fitZX)$coef[-1,2]
    ZXest<-c( ZXest , as.numeric(bZX))
    SEs_ZX<-c( SEs_ZX , as.numeric(seZX) )
    # ZY
    fitZY<-lm(    selected_dat$Y ~  selected_dat$Z  );bZY<-summary(fitZY)$coef[-1,1] ;seZY<- summary(fitZY)$coef[-1,2]
    ZYest<-c( ZYest , as.numeric(bZY))
    SEs_ZY<-c( SEs_ZY , as.numeric(seZY) )
    # weight function
    s_weight_values_over_tpoints<-weight_function( tpoints , selected_dat$Z , selected_dat$X   )
    Weight_values_over_tpoints<-rbind(Weight_values_over_tpoints , s_weight_values_over_tpoints  )
  }

  RES$stratum_summary <- data.frame(  IV_exposure_est = ZXest ,
                                      IV_outcome_est = ZYest,
                                      IV_exposure_se = SEs_ZX,
                                      IV_outcome_se = SEs_ZY,
                                      exposure_mean = Xmean )

  Covariates<-c()
  for(pp in pos){
    S_Covariates<-c()  #single covariate value for the pp-th pos
    for(s in 1:max(rdat$DRstratum)){
      s_Covariates <- my_int_function( tpoints, Weight_values_over_tpoints[s,]*( tpoints   >= pp )    )$sum
      S_Covariates<-c( S_Covariates , s_Covariates   )
    }
    Covariates<-cbind(  Covariates , S_Covariates   )
  }
  #add the 'intercept' term (i.e. the -inf change-point)  # do we really the intercept term now if we use all exposure quantile candidates?

  dim(Covariates)





  Xmat <- Covariates                                      # X design matrix
  Sig  <-  diag(  SEs_ZY^2/ZXest^2    )                   # variance matrix of stratum-specific IV (ratio) estimator
  Beta <-  ZYest/ZXest                                    # response vector of stratum-specific IV (ratio) estimates

  #Frequentist inference
  bhat  <-    solve(   t(Xmat) %*% solve(Sig) %*%  Xmat       ) %*% t(Xmat) %*% solve(Sig) %*% Beta
  Var_bhat <- solve(   t(Xmat) %*% solve(Sig) %*%  Xmat       )

  RES$bhat <- bhat
  RES$Var_bhat <- Var_bhat

  #or direct via lm(); 使用权重进行加权回归
  #wls_model <- lm(Beta ~ Xmat - 1, weights = 1 / diag(Sig) ) # "-1" 表示不添加默认截距
  #summary(wls_model)

  time_points_for_prediction<- tpoints

  if( is.na(x_baseline_used)  ){ x_baseline_used <- mean(  X )  }
  x_baseline<-x_baseline_used

  hx_interval_estimate<-c()
  confirmed_timepoints<- pos
  for(x_current in time_points_for_prediction){
    result_vector <- outer( x_current , confirmed_timepoints, function(x, t){   (x_current -t)*((x_current -t)>0)  - (x_baseline-t)*((x_baseline-t)>0)   }       )   # F matrix/vector
    result_vector<- as.vector(result_vector  )
    # estimate at x_current and 95%  CI at x_current according to bhat and Var_bhat

    hx_interval_estimate<-rbind( hx_interval_estimate,
                                 c( sum( result_vector*bhat  )    ,
                                    sum( result_vector*bhat  )  - 1.96* sqrt(   t(result_vector) %*% Var_bhat %*% result_vector )  ,
                                    sum( result_vector*bhat  )  + 1.96* sqrt(   t(result_vector) %*% Var_bhat %*% result_vector )   )
                                )
  }

  dim( hx_interval_estimate )  # length(tpoints)  3

  gg_data2<-data.frame(  xpositions = time_points_for_prediction  ,
                         est = hx_interval_estimate[,1] ,
                         CI_low = hx_interval_estimate[,2] ,
                         CI_up = hx_interval_estimate[,3]    )
  p<- ggplot(gg_data2, aes(x = xpositions)) +
    geom_ribbon(aes(ymin = CI_low, ymax = CI_up ), fill = "grey") +
    geom_vline(xintercept = x_baseline,  color = "grey77",linetype = "dashed", linewidth = 1) +
    geom_line(aes(y = est), color = "black", linewidth = 1) +
    labs( x = "Exposure levels",  y = paste0("Estimate of h(x) with baseline x = ", x_baseline) ) +
    theme_bw() +theme(legend.position = "none")
  p

  RES$hx_ggdata<-gg_data2  #the complete information for ggplot
  RES$hx<-p  #inferred h(x) shape

  return(RES)
}


# # example
# N<-50000
# set.seed(100)
# Z <- rbinom( N , 1 , 0.5 )
# U <- rnorm( N , 0 , 1 )
# Ex <- rnorm( N , 0 , 1 )
# alphaZ<- 0.15
# X <-  alphaZ*Z + U + Ex    # weak instrument
# summary( lm(  X  ~ Z )  )$r.squared   #  < 0.01
# Ey <- rnorm( N , 0 , 1 )
# Y <- 2*(X+1)*(X>-1) -2*(X-1)*(X>1)   + U + Ey  # one change-point with effect size 1.0
#
# SSS_res<- SSS(Z,X,Y,x_baseline_used = 0)
# SSS_res$posterior_mean
#
# SS_res <- SS(Z,X,Y,x_baseline_used = 0, pos = SSS_res$posterior_mean )
# SS_res$hx




