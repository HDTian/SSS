



# Core function: complete SSS estimation



SSS<-function( Z ,
               X ,
               Y ,
               Ns_used=100, # the number of strata
               x_baseline_used = NA ,  # the baseline exposure level; default is mean(X)
               tpoints_used = 'quantile', # or 'uniform'  # the tpoints style
               boundary_used = NULL, # the user-specific vector for the left and right boundary for tpoints
               precision =100    # the precision number used in the exposure quantile region; default is 100
               ){


  RES<-list()

  if(tpoints_used == 'quantile'){
    if(is.null(boundary_used)){
      tpoints<- quantile(X, seq(0,1, length=precision+1)    ) # tpoints are used for weight function and its integration
    }else{
      if (!(is.numeric(boundary_used) &&
            length(boundary_used) == 2 &&
            all(is.finite(boundary_used)) &&
            all(boundary_used > 0 & boundary_used < 1) &&
            isTRUE(diff(boundary_used) > 0))) {
        stop("`boundary_used` must be numeric length-2 with 0 < boundary_used[1] < boundary_used[2] < 1.")
      }
      tpoints<- quantile(X, seq(boundary_used[1],boundary_used[2], length=precision+1)    )
    }

  }



  if(tpoints_used == 'uniform'){
    if(is.null(boundary_used)){
      tpoints<- seq(  min(X)  ,max(X), length=precision+1)
    }else{
      if (!(is.numeric(boundary_used) &&
            length(boundary_used) == 2 &&
            all(is.finite(boundary_used)) &&
            isTRUE(diff(boundary_used) > 0) &&
            all(boundary_used > min(X, na.rm = TRUE)) &&
            all(boundary_used < max(X, na.rm = TRUE)))) {
        stop("`boundary_used` must be numeric length-2 with min(X) <boundary_used[1] < boundary_used[2] < 1 max(X).")
      }
      tpoints<- seq(  boundary_used[1]  ,boundary_used[2], length=precision+1)
    }

  }

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
  for(s in 1:max(rdat$DRstratum)){
    s_areas <- my_int_function( tpoints, Weight_values_over_tpoints[s,]    )$areas  #length = 1000
    cumsum_reverse <- rev(cumsum(rev(s_areas))) # 累积和 从后往前计算(因为要模拟change-point之后的积分)
    #s_Covariates <- cumsum_reverse[6:95] #很关键的量！直接决定high-dim model是什么样！ 获取第某个项到第某个项的累积和
    s_Covariates <- cumsum_reverse[1:length(cumsum_reverse)]  # now use all exposure quantile candidates
    Covariates<-rbind( Covariates , s_Covariates   )
  }
  #add the 'intercept' term (i.e. the -inf change-point)  # do we really the intercept term now if we use all exposure quantile candidates?
  Covariates<-cbind(  Covariates ) # first 1 responses to the '-Inf', i.e. quantile(X,0)

  dim(Covariates)

  RES$Covariates <- Covariates

  # weighted_covariates according to s.e. of the error term in SoF regression
  responses <- ZYest/ZXest
  se_error<- SEs_ZY/ZXest
  weighted_responses <- responses/se_error
  weighted_Covariates <-Covariates / se_error

  ### SiSiE fitting ---------------------------------------------------------------------------------------------------------
  res <- susie(weighted_Covariates,weighted_responses ,L=10, intercept = FALSE,standardize = FALSE,max_iter=500,
               estimate_residual_variance = FALSE, residual_variance=1 )

  res$converged  #TRUE or FALSE indicating whether the IBSS converged to a solution within the chosen tolerance level
  res$sigma2  # plot( quantile(X, seq(0.06 , 0.95 , by=0.01) ) ,  res$pip[-1] ) #see the posterior inclusion plot

  cs_list<-susie_get_cs(res)


  ### results storage and visulization ----------------------------------------------------------------------------------------
  tpoints<-tpoints[-length(tpoints)] # as   tpoints<- quantile(X, seq(0,1, length=101)    )  包含了最右侧的点 需要去掉 使得PIP和candidate exposure points长度一致

  ### predicted change-point
  L_used<- length(   cs_list$cs )
  RES$L_used <- L_used

  PIPres<- res$alpha[1:L_used,,drop = FALSE ]  # posterior inclusion probabilities
  RES$PIPres <- PIPres

  # posterior density plot
  mode_index<- apply(  PIPres , 1, which.max ); posterior_mode<-tpoints[mode_index]
  posterior_mean <- apply(   t(PIPres)*tpoints  , 2, sum      )

  RES$posterior_mode<-posterior_mode
  RES$posterior_mean<-posterior_mean

  ggdata<-data.frame(    positions = rep(tpoints, times = L_used )   ,
                         PIP = as.vector( t(PIPres)  )  ,
                         posterior_mode = rep(posterior_mode  ,each = ncol(PIPres)  )  ,
                         posterior_mean = rep(posterior_mean  ,each = ncol(PIPres)  )  ,
                         KK  = rep(   paste0(  'Detected Change-point #' , 1:L_used )   ,each=ncol(PIPres  )     )     )

  p<-ggplot(ggdata, aes(x = positions, y = PIP)) +
    geom_point()+
    geom_line() +  # 直接连接点
    #geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs", k = 60), color='black',se = FALSE) +
    geom_vline(aes(xintercept = posterior_mode), linetype = "dashed") +
    geom_vline(aes(xintercept = posterior_mean), linetype = "dotted") +
    labs(x = "Exposure level", y = "Posterior inclusion probability") +
    facet_grid( ~ KK ) +  # 根据KK分面  # now KK: the detected change-point
    scale_x_continuous(breaks = extended_breaks(n = 9)) +  # 让刻度尽量均匀
    theme_bw() +theme(legend.position = "none")

  RES$PIP_ggdata<-ggdata  #the complete information for ggplot
  RES$PIP_plot<-p  #posterior inclusion probability plot  # dashed is the mode


  ### h(x) prediction (must matched with Covariate) ----------------------------
  time_points_for_prediction<- tpoints

  if( is.na(x_baseline_used)  ){ x_baseline_used <- mean(  X )  }
  x_baseline<-x_baseline_used


  hx_point_estimate<-c()
  for(x_current in time_points_for_prediction){
    result_matrix <- outer( x_current , time_points_for_prediction, function(x, t){   (x_current -t)*((x_current -t)>0)  - (x_baseline-t)*((x_baseline-t)>0)   }       )   # F matrix/vector
    #dim(result_matrix)  #1  91=length( time_points_for_prediction )
    posterior_mean_matrix <- res$mu[1:L_used, ,drop = FALSE ]     #dim  # L_used 91  or 10 91
    inclusion_prob_matrix <- res$alpha[1:L_used, ,drop = FALSE ]   #dim  # L_used 91  or 10 91
    poterior_mean_at_x<-result_matrix %*%  rowSums(  t( inclusion_prob_matrix * posterior_mean_matrix )  )
    hx_point_estimate<-c(hx_point_estimate, poterior_mean_at_x)
  }





  ### h(x) interval estimation -------------------------------------------------
  hx_interval_estimate<-c()
  set.seed(1123)
  N_poterior_sim<-10000  #Monte Carlo
  for(x_current in time_points_for_prediction){
    result_matrix <- outer( x_current , time_points_for_prediction, function(x, t){   (x_current -t)*((x_current -t)>0)  - (x_baseline-t)*((x_baseline-t)>0)   }       )   # F matrix/vector
    #dim(result_matrix)  #1  91=length( time_points_for_prediction )
    layer_coefficient <- result_matrix


    posterior_mean_matrix <- res$mu[1:L_used, ,drop = FALSE ]              #dim  # L_used 91  or 10 91
    posterior_2nd_moments_matrix <- res$mu2[1:L_used, ,drop = FALSE ]      #dim  # L_used 91  or 10 91
    posterior_variance_matrix <-  res$mu2[1:L_used, ,drop = FALSE ]  - ( res$mu[1:L_used, ,drop = FALSE ] )^2  #dim  # L_used 91  or 10 91

    inclusion_prob_vector <- res$alpha[1:L_used, ,drop = FALSE ]           #dim  # L_used 91  or 10 91

    increment_vectors_at_x<-c()  #length to be L_used
    for(lll in 1:L_used){
      normmix_samples_at_x <- rnormmix(N_poterior_sim, lambda = inclusion_prob_vector[lll, ]   ,
                                       mu = as.vector(   layer_coefficient * posterior_mean_matrix[lll,]            ) ,
                                       sigma = as.vector(   abs(layer_coefficient) * sqrt(posterior_variance_matrix[lll,])  )   )
      increment_vectors_at_x<-rbind(   increment_vectors_at_x ,  normmix_samples_at_x       )
    }
    posterior_sim_at_x <- apply(  increment_vectors_at_x  , 2 , sum    )  #the h(x) posterior samples at x_current

    hx_interval_estimate<-rbind( hx_interval_estimate,
                                 c(  quantile( posterior_sim_at_x, 0.025 ) ,  quantile( posterior_sim_at_x, 0.975)  )
    )
  }

  dim( hx_interval_estimate )  # length(tpoints)  2

  gg_data1<-data.frame(  xpositions = time_points_for_prediction  ,
                         posterior_sim_mean = hx_point_estimate ,
                         posterior_sim_low = hx_interval_estimate[,1] ,
                         posterior_sim_up = hx_interval_estimate[,2]    )
  p<- ggplot(gg_data1, aes(x = xpositions)) +
    geom_ribbon(aes(ymin = posterior_sim_low, ymax = posterior_sim_up), fill = "grey") +
    geom_vline(xintercept = x_baseline,  color = "grey77",linetype = "dashed", linewidth = 1) +
    geom_line(aes(y = posterior_sim_mean), color = "black", linewidth = 1) +
    labs( x = "Exposure levels",  y = paste0("Posterior estimate of h(x) with baseline x = ", x_baseline),  title = "Posterior with 95% credible interval" ) +
    theme_bw() +theme(legend.position = "none")
  p

  RES$hx_ggdata<-gg_data1  #the complete information for ggplot
  RES$hx<-p  #inferred h(x) shape

  return(RES)
}



# ### Example:
# N<-50000
# set.seed(100)
# Z <- rbinom( N , 1 , 0.5 )
# U <- rnorm( N , 0 , 1 )
# Ex <- rnorm( N , 0 , 1 )
# alphaZ<- 0.15
# X <-  alphaZ*Z + U + Ex    # weak instrument
# summary( lm(  X  ~ Z )  )$r.squared   #  < 0.01
# Ey <- rnorm( N , 0 , 1 )
# Y <- 1.0*(X-1)*(X>1)  + U + Ey  # one change-point with effect size 1.0
#
# SSS_res<- SSS(Z,X,Y,x_baseline_used = 0)
# SSS_res$PIP_plot
# SSS_res$hx
#
#
# Y <- 2*(X+1)*(X>-1) -2*(X-1)*(X>1)   + U + Ey  # one change-point with effect size 1.0
# SSS_res<- SSS(Z,X,Y,x_baseline_used = 0)
# SSS_res$PIP_plot
# SSS_res$hx




