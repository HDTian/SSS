

### this script is to reproduce the simulation results of the paper:
### Stratification-based Instrumental Variable Analysis Framework for Nonlinear Effect Analysis
### Haodong Tian, Ashish Patel, and Stephen Burgess

# Note that you do not need to install the SSS package
# Some Chinese characters may be used as comments - please bear with me

library(tidyverse)
library(data.table)
library(MendelianRandomization)
library(mixtools)  # mixtrue of normal distribution
library(DRMR)      # devtools::install_github("HDTian/DRMR")
library(susieR)
library(scales)    # for pretty_breaks


### two required functions
# integer function
my_int_function<-function(xvec, yvec1, yvec2=1  ){
  RES<-list()
  yvec<-  yvec1*yvec2
  areas<-(yvec[-length(yvec)] + yvec[-1])*(xvec[-1]  - xvec[-length(xvec)]      ) /2
  RES$sum<-sum(areas)
  RES$areas<-areas
  return(RES   )
}
#my_int_function(   seq(-10,10,length=333) ,  dnorm( seq(-10,10,length=333) , 2, 0.5  )     )$sum  #checked

# weight funciton
weight_function<-function(coords,Z , X ){    #get the smoothed and re-weighted weight function values over coords (e.g. tpoints)   #also auto re-weighted
  matrix_I <- sapply(coords, function(j) X > j)
  # get the covariance for mutiple exposure values
  mean_Z <- mean(Z)
  cov_values <- colMeans((matrix_I * 1) * Z) - colMeans(matrix_I * 1) * mean_Z
  weight_values<- as.numeric(cov_values)/cov(   Z ,X )  #length(weight_values  )  #1001

  #smoothing the weight function
  #also controlled to be positive
  weight_values[weight_values<0  ] <- 0  # all forced to be non-negative
  epsilon <- 1e-6  # small disturbance
  log_x<- log(weight_values+epsilon)  # log stransfromation
  smoothed_log_x <- stats::ksmooth(seq_along(log_x), log_x, kernel = "normal", bandwidth = 10)$y
  #seq_along(log_x) + bandwidth is more easy to control the smoothness than coords
  smoothed_weight_values<- exp(smoothed_log_x)  # back to exp; and the smoothed weight function > 0 surely

  #re-standarized
  cum_value<-my_int_function( coords , smoothed_weight_values   )$sum
  reweight_values<-smoothed_weight_values/cum_value

  #my_int_function( coords , reweight_values   )$sum #1 #checked; integral of the weight function is 1
  return(  reweight_values  )
}




### Fig 1 ######################################################################
### Changepoint prediction boxpolot for 4 Z-X model scenarios (2 IVtype * 2 Xtype) and one change-point model

N<-50000  #as real size is 143 000
Ns_used<-100 #for better change-point capture ability
simtimes<-1000


Posterior_predict<-array(NA, dim = c(4, 2, simtimes))  # 4 scenarios * posterior mean and mode * simtimes



for(iii in 1:simtimes){

  posterior_predict<-c()
  for(scenarios_used in 1:4){
    set.seed(iii)


    ### individual-level data
    U <- rnorm( N , 0 , 1 )
    Ex <- rnorm( N , 0 , 1 )
    Ey <- rnorm( N , 0 , 1 )
    if(scenarios_used == 1){   Z <- rbinom( N, 1, 0.5)-0.5;  alphaZ<- 0.15 ;  X <-        alphaZ*Z + U + Ex          ; Y <- 1.0*(X-0)*(X>0)  + U + Ey  }
    if(scenarios_used == 2){   Z <- rnorm( N , 0, 1 )  ;     alphaZ<- 0.15 ;  X <-        alphaZ*Z + U + Ex          ; Y <- 1.0*(X-0)*(X>0)  + U + Ey  }
    if(scenarios_used == 3){   Z <- rbinom( N , 1, 0.5 )  ;  alphaZ<- 0.3  ;  X <-  exp(  alphaZ*Z + U + Ex        ) ; Y <- 1.0*(X-2.5)*(X>2.5)  + U + Ey   }
    if(scenarios_used == 4){   Z <- rnorm( N , 0, 1 )  ;     alphaZ<- 0.3  ;  X <-  exp(  alphaZ*Z + U + Ex        ) ; Y <- 1.0*(X-2.5)*(X>2.5)  + U + Ey  }




    #if(yyy ==1){Y <- 1.0*(X-0.5)*(X>0.5)  + U + Ey}  # one change-point
    #if(yyy ==2){Y <- 1.0*(X-0.5)*(X>0.5)  - 0.5*(X-2)*(X>2) + U + Ey}  # two change-point
    #if(yyy ==3){Y <- 1.0*(X-0.5)*(X>0.5)  - 0.5*(X-1)*(X>1)  - 0.5*(X-2)*(X>2)  + U + Ey}  # three change-point

    ###
    tpoints<- quantile(X, seq(0,1, length=101)    ) # tpoints are used for weight function and its integration

    ### DR stratification with strata number = Ns_used
    dat_<-data.frame(  Z=Z, X=X, Y=Y )   # X:= alcohol   with outcome SBP
    rdat<-Stratify(dat_, Ns = Ns_used) #DR stratification

    ### stratum-specific results
    ZXest<-c()  #stratum-specific IV-exposure association
    ZYest<-c()  #stratum-specific IV-outcome association
    SEs_ZY<-c() #the s.e. of IV-outcome association  (so, the SoF model residual term uses first-order error)
    SEs_ZX<-c() #the s.e. of IV-exposure association; mainly for Q test for linearity
    Xmean<-c()  #stratum-specific eposure mean
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

    ### Covariate matrix; using 0 + 95 or 91 quantile levels as the candidates
    Covariates<-c()
    for(s in 1:max(rdat$DRstratum)){
      s_areas <- my_int_function( tpoints, Weight_values_over_tpoints[s,]    )$areas  #length = 1000
      cumsum_reverse <- rev(cumsum(rev(s_areas))) # 累积和 从后往前计算(因为要模拟change-point之后的积分)
      if(scenarios_used %in% c(3,4)){  s_Covariates <- cumsum_reverse[1:95]      }
      if(scenarios_used %in% c(1,2)){  s_Covariates <- cumsum_reverse[6:95]      } #很关键的量！直接决定high-dim model是什么样！ 获取第某个项到第某个项的累积和
      Covariates<-rbind( Covariates , s_Covariates   )
    }
    #add the 'intercept' term (i.e. the -inf chenge-point)
    Covariates<-cbind( 1, Covariates ) # looks good - all positive and all decreasing
    dim(Covariates)  # Ns   96 = 1 + 95 quantiles over [0.0th-quantial , 0.95th-quantile ]

    # weighted_covariates according to s.e. of the error term in SoF regression
    responses <- ZYest/ZXest
    se_error<- SEs_ZY/ZXest
    weighted_responses <- responses/se_error
    weighted_Covariates <-Covariates / se_error

    ### SiSiE fitting
    res <- susie(weighted_Covariates,weighted_responses ,L=10, intercept = FALSE,standardize = FALSE,
                 estimate_residual_variance = FALSE, residual_variance=1 )

    res$converged  #TRUE or FALSE indicating whether the IBSS converged to a solution within the chosen tolerance level
    res$sigma2

    cs_list<-susie_get_cs(res)
    cs_list$cs


    ### record the posterior mode and mean for each possible detected parameter
    if(scenarios_used %in% c(3,4)){      possible_points<-c(   0 , quantile( X, seq(0.01 , 0.95 , by=0.01) )  )       }     #1 + 95  candidate exposure levels
    if(scenarios_used %in% c(1,2)){      possible_points<-c(   quantile(X,0) , quantile( X, seq(0.06 , 0.95 , by=0.01) )  )       }  #1 + 90  candidate exposure levels


    predict_mode<-possible_points[ which.max( res$alpha[1,] )  ]   # posterior mode prediction for the first (detected) parameter
    predict_mean<-sum(   possible_points * (  res$alpha[1,] )  )   # posterior mean prediction for the first (detected) parameter

    # posterior_predict<-matrix(  NA,4,2   )
    # for(ppp in 1:(  min( length(cs_list$cs), 3 )  )  ){
    #   posterior_predict[ppp,1] <- possible_points[ which.max( res$alpha[ppp,] )  ]   # posterior mode prediction
    #   posterior_predict[ppp,2] <- sum(   possible_points * (  res$alpha[ppp,] )  )   # posterior mean prediction
    # } #remaining enteries left to be 0
    #
    # if(yyy ==1){   Posterior_predict[,,iii]  <- posterior_predict   }
    # if(yyy ==2){   Posterior_predict2[,,iii] <- posterior_predict   }
    # #if(yyy ==3){   Posterior_predict3[,,iii] <- posterior_predict   }


    posterior_predict<-rbind( posterior_predict , c(predict_mode , predict_mean)  )
  }  # end of 4 scenarios loop

  Posterior_predict[,,iii] <- posterior_predict
}

dim( Posterior_predict )  #4   2 100
#View(  Posterior_predict[1,,]   )

saveRDS(Posterior_predict, file = "/Users/haodongtian/Documents/Project SoF/code/Posterior_predict.rds")
#Posterior_predict <- readRDS("/Users/haodongtian/Documents/Project SoF/code/Posterior_predict.rds")
#Posterior_predict <- readRDS('/Users/htian/Library/Mobile Documents/com~apple~CloudDocs/Documents/Project SoF/code/Posterior_predict.rds')




predicte_values  <- as.vector( Posterior_predict )  #order is: array[1->4, 1 , 1]  array[1->4, 2 ,1  ] array[1->4, 1 , 2] array[1->4, 2 ,2  ] ...
df<- data.frame(  predicte_values =predicte_values,
                  Scenario = c("Scenario 1", "Scenario 2", "Scenario 3", "Scenario 4")  ,
                  ResultType =  rep(  c(  'Mode' , 'Mean'   ) , each = 4    )
)

yrange1<- qmnorm(  c( 0.25 ,0.5 , 0.75 )        )                        ; yrange1 <-round(yrange1,1)
yrange2<- qnorm( c( 0.25 ,0.5 , 0.75 )   , 0 , sqrt(0.15^2 + 1 +1)  )    ; yrange2 <-round(yrange2,1)
yrange3<- qmlnorm(  c( 0.25 ,0.5 , 0.75 )        )                       ; yrange3 <-round(yrange3,1)
yrange4<- qlnorm( c( 0.25 ,0.5 , 0.75 )   , 0 , sqrt( 0.3^2 + 1 +1)  )   ; yrange4 <-round(yrange4,1)


ggdata<- df[df$Scenario == "Scenario 1" ,  ]
p1<-ggplot(ggdata, aes(x = ResultType, y = predicte_values )) +
  geom_hline(yintercept = 0, color = "lightblue", linewidth = 1) +
  geom_boxplot(fill='grey') +
  labs(x='', y = "Predicted change-point",title='Scenario 1') +
  coord_cartesian(ylim = c(min(yrange1),  max(yrange1)))  +
  scale_y_continuous(breaks = yrange1) +  # 使用 yrange1 定义 y 轴刻度
  theme_bw()+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
p1

ggdata<- df[df$Scenario == "Scenario 2" ,  ]
p2<-ggplot(ggdata, aes(x = ResultType, y = predicte_values )) +
  geom_hline(yintercept = 0, color = "lightblue", linewidth = 1) +
  geom_boxplot(fill='grey') +
  labs(x='', y = " ",title='Scenario 2') +
  coord_cartesian(ylim = c(min(yrange2),  max(yrange2)))  +
  scale_y_continuous(breaks = yrange2) +  # 使用 yrange1 定义 y 轴刻度
  theme_bw()+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
p2


ggdata<- df[df$Scenario == "Scenario 3" ,  ]
p3<-ggplot(ggdata, aes(x = ResultType, y = predicte_values )) +
  geom_hline(yintercept = 2.5, color = "lightblue", linewidth = 1) +
  geom_boxplot(fill='grey') +
  labs(x='', y = " ",title='Scenario 3') +
  coord_cartesian(ylim = c(min(yrange3),  max(yrange3)))  +
  scale_y_continuous(breaks = yrange3) +  # 使用 yrange1 定义 y 轴刻度
  theme_bw()+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
p3

ggdata<- df[df$Scenario == "Scenario 4" ,  ]
p4<-ggplot(ggdata, aes(x = ResultType, y = predicte_values )) +
  geom_hline(yintercept = 2.5, color = "lightblue", linewidth = 1) +
  geom_boxplot(fill='grey') +
  labs(x='', y = " ",title='Scenario 4') +
  coord_cartesian(ylim = c(min(yrange4),  max(yrange4)))  +
  scale_y_continuous(breaks = yrange4) +  # 使用 yrange1 定义 y 轴刻度
  theme_bw()+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
p4
library(cowplot)
plot_grid(p1, p2, p3, p4, ncol = 4)
#sim_Fig1



### Fig 2 ######################################################################
### boxplot of the point-wise predicted h(x)
###
#   IV-X case:  Case 1  Case 2  Case 3  Case 4
# h scenario 1
# h scenario 2
# h scenario 3
# h scenario 4


### nonnegative exposure scenairos (lognormal and mix of lognormal) ------------

#  Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
#0.0019   0.3783   1.0048   2.8710   2.6572 375.4682
#0.0016   0.4500   1.1620   3.2028   3.0577 782.5293

N<-50000  #as real size is 143 000
simtimes<-1000

Ns_used<-100 #for better change-point capture ability；Ns=10 is also okay, only move to appendix

hx_point_Estimate_total<-c()
for(hx_used in 1:4){
  if(hx_used==1){ hx<-function(X){1.0*(X-0)*(X>0)}  }                                        #[h scenario 1] linear scenario
  #hx<-function(X){1.0*(X-0.5)*(X>0.5)}                                   #nonlinear with one change-point
  if(hx_used==2){ hx<-function(X){1.0*(X-2.5)*(X>2.5)}  }                                    #[h scenario 2]
  if(hx_used==3){ hx<-function(X){0.5*X + 0.5*(X-0.5)*(X>0.5)  + 0.5*(X-2.5)*(X>2.5) }  }    #[h scenario 3] nonlinear with two change-point (3 parameters); good in both continuous and binary IV
  if(hx_used==4){ hx<-function(X){ -2*X + 0.5*X^2   }  }                                     #[h scenario 4] other nonlinear case


  for(Ztype in c('continuous' , 'binary' )){

    hx_point_Estimate<-c()  # simtimes  96

    for(iii in 1:simtimes){

      set.seed(iii)
      ### individual-level data
      if(Ztype == 'continuous'){Z <- rnorm( N , 0 , 1 ) } #continuous IV
      if(Ztype == 'binary'){    Z <- rbinom(N, 1, 0.5 ) } #binary IV
      U <- rnorm( N , 0 , 1 )
      Ex <- rnorm( N , 0 , 1 )
      ### lognormal expoure
      alphaZ<- 0.3  #same for any IV type
      X <-  exp(  alphaZ*Z + U + Ex        )
      Ey <- rnorm( N , 0 , 1 )
      Y <- hx(X) + U + Ey


      tpoints<- quantile(X, seq(0,1, length=101)    ) # tpoints are used for weight function and its integration

      ### DR stratification with strata number = Ns_used
      dat_<-data.frame(  Z=Z, X=X, Y=Y )   # X:= alcohol   with outcome SBP
      rdat<-Stratify(dat_, Ns = Ns_used) #DR stratification


      ### stratum-specific results
      ZXest<-c()  #stratum-specific IV-exposure association
      ZYest<-c()  #stratum-specific IV-outcome association
      SEs_ZY<-c() #the s.e. of IV-outcome association  (so, the SoF model residual term uses first-order error)
      SEs_ZX<-c() #the s.e. of IV-exposure association; mainly for Q test for linearity
      Xmean<-c()  #stratum-specific eposure mean
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

      ### Covariate matrix; using 0 + 95 quantile levels as the candidates
      Covariates<-c()
      for(s in 1:max(rdat$DRstratum)){
        s_areas <- my_int_function( tpoints, Weight_values_over_tpoints[s,]    )$areas  #length = 1000
        cumsum_reverse <- rev(cumsum(rev(s_areas))) # 累积和 从后往前计算(因为要模拟change-point之后的积分)
        s_Covariates <- cumsum_reverse[1:95] #很关键的量！直接决定high-dim model是什么样！ 获取第某个项到第某个项的累积和
        # length(s_Covariates) # 801   # and also decreasing
        Covariates<-rbind( Covariates , s_Covariates   )
      }
      #add the 'intercept' term (i.e. the -inf chenge-point)
      Covariates<-cbind( 1, Covariates ) # looks good - all positive and all decreasing
      dim(Covariates)  # Ns   96 = 1 + 95 quantiles over [0.0th-quantial , 0.95th-quantile ]

      # weighted_covariates according to s.e. of the error term in SoF regression
      responses <- ZYest/ZXest
      se_error<- SEs_ZY/ZXest
      weighted_responses <- responses/se_error
      weighted_Covariates <-Covariates / se_error

      ### SuSiE fitting
      res <- susie(weighted_Covariates,weighted_responses ,L=10, intercept = FALSE,standardize = FALSE,max_iter=500,
                   estimate_residual_variance = FALSE, residual_variance=1 )

      res$converged  #TRUE or FALSE indicating whether the IBSS converged to a solution within the chosen tolerance level
      res$sigma2

      cs_list<-susie_get_cs(res)

      ### h(x) prediction
      time_points_for_prediction<- c(   0 , quantile(X, seq(0.01 , 0.95 , by=0.01) )  ) #1 + 95 quantiles over [0.0th-quantial , 0.95th-quantile ]

      x_baseline<-0
      L_used<-length( cs_list$cs   )
      hx_point_estimate<-c()  #length = 5
      for(x_current in time_points_for_prediction){
        result_matrix <- outer( x_current , time_points_for_prediction, function(x, t){   (x_current -t)*((x_current -t)>0)  - (x_baseline-t)*((x_baseline-t)>0)   }       )   # F matrix/vector
        #dim(result_matrix)  #1  91=length( time_points_for_prediction )
        posterior_mean_matrix <- res$mu[1:L_used, ,drop = FALSE ]     #dim  # L_used 91  or 10 91
        inclusion_prob_matrix <- res$alpha[1:L_used, ,drop = FALSE ]   #dim  # L_used 91  or 10 91
        poterior_mean_at_x<-result_matrix %*%  rowSums(  t( inclusion_prob_matrix * posterior_mean_matrix )  )
        hx_point_estimate<-c(hx_point_estimate, poterior_mean_at_x)
      }

      hx_point_Estimate<-rbind(hx_point_Estimate , hx_point_estimate   )

    }  #end of hx_point_Estimate
    print(  paste0(  'finihed: ' , hx_used , '  ' , Ztype  ) )
    dim(hx_point_Estimate)# simtimes 96= 1 + 95 quantiles over [0.0th-quantial , 0.95th-quantile ]

    hx_point_Estimate <- cbind( hx_point_Estimate , hx_used ,  Ztype  )

    hx_point_Estimate_total <- rbind( hx_point_Estimate_total , hx_point_Estimate )
  }


}

dim( hx_point_Estimate_total )  #  8*simtimes    98 = 1 + 95 + hx_used +  Ztype

saveRDS(hx_point_Estimate_total, file = "/Users/haodongtian/Documents/Project SoF/code/hx_point_Estimate_total_lognorm.rds")
#hx_point_Estimate_total <- readRDS("/Users/haodongtian/Documents/Project SoF/code/hx_point_Estimate_total_lognorm.rds")
#hx_point_Estimate_total<- readRDS( '/Users/htian/Library/Mobile Documents/com~apple~CloudDocs/Documents/Project SoF/code/hx_point_Estimate_total_lognorm.rds'  )

# 混合对数正态的累积分布函数
F_X_log <- function(q) {
  0.5 * plnorm(q, meanlog = 0, sdlog = sqrt(2)) +
    0.5 * plnorm(q, meanlog = 0.3, sdlog = sqrt(2))
}
# 求分位点函数
quantile_mixed_log <- function(p) {
  uniroot(function(q) F_X_log(q) - p, interval = c(0, 1e4))$root
}
qmlnorm<-function(vec){  sapply(  vec    , quantile_mixed_log   )       }

GGdata<-c()
for(hx_used in 1:4){
  if(hx_used==1){ hx<-function(X){1.0*(X-0)*(X>0)}  }                                        #[h scenario 1] linear scenario
  #hx<-function(X){1.0*(X-0.5)*(X>0.5)}     #nonlinear with one change-point
  if(hx_used==2){ hx<-function(X){1.0*(X-2.5)*(X>2.5)}  }                                    #[h scenario 2]
  if(hx_used==3){ hx<-function(X){0.5*X + 0.5*(X-0.5)*(X>0.5)  + 0.5*(X-2.5)*(X>2.5) }  }    #[h scenario 3] nonlinear with two change-point (3 parameters); good in both continuous and binary IV
  if(hx_used==4){ hx<-function(X){ -2*X + 0.5*X^2   }  }                                     #[h scenario 4] other nonlinear case

  for(Ztype in c('continuous' , 'binary' )){
    current_index <- (hx_point_Estimate_total[ ,97 ] == hx_used)&( hx_point_Estimate_total[ ,98 ] == Ztype    )
    hx_point_Estimate <- hx_point_Estimate_total[ current_index  ,  1:96 ]
    hx_point_Estimate <- apply(hx_point_Estimate, 2, as.numeric)

    if(Ztype == 'continuous'){    theoretical_quantiles<- qlnorm( seq(0, 0.95, by=0.01) , 0 , sqrt(alphaZ^2 + 1 +1)  )     }
    if(Ztype == 'binary'){    theoretical_quantiles<- qmlnorm( seq(0, 0.95, by=0.01) )                                 }
    ggdata<-data.frame(  x = theoretical_quantiles  ,
                         mean = apply(  hx_point_Estimate  , 2, mean    )  ,
                         lower =  apply(  hx_point_Estimate  , 2, function(vec){  quantile(vec , 0.025) }    ) ,
                         upper =  apply(  hx_point_Estimate  , 2, function(vec){  quantile(vec , 0.975) }    ) ,
                         hx_used = hx_used,
                         Ztype  = Ztype
    )
    ggdata$true <- hx(  ggdata$x   )
    GGdata<- rbind(  GGdata, ggdata)
  }
}
dim(GGdata)

GGdata$Ztype[GGdata$Ztype=='continuous' ]<-'Continuous IV'
GGdata$Ztype[GGdata$Ztype=='binary' ]<-'Binary IV'

GGdata$hx_used <- sapply( GGdata$hx_used  , function(x){   paste0( 'Case ' , x    ) }    )

#only show with certain lower exposure levels
GGdata<- GGdata[ GGdata$x <=6  , ]

p<-ggplot(GGdata, aes(x = x)) +
  geom_hline(yintercept = 0 , color = "grey", linetype = "dashed", linewidth = 0.7) +
  geom_vline(xintercept = 0.5 , color = "grey", linetype = "dashed", linewidth = 0.7) +
  geom_vline(xintercept = 2.5 , color = "grey", linetype = "dashed", linewidth = 0.7) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "grey60") +  # 95% 区域
  geom_line(aes(y = mean), color = "black", linewidth = 1) +  # 模拟均值曲线
  geom_line(aes(y = true), color = "lightblue", linewidth = 1) +  # 模拟均值曲线
  facet_grid( hx_used ~ Ztype, scales = "free") +  # 根据KK分面
  labs(
    x = "Exposure level",
    y = "Predicted effect",
    title ='Lognormal exposure'
  ) +
  scale_x_continuous(breaks = pretty_breaks(n = 5)) + # 自动设置 20 个主要刻度
  theme_bw()+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
p
#sim_Fig2_right  400 450

GGdata_34<- GGdata #for the whole Part III Fig

### summary:
# when # changepoint = 1, the predicted h(x) works very well
# also when Ns is larger, the prediction performance is improved

# when  # changepoint = 1, the predicted h(x) works litter worse, but still accepatable






### continuous exposure scenairos (normal and mix of normal) ------------
###
###
###
###

#    Min.   1st Qu.    Median      Mean   3rd Qu.      Max.
#-7.039574 -0.941549  0.002219  0.001284  0.942582  5.691923
#-6.45705 -0.86894  0.07750  0.07917  1.03218  6.66253

N<-50000  #as real size is 143 000
simtimes<-1000

Ns_used<-100 #for better change-point capture ability；Ns=10 is also okay, only move to appendix

hx_point_Estimate_total<-c()
for(hx_used in 1:4){     #hx(0) = 0 must hold
  if(hx_used==1){ hx<-function(X){1.0*X}  }                                          #[h scenario 1] linear scenario
  if(hx_used==2){ hx<-function(X){1.0*(X-0)*(X>0)}  }                                #[h scenario 2]
  if(hx_used==3){ hx<-function(X){0.5*X + 0.5*(X+0.5)*(X>-0.5)-0.25  + 0.5*(X-0.5)*(X>0.5) }  }    #[h scenario 3] nonlinear with two change-point (3 parameters); good in both continuous and binary IV
  if(hx_used==4){ hx<-function(X){ -1*X + 0.5*X^2   }  }                             #[h scenario 4] other nonlinear case


  for(Ztype in c('continuous' , 'binary' )){

    hx_point_Estimate<-c()  # simtimes  96

    for(iii in 1:simtimes){

      set.seed(iii)
      ### individual-level data
      if(Ztype == 'continuous'){Z <- rnorm( N , 0 , 1 )      } #continuous IV
      if(Ztype == 'binary'){    Z <- rbinom(N, 1, 0.5 ) -0.5 } #binary IV  (for continuous X, centerize Z so that basline exposure level 0 is interpretable)
      U <- rnorm( N , 0 , 1 )
      Ex <- rnorm( N , 0 , 1 )
      ### normal expoure
      alphaZ<- 0.15  #same for any IV type
      X <-  alphaZ*Z + U + Ex
      Ey <- rnorm( N , 0 , 1 )
      Y <- hx(X) + U + Ey


      tpoints<- quantile(X, seq(0,1, length=101)    ) # tpoints are used for weight function and its integration

      ### DR stratification with strata number = Ns_used
      dat_<-data.frame(  Z=Z, X=X, Y=Y )   # X:= alcohol   with outcome SBP
      rdat<-Stratify(dat_, Ns = Ns_used) #DR stratification


      ### stratum-specific results
      ZXest<-c()  #stratum-specific IV-exposure association
      ZYest<-c()  #stratum-specific IV-outcome association
      SEs_ZY<-c() #the s.e. of IV-outcome association  (so, the SoF model residual term uses first-order error)
      SEs_ZX<-c() #the s.e. of IV-exposure association; mainly for Q test for linearity
      Xmean<-c()  #stratum-specific eposure mean
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

      ### Covariate matrix; using 0 + 95 quantile levels as the candidates
      Covariates<-c()
      for(s in 1:max(rdat$DRstratum)){
        s_areas <- my_int_function( tpoints, Weight_values_over_tpoints[s,]    )$areas  #length = 1000
        cumsum_reverse <- rev(cumsum(rev(s_areas))) # 累积和 从后往前计算(因为要模拟change-point之后的积分)
        s_Covariates <- cumsum_reverse[6:95] #很关键的量！直接决定high-dim model是什么样！ 获取第某个项到第某个项的累积和
        # length(s_Covariates) # 801   # and also decreasing
        Covariates<-rbind( Covariates , s_Covariates   )
      }
      #add the 'intercept' term (i.e. the -inf chenge-point)
      Covariates<-cbind( 1, Covariates ) # first 1 responses to the '-Inf', i.e. quantile(X,0)
      dim(Covariates)  # Ns   91 = 1 + 90 quantiles over [0.05th-quantial , 0.95th-quantile ]

      # weighted_covariates according to s.e. of the error term in SoF regression
      responses <- ZYest/ZXest
      se_error<- SEs_ZY/ZXest
      weighted_responses <- responses/se_error
      weighted_Covariates <-Covariates / se_error

      ### SiSiE fitting
      res <- susie(weighted_Covariates,weighted_responses ,L=10, intercept = FALSE,standardize = FALSE,max_iter=500,
                   estimate_residual_variance = FALSE, residual_variance=1 )

      res$converged  #TRUE or FALSE indicating whether the IBSS converged to a solution within the chosen tolerance level
      res$sigma2

      cs_list<-susie_get_cs(res)

      ### h(x) prediction (must mached with Covariate)
      time_points_for_prediction<- c(  quantile(X,0)  ,  quantile(X, seq(0.06 , 0.95 , by=0.01) )  ) #1 + 91 quantiles over [0.05th-quantial , 0.95th-quantile ]

      x_baseline<-0
      L_used<-length( cs_list$cs   )
      hx_point_estimate<-c()  #length = 5
      for(x_current in time_points_for_prediction){
        result_matrix <- outer( x_current , time_points_for_prediction, function(x, t){   (x_current -t)*((x_current -t)>0)  - (x_baseline-t)*((x_baseline-t)>0)   }       )   # F matrix/vector
        #dim(result_matrix)  #1  91=length( time_points_for_prediction )
        posterior_mean_matrix <- res$mu[1:L_used, ,drop = FALSE ]     #dim  # L_used 91  or 10 91
        inclusion_prob_matrix <- res$alpha[1:L_used, ,drop = FALSE ]   #dim  # L_used 91  or 10 91
        poterior_mean_at_x<-result_matrix %*%  rowSums(  t( inclusion_prob_matrix * posterior_mean_matrix )  )
        hx_point_estimate<-c(hx_point_estimate, poterior_mean_at_x)
      }

      hx_point_Estimate<-rbind(hx_point_Estimate , hx_point_estimate   )

    }  #end of hx_point_Estimate
    print(  paste0(  'finihed: ' , hx_used , '  ' , Ztype  ) )
    dim(hx_point_Estimate)# simtimes 91= 1 + 90 quantiles over [0.05th-quantial , 0.95th-quantile ]


    hx_point_Estimate <- cbind( hx_point_Estimate , hx_used ,  Ztype  )
    hx_point_Estimate_total <- rbind( hx_point_Estimate_total , hx_point_Estimate )
  }


}

dim( hx_point_Estimate_total )  #  8*simtimes    93 = 1 + 90 + hx_used +  Ztype

saveRDS(hx_point_Estimate_total, file = "/Users/haodongtian/Documents/Project SoF/code/hx_point_Estimate_total_norm.rds")
#hx_point_Estimate_total <- readRDS("/Users/haodongtian/Documents/Project SoF/code/hx_point_Estimate_total_norm.rds")
#hx_point_Estimate_total<-readRDS( '/Users/htian/Library/Mobile Documents/com~apple~CloudDocs/Documents/Project SoF/code/hx_point_Estimate_total_norm.rds' )

# 混合正态的累积分布函数
F_X <- function(q) {
  0.5 * pnorm(q, mean = -0.5, sd = sqrt(2)) +
    0.5 * pnorm(q, mean = 0.5, sd = sqrt(2))
}
# 求分位点函数
quantile_mixed <- function(p) {
  if(p==0){ -Inf  }else{
    uniroot(function(q) F_X(q) - p, interval = c(-1e4, 1e4))$root
  }
}
qmnorm<-function(vec){  sapply(  vec    , quantile_mixed   )       }

GGdata<-c()
for(hx_used in 1:4){
  if(hx_used==1){ hx<-function(X){1.0*X}  }                                          #[h scenario 1] linear scenario
  if(hx_used==2){ hx<-function(X){1.0*(X-0)*(X>0)}  }                                #[h scenario 2]
  if(hx_used==3){ hx<-function(X){0.5*X + 0.5*(X+0.5)*(X>-0.5)-0.25  + 0.5*(X-0.5)*(X>0.5) }  }  #[h scenario 3] nonlinear with two change-point (3 parameters); good in both continuous and binary IV
  if(hx_used==4){ hx<-function(X){ -1*X + 0.5*X^2    }  }                             #[h scenario 4] other nonlinear case

  for(Ztype in c('continuous' , 'binary' )){
    current_index <- (hx_point_Estimate_total[ ,(ncol(hx_point_Estimate_total)-1) ] == hx_used)&( hx_point_Estimate_total[ ,ncol(hx_point_Estimate_total) ] == Ztype    )
    hx_point_Estimate <- hx_point_Estimate_total[ current_index  ,  1:91 ]
    hx_point_Estimate <- apply(hx_point_Estimate, 2, as.numeric)

    if(Ztype == 'continuous'){    theoretical_quantiles<- qnorm( c(0,  seq(0.06, 0.95, by=0.01) ) , 0 , sqrt(alphaZ^2 + 1 +1)  )     }
    if(Ztype == 'binary'){        theoretical_quantiles<- qmnorm( c(0,  seq(0.06, 0.95, by=0.01) )  )                   }
    ggdata<-data.frame(  x = theoretical_quantiles  ,
                         mean = apply(  hx_point_Estimate  , 2, mean    )  ,
                         lower =  apply(  hx_point_Estimate  , 2, function(vec){  quantile(vec , 0.025) }    ) ,
                         upper =  apply(  hx_point_Estimate  , 2, function(vec){  quantile(vec , 0.975) }    ) ,
                         hx_used = hx_used,
                         Ztype  = Ztype
    )
    ggdata$true <- hx(  ggdata$x   )
    GGdata<- rbind(  GGdata, ggdata)
  }
}
dim(GGdata)

GGdata$Ztype[GGdata$Ztype=='continuous' ]<-'Continuous IV'
GGdata$Ztype[GGdata$Ztype=='binary' ]<-'Binary IV'

GGdata$hx_used <- sapply( GGdata$hx_used  , function(x){   paste0( 'Case ' , x    ) }    )

#only show with appripriate exposure region
GGdata<- GGdata[ (GGdata$x <= 1.5)&(GGdata$x >= -1.5  )  , ]

p<-ggplot(GGdata, aes(x = x)) +
  geom_hline(yintercept = 0 , color = "grey", linetype = "dashed", linewidth = 0.7) +
  geom_vline(xintercept = 0 , color = "grey", linetype = "dashed", linewidth = 0.7) +
  geom_vline(xintercept = -0.5 , color = "grey", linetype = "dashed", linewidth = 0.7) +
  geom_vline(xintercept = 0.5 , color = "grey", linetype = "dashed", linewidth = 0.7) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "grey60") +  # 95% 区域
  geom_line(aes(y = mean), color = "black", linewidth = 1) +  # 模拟均值曲线
  geom_line(aes(y = true), color = "lightblue", linewidth = 1) +  # 模拟均值曲线
  facet_grid( hx_used ~ Ztype, scales = "free") +  # 根据KK分面
  labs(
    x = "Exposure level",
    y = "Predicted effect",
    title = 'Normal exposure'
  ) +
  scale_x_continuous(breaks = pretty_breaks(n = 10)) + # 自动设置 20 个主要刻度
  theme_bw()+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
p
#sim_Fig2_left  400 450


GGdata_12<- GGdata #for the whole Part III Fig


GGdata_12$Scenario <-  paste0(  'Scenario '  ,   as.character(  ( GGdata_12$Ztype == 'Continuous IV' ) + 1  )   )
GGdata_34$Scenario <-  paste0(  'Scenario '  , as.character(  ( GGdata_34$Ztype == 'Continuous IV' ) + 1 +2  )  )

GGdata_1234 <-  rbind(GGdata_12  ,  GGdata_34  )

p_h<-ggplot(GGdata_1234, aes(x = x)) +
  geom_hline(yintercept = 0 , color = "grey", linetype = "dashed", linewidth = 0.7) +
  geom_vline(data = df %>% filter(Scenario %in% c("Scenario 1", "Scenario 2")),
             aes(xintercept = 0), color = "grey", linetype = "dashed", linewidth = 0.7)+
  geom_vline(data = df %>% filter(Scenario %in% c("Scenario 1", "Scenario 2")),
             aes(xintercept = -0.5), color = "grey", linetype = "dashed", linewidth = 0.7)+
  geom_vline(data = df %>% filter(Scenario %in% c("Scenario 1", "Scenario 2")),
             aes(xintercept = 0.5), color = "grey", linetype = "dashed", linewidth = 0.7)+
  geom_vline(data = df %>% filter(Scenario %in% c("Scenario 3", "Scenario 4")),
             aes(xintercept = 0), color = "grey", linetype = "dashed", linewidth = 0.7)+
  geom_vline(data = df %>% filter(Scenario %in% c("Scenario 3", "Scenario 4")),
             aes(xintercept = 0.5), color = "grey", linetype = "dashed", linewidth = 0.7)+
  geom_vline(data = df %>% filter(Scenario %in% c("Scenario 3", "Scenario 4")),
             aes(xintercept = 2.5), color = "grey", linetype = "dashed", linewidth = 0.7)+
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "grey60") +  # 95% 区域
  geom_line(aes(y = mean), color = "black", linewidth = 1) +  # 模拟均值曲线
  geom_line(aes(y = true), color = "lightblue", linewidth = 1) +  # 模拟均值曲线
  facet_grid( hx_used ~ Scenario, scales = "free") +  # 根据KK分面
  labs(
    x = "Exposure level",
    y = "Predicted effect"
  ) +
  scale_x_continuous(breaks = extended_breaks(n = 7)) +  # 让刻度尽量均匀
  #scale_x_continuous(breaks = pretty_breaks(n = 10)) + # 自动设置 20 个主要刻度
  theme_bw()+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
p_h
#sim_Fig2_all  700 550



### combined both plots in the sim
#library(cowplot)
#plot_grid(p1, p2, p3, p4, ncol = 4)


library(grid)

grid.newpage()
pushViewport(viewport(layout = grid.layout(100 , 100 )))
print(p1, vp = viewport(layout.pos.row = 5:30, layout.pos.col = 10:30))
print(p2, vp = viewport(layout.pos.row = 5:30, layout.pos.col = 31:50))
print(p3, vp = viewport(layout.pos.row = 5:30, layout.pos.col = 51:70))
print(p4, vp = viewport(layout.pos.row = 5:30, layout.pos.col = 71:90))
print(p_h, vp = viewport(layout.pos.row = 30:100, layout.pos.col = 10:90))
#780 720



