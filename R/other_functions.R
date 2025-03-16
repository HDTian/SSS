

# other functions


my_int_function<-function(xvec, yvec1, yvec2=1  ){
  RES<-list()
  yvec<-  yvec1*yvec2
  areas<-(yvec[-length(yvec)] + yvec[-1])*(xvec[-1]  - xvec[-length(xvec)]      ) /2
  RES$sum<-sum(areas)
  RES$areas<-areas
  return(RES   )
}
#my_int_function(   seq(-10,10,length=333) ,  dnorm( seq(-10,10,length=333) , 2, 0.5  )     )$sum  #checked


weight_function<-function(coords,Z , X ){    #get the smoothed and re-weighted weight function values over coords (e.g. tpoints)   #also auto re-weighted 
  matrix_I <- sapply(coords, function(j) X > j)
  # 计算协方差
  mean_Z <- mean(Z)
  cov_values <- colMeans((matrix_I * 1) * Z) - colMeans(matrix_I * 1) * mean_Z
  weight_values<- as.numeric(cov_values)/cov(   Z ,X )  #length(weight_values  )  #1001
  
  #smoothed #also controlled to be positive (!)
  weight_values[weight_values<0  ] <- 0  # all forced to be non-negative
  epsilon <- 1e-6  # 小偏移量
  log_x<- log(weight_values+epsilon)  # 对数变换，x 必须严格大于 0
  smoothed_log_x <- stats::ksmooth(seq_along(log_x), log_x, kernel = "normal", bandwidth = 10)$y  #seq_along(log_x)配合bandwidth相比coords更好控制平滑度一些
  smoothed_weight_values<- exp(smoothed_log_x)  # 反变换，确保结果严格大于 0
  
  #re-standarized
  cum_value<-my_int_function( coords , smoothed_weight_values   )$sum
  reweight_values<-smoothed_weight_values/cum_value
  
  #my_int_function( coords , reweight_values   )$sum #1 #checked
  # 查看结果
  return(  reweight_values  )
}