


### simulation study to address Comment (B1) in revision:
### Compare eyeballing stratum-specific effects to the rigor SSS inference


devtools::install_github("HDTian/SSS")
library(SSS)

### Simulation
N <- 50000
set.seed(100)
Z <- rbinom( N , 1 , 0.5 )
U <- rnorm( N , 0 , 1 )
Ex <- rnorm( N , 0 , 1 )
alphaZ<- 0.15
X <-  alphaZ*Z + U + Ex
summary( lm(  X  ~ Z )  )$r.squared
Ey <- rnorm( N , 0 , 1 )
Y <- 1.0*X  - 2.0*(X-0)*(X>0)  + U + Ey  # one change-point located at x=-Inf and x=0


### SSS analysis
SSS_res <- SSS(Z,X,Y,Ns_used=10,x_baseline_used = 0,
               tpoints_used = 'quantile', boundary_used=c(0.05,0.95) )


### change-pooint PIP plot
ggdata<-SSS_res$PIP_ggdata
p1<-ggplot(ggdata, aes(x = positions, y = PIP)) +
  #geom_point()+
  geom_line() +
  #geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs", k = 60), color='black',se = FALSE) +
  geom_vline(aes(xintercept = posterior_mode), linetype = "dashed") +
  #geom_vline(aes(xintercept = posterior_mean), linetype = "dotted") +
  labs(x = "Exposure level", y = "Posterior inclusion probability") +
  facet_grid( KK ~. ,scales = "free_y" ) +
  scale_x_continuous(breaks = extended_breaks(n = 9)) +  # 让刻度尽量均匀
  theme_bw() +theme(legend.position = "none")
p1


### LACE plot
stratum_summary<-SSS_res$stratum_summary
stratum_summary$MR_est<-stratum_summary$IV_outcome_est/stratum_summary$IV_exposure_est
stratum_summary$MR_se <- stratum_summary$IV_outcome_se/stratum_summary$IV_exposure_est
p2<-ggplot(stratum_summary, aes(exposure_mean, MR_est)) +
  annotate("segment",x = -Inf, xend = 0,y = 1, yend = 1, colour = "grey70", linewidth = 0.6) +
  annotate("segment",x = 0, xend = Inf,y = -1, yend = -1, colour = "grey70", linewidth = 0.6) +
  annotate("segment",x = 0, xend = 0,y = -1, yend = 1, colour = "grey70", linewidth = 0.6) +
  #geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
  geom_pointrange(aes(ymin = MR_est - 1.96*MR_se, ymax = MR_est + 1.96*MR_se)) +
  labs(x = "Exposure level", y = "Stratum-specific IV estimate")+theme_bw()
p2


### joint plot
library(grid)
grid.newpage()
pushViewport(viewport(layout = grid.layout(60 , 100 )))
print(p2, vp = viewport(layout.pos.row = 1:60, layout.pos.col = 2:54))
print(p1, vp = viewport(layout.pos.row = 1:60, layout.pos.col = 55:99))
#suggested size: 830 360

