
###############################################################

# Initialization

error_round=0
warn_round=0

error_count=0
warn_count=0

firsterror_round=0
firstwarn_round=0

firsterror_count=0
firstwarn_count=0

effect_posterior_mean=NULL
discount_parameter_posterior_mean=NULL
response_rate_hdi_coverage_rate=NULL
response_rate_hdi_length=NULL

firstresponse_rate_hdi_coverage_rate=NULL
firstresponse_rate_hdi_length=NULL


firsteffect_posterior_mean=NULL


sim_count=0
set.seed(11)



for(i in 1:SIMULATION_N){
  # generate data for the simulation
  sim_count=sim_count+1
  
  #T1
  DATA_TREAT_STAGE_I<-c(rep(1,SAMPLE_SIZE_ARM_P),rep(2,SAMPLE_SIZE_ARM_L),rep(3,SAMPLE_SIZE_ARM_H))
  
  #Y1
  P <- rnorm(SAMPLE_SIZE_ARM_P, EFFECT_STAGE_I_P, EFFECT_STAGE_I_SD)
  L <- rnorm(SAMPLE_SIZE_ARM_L, EFFECT_STAGE_I_L, EFFECT_STAGE_I_SD)
  H <- rnorm(SAMPLE_SIZE_ARM_H, EFFECT_STAGE_I_H, EFFECT_STAGE_I_SD)
  
  Y1 <- c(P, L,H)
  
  #Intermediate response that determines the rerandomization in stage 2
  Z <- Y1>0
  
  #T2
  DATA_TREAT_STAGE_II <- rbinom(SAMPLE_SIZE_TOTAL,1,0.5) +2
  DATA_TREAT_STAGE_II[Z==FALSE & DATA_TREAT_STAGE_I==3] <- 3
  
  #Y2
  mu2 <- EFFECT_STAGE_I[DATA_TREAT_STAGE_II] + ALPHA[DATA_TREAT_STAGE_I] + BETA[DATA_TREAT_STAGE_I] * (Y1-EFFECT_STAGE_I[DATA_TREAT_STAGE_I])
  Y2 <- rnorm(n=SAMPLE_SIZE_TOTAL,mean=mu2, sd=EFFECT_STAGE_I_SD)  
  
  # First stage model
  tryCatch({
    firstjags <- jags.model(firststage_bugpath,
                            data = list(overall_sample_size=SAMPLE_SIZE_TOTAL,
                                        data_effect_stageI=Y1,
                                        treatment_stageI=DATA_TREAT_STAGE_I),
                            n.chains = NUM_MCMC_CHAIN) 
    firstposterior_sample <- coda.samples(firstjags,
                                          c('mu'),
                                          #c('mu'),
                                          MCMC_SAMPLE)
  },
  firsterror = function(c) {rbind(firsterror_round,i)
    firstposterior_sample_burn<<-window(firstposterior_sample,start=BURNING, end=MCMC_SAMPLE)
    firstposterior_sample_cmb<<-do.call(rbind, firstposterior_sample_burn)
    firsterror_round<<-rbind(firsterror_round,i)
    firsterror_count<<-firsterror_count+1
  },
  firstwarning = function(c) {firstwarn_round<<-rbind(firstwarn_round,i)
  firstwarn_count<<-firstwarn_count+1},
  finally = {  
    firstposterior_sample_burn<-window(firstposterior_sample,start=BURNING, end=MCMC_SAMPLE)
    firstposterior_sample_cmb<-do.call(rbind, firstposterior_sample_burn)
  }
  )
  
  
  ###############################################################  
  
  firstposterior_sample_cmb <- cbind(firstposterior_sample_cmb, firstposterior_sample_cmb[,2]-firstposterior_sample_cmb[,1])
  firstposterior_sample_cmb <- cbind(firstposterior_sample_cmb, firstposterior_sample_cmb[,3]-firstposterior_sample_cmb[,1])
  colnames(firstposterior_sample_cmb) <- c("mu[1]", "mu[2]", "mu[3]", "mu[4]","mu[5]")
  
  tmp_effect_posterior_mean=colMeans(firstposterior_sample_cmb[,c("mu[1]","mu[2]","mu[3]", "mu[4]","mu[5]")])
  
  ###############################################################  
  tmp_hdi=hdi(firstposterior_sample_cmb, COVERAGE_RATE)
  
  tmp_response_rate_hdi_coverage_rate_mu1=as.numeric(tmp_hdi["lower","mu[1]"]<=EFFECT_STAGE_I_P & tmp_hdi["upper","mu[1]"]>=EFFECT_STAGE_I_P)
  tmp_response_rate_hdi_coverage_rate_mu2=as.numeric(tmp_hdi["lower","mu[2]"]<=EFFECT_STAGE_I_L & tmp_hdi["upper","mu[2]"]>=EFFECT_STAGE_I_L)
  tmp_response_rate_hdi_coverage_rate_mu3=as.numeric(tmp_hdi["lower","mu[3]"]<=EFFECT_STAGE_I_H & tmp_hdi["upper","mu[3]"]>=EFFECT_STAGE_I_H)
  tmp_response_rate_hdi_coverage_rate_mu4=as.numeric(tmp_hdi["lower","mu[4]"]<=(EFFECT_STAGE_I_L-EFFECT_STAGE_I_P) & tmp_hdi["upper","mu[4]"]>=(EFFECT_STAGE_I_L-EFFECT_STAGE_I_P))
  tmp_response_rate_hdi_coverage_rate_mu5=as.numeric(tmp_hdi["lower","mu[5]"]<=(EFFECT_STAGE_I_H-EFFECT_STAGE_I_P) & tmp_hdi["upper","mu[5]"]>=(EFFECT_STAGE_I_H-EFFECT_STAGE_I_P))
  tmp_response_rate_hdi_coverage_rate=cbind(tmp_response_rate_hdi_coverage_rate_mu1,tmp_response_rate_hdi_coverage_rate_mu2,tmp_response_rate_hdi_coverage_rate_mu3,
                                            tmp_response_rate_hdi_coverage_rate_mu4, tmp_response_rate_hdi_coverage_rate_mu5)
  firstresponse_rate_hdi_coverage_rate=rbind(firstresponse_rate_hdi_coverage_rate,tmp_response_rate_hdi_coverage_rate)
  
  tmp_response_rate_hdi_length_mu1=abs(tmp_hdi["lower","mu[1]"]-tmp_hdi["upper","mu[1]"])
  tmp_response_rate_hdi_length_mu2=abs(tmp_hdi["lower","mu[2]"]-tmp_hdi["upper","mu[2]"])
  tmp_response_rate_hdi_length_mu3=abs(tmp_hdi["lower","mu[3]"]-tmp_hdi["upper","mu[3]"])
  tmp_response_rate_hdi_length_mu4=abs(tmp_hdi["lower","mu[4]"]-tmp_hdi["upper","mu[4]"])
  tmp_response_rate_hdi_length_mu5=abs(tmp_hdi["lower","mu[5]"]-tmp_hdi["upper","mu[5]"])
  tmp_response_rate_hdi_length=cbind(tmp_response_rate_hdi_length_mu1,tmp_response_rate_hdi_length_mu2,tmp_response_rate_hdi_length_mu3,
                                     tmp_response_rate_hdi_length_mu4,tmp_response_rate_hdi_length_mu5)
  firstresponse_rate_hdi_length=rbind(firstresponse_rate_hdi_length,tmp_response_rate_hdi_length)
  
  ###############################################################
  firsteffect_posterior_mean=rbind(firsteffect_posterior_mean,tmp_effect_posterior_mean)
  
  
  # Joint stage model
  
  tryCatch({
    jags <- jags.model(bugpath,
                       data = list(overall_sample_size=SAMPLE_SIZE_TOTAL,
                                   data_effect_stageI=Y1,
                                   data_effect_stageII=Y2,
                                   treatment_stageI=DATA_TREAT_STAGE_I,
                                   treatment_stageII=DATA_TREAT_STAGE_II
                       ),
                       n.chains = NUM_MCMC_CHAIN) 
    posterior_sample <- coda.samples(jags,
                                     c('mu','alpha','beta'),
                                     MCMC_SAMPLE)
  },
  error = function(c) {rbind(error_round,i)
    posterior_sample_burn<<-window(posterior_sample,start=BURNING, end=MCMC_SAMPLE)
    posterior_sample_cmb<<-do.call(rbind, posterior_sample_burn)
    error_round<<-rbind(error_round,i)
    error_count<<-error_count+1
  },
  warning = function(c) {warn_round<<-rbind(warn_round,i)
  warn_count<<-warn_count+1},
  finally = {  
    posterior_sample_burn<-window(posterior_sample,start=BURNING, end=MCMC_SAMPLE)
    posterior_sample_cmb<-do.call(rbind, posterior_sample_burn)
  }
  )
  #####################################################################################
  
  
  posterior_sample_cmb <- cbind(posterior_sample_cmb, posterior_sample_cmb[,c("mu[2]")]-posterior_sample_cmb[,c("mu[1]")])
  posterior_sample_cmb <- cbind(posterior_sample_cmb, posterior_sample_cmb[,c("mu[3]")]-posterior_sample_cmb[,c("mu[1]")])
  colnames(posterior_sample_cmb) <- c("alpha","beta","mu[1]", "mu[2]", "mu[3]", "mu[4]","mu[5]")
  
  tmp_effect_posterior_mean=colMeans(posterior_sample_cmb[,c("mu[1]","mu[2]","mu[3]", "mu[4]","mu[5]")])
  
  tmp_effect_posterior_mean[4] = mean(posterior_sample_cmb[,c("mu[2]")]-posterior_sample_cmb[,c("mu[1]")])
  tmp_effect_posterior_mean[5] = mean(posterior_sample_cmb[,c("mu[3]")]-posterior_sample_cmb[,c("mu[1]")])
  
  tmp_discount_parameter_posterior_mean=colMeans(posterior_sample_cmb[,c("alpha","beta")])
  
  ###############################################################  
  tmp_hdi=hdi(posterior_sample_cmb, COVERAGE_RATE)
  
  tmp_response_rate_hdi_coverage_rate_mu1=as.numeric(tmp_hdi["lower","mu[1]"]<=EFFECT_STAGE_I_P & tmp_hdi["upper","mu[1]"]>=EFFECT_STAGE_I_P)
  tmp_response_rate_hdi_coverage_rate_mu2=as.numeric(tmp_hdi["lower","mu[2]"]<=EFFECT_STAGE_I_L & tmp_hdi["upper","mu[2]"]>=EFFECT_STAGE_I_L)
  tmp_response_rate_hdi_coverage_rate_mu3=as.numeric(tmp_hdi["lower","mu[3]"]<=EFFECT_STAGE_I_H & tmp_hdi["upper","mu[3]"]>=EFFECT_STAGE_I_H)
  tmp_response_rate_hdi_coverage_rate_mu4=as.numeric(tmp_hdi["lower","mu[4]"]<=(EFFECT_STAGE_I_L-EFFECT_STAGE_I_P) & tmp_hdi["upper","mu[4]"]>=(EFFECT_STAGE_I_L-EFFECT_STAGE_I_P))
  tmp_response_rate_hdi_coverage_rate_mu5=as.numeric(tmp_hdi["lower","mu[5]"]<=(EFFECT_STAGE_I_H-EFFECT_STAGE_I_P) & tmp_hdi["upper","mu[5]"]>=(EFFECT_STAGE_I_H-EFFECT_STAGE_I_P))
  tmp_response_rate_hdi_coverage_rate=cbind(tmp_response_rate_hdi_coverage_rate_mu1,tmp_response_rate_hdi_coverage_rate_mu2,tmp_response_rate_hdi_coverage_rate_mu3,
                                            tmp_response_rate_hdi_coverage_rate_mu4, tmp_response_rate_hdi_coverage_rate_mu5)
  response_rate_hdi_coverage_rate=rbind(response_rate_hdi_coverage_rate,tmp_response_rate_hdi_coverage_rate)
  
  tmp_response_rate_hdi_length_mu1=abs(tmp_hdi["lower","mu[1]"]-tmp_hdi["upper","mu[1]"])
  tmp_response_rate_hdi_length_mu2=abs(tmp_hdi["lower","mu[2]"]-tmp_hdi["upper","mu[2]"])
  tmp_response_rate_hdi_length_mu3=abs(tmp_hdi["lower","mu[3]"]-tmp_hdi["upper","mu[3]"])
  tmp_response_rate_hdi_length_mu4=abs(tmp_hdi["lower","mu[4]"]-tmp_hdi["upper","mu[4]"])
  tmp_response_rate_hdi_length_mu5=abs(tmp_hdi["lower","mu[5]"]-tmp_hdi["upper","mu[5]"])
  tmp_response_rate_hdi_length=cbind(tmp_response_rate_hdi_length_mu1,tmp_response_rate_hdi_length_mu2,tmp_response_rate_hdi_length_mu3,
                                     tmp_response_rate_hdi_length_mu4, tmp_response_rate_hdi_length_mu5)
  response_rate_hdi_length=rbind(response_rate_hdi_length,tmp_response_rate_hdi_length)
  
  ###############################################################
  
  effect_posterior_mean=rbind(effect_posterior_mean,tmp_effect_posterior_mean)
  discount_parameter_posterior_mean=rbind(discount_parameter_posterior_mean,tmp_discount_parameter_posterior_mean)
  
}
# Mu 
########### First stage model Bias ########### 
firstbias_effect_posterior_mean=colMeans(firsteffect_posterior_mean)-c(EFFECT_STAGE_I_P,EFFECT_STAGE_I_L,EFFECT_STAGE_I_H,EFFECT_STAGE_I_L-EFFECT_STAGE_I_P,EFFECT_STAGE_I_H-EFFECT_STAGE_I_P)

########### First stage model rMSE ########### 
firstrmse_effect_posterior_mean=sqrt(firstbias_effect_posterior_mean^2+diag(var(firsteffect_posterior_mean)))

########### Joint stage model Bias ########### 
bias_effect_posterior_mean=colMeans(effect_posterior_mean)-c(EFFECT_STAGE_I_P,EFFECT_STAGE_I_L,EFFECT_STAGE_I_H, EFFECT_STAGE_I_L-EFFECT_STAGE_I_P, EFFECT_STAGE_I_H-EFFECT_STAGE_I_P)

########### Joint stage model rMSE ########### 
rmse_effect_posterior_mean=sqrt(bias_effect_posterior_mean^2+diag(var(effect_posterior_mean)))

#Alpha and Beta

bias_discount_parameter_posterior_mean=rep(colMeans(discount_parameter_posterior_mean),c(3,3))-c(ALPHA, BETA)
rmse_discount_parameter_posterior_mean=sqrt(bias_discount_parameter_posterior_mean^2+diag(var(discount_parameter_posterior_mean)))


## COVERAGE_RATE
###############################################################
# HDI
cr_response_rate=colMeans(response_rate_hdi_coverage_rate)
firstcr_response_rate=colMeans(firstresponse_rate_hdi_coverage_rate)

###############################################################

###############################################################
## HDI length
hdi_length_response_rate=colMeans(response_rate_hdi_length)                               
firsthdi_length_response_rate=colMeans(firstresponse_rate_hdi_length)  

###############################################################



# save posterior mean from simulation on disk
#save the result list to an excel, each item as a workbook
save.xlsx <- function (file, ...)
{
  require(xlsx, quietly = TRUE)
  objects <- list(...)
  fargs <- as.list(match.call(expand.dots = TRUE))
  objnames <- as.character(fargs)[-c(1, 2)]
  nobjects <- length(objects)
  for (i in 1:nobjects) {
    if (i == 1)
      write.xlsx(objects[[i]], file, sheetName = objnames[i])
    else write.xlsx(objects[[i]], file, sheetName = objnames[i],
                    append = TRUE)
  }
  print(paste("Workbook", file, "has", nobjects, "worksheets."))
}


