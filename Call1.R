library(coda)
library(rjags)
library(rJava)

library(xlsx)
library(HDInterval)
load.module("mix")


snSMART_PLH_cont_sim <- function(EFFECT_STAGE_I=c(-75,-75,-75),
                                 EFFECT_STAGE_I_SD=25,
                                 SAMPLE_SIZE_TOTAL=60,
                                 SIMULATION_N=5,
                                 ALPHA=0, 
                                 BETA=1, 
                                 bugpath=c("JointStageBayes_mixture.bug"),
                                 firststage_bugpath=c("FirstStageBayes_mixture.bug"),
                                 sourcefile=c("Simulation codes.R"),
                                 outfile=c("test.xlsx")
){
  ###############################################################
  
  # Parameters settings
  EFFECT_STAGE_I_P=EFFECT_STAGE_I[1]
  EFFECT_STAGE_I_L=EFFECT_STAGE_I[2]
  EFFECT_STAGE_I_H=EFFECT_STAGE_I[3]
  
  
  SAMPLE_SIZE_ARM_P=SAMPLE_SIZE_TOTAL/3
  SAMPLE_SIZE_ARM_L=SAMPLE_SIZE_TOTAL/3
  SAMPLE_SIZE_ARM_H=SAMPLE_SIZE_TOTAL/3
  
  ###############################################################
  COVERAGE_RATE=0.95
  NUM_ARMS=3
  
  MCMC_SAMPLE=20000
  BURNING=10000
  NUM_MCMC_CHAIN=2
  
  ###############################################################
  # source simulation
  source(sourcefile,local=TRUE, encoding = "UTF-8",echo = T)
  
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
  
  
  
  save.xlsx(outfile,
            
            round(bias_effect_posterior_mean,2),
            round(rmse_effect_posterior_mean,2),
            
            round(firstbias_effect_posterior_mean,2),
            round(firstrmse_effect_posterior_mean,2),
            
            round(bias_discount_parameter_posterior_mean,2),
            round(rmse_discount_parameter_posterior_mean,2),
            ###############################################################
            round(cr_response_rate,2),
            round(hdi_length_response_rate,2),
            round(firstcr_response_rate,2),
            round(firsthdi_length_response_rate,2),
            ###############################################################
            error_round,
            warn_round,
            error_count,
            warn_count,
            sim_count
  )
}

snSMART_PLH_cont_sim(
  SAMPLE_SIZE_TOTAL=60,
  EFFECT_STAGE_I=c(-60,-60,-60),
  EFFECT_STAGE_I_SD=25,
  SIMULATION_N=2000, 
  ALPHA=0, 
  BETA=1, 
  bugpath=c("/home/sqrfang/ContsnSMART/RerunZ=0/JointStageBayes_mixture.bug"), # Mixture priors for joint stage model
  firststage_bugpath=c("/home/sqrfang/ContsnSMART/RerunZ=0/FirstStageBayes_mixture.bug"), # Mixture priors for first stage model
  sourcefile=c("/home/sqrfang/ContsnSMART/RerunZ=0/Simulation codes.R"),
  outfile=c("/home/sqrfang/ContsnSMART/RerunZ=0/N=60_m60_mix.xlsx")
)

