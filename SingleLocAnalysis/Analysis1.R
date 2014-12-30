library(dplyr)

cases = filter(select(read.csv("../Data/DataProcessing/processedData.csv"),
                      YEAR, WEEK, CASES_WNC),
               YEAR >= 2007 & YEAR <= 2010)
cases = rename(cases, CASES = CASES_WNC)

N = as.numeric(filter(read.csv("../Data/DataProcessing/summaryPopulation.csv"), 
                      Region == "WNC")$pop)
N = matrix(N, ncol=1, nrow=nrow(cases))

timeIndex = (cases$WEEK )/52*2*pi


# Create Co-variates
X = matrix(1)
Z = cbind(sin(timeIndex), cos(timeIndex), sin(timeIndex)*cos(timeIndex))

X_prs = cbind(matrix(1, nrow = nrow(cases), ncol = 1),
              cases$YEAR == 2009 & cases$WEEK == 16) # Indicator for H1N1 introduction

modelComponents = list(I_star=matrix(cases$CASES, ncol = 1),
                       N=N,
                       X=X,
                       X_RS=X_prs,
                       Z=Z,
                       beta_SE=c(-4, rep(0, ncol(Z))),
                       beta_RS=c(-4, 0),
                       gamma_ei=4, # incubation is 1-4 days, so weekly data transition prob is almost 1
                       gamma_ir=2, # 5-7 day infectious period
                       effectiveTransitionSampleSize=1000
                       )

modelParams = list(list(seed=1901924, 
                        outFileName="./chain_output_influenza_1.txt", 
                        modelComponents=modelComponents),
                   list(seed=1555924, 
                        outFileName="./chain_output_influenza_2.txt", 
                        modelComponents=modelComponents),
                   list(seed=6612314, 
                        outFileName="./chain_output_influenza_3.txt", 
                        modelComponents=modelComponents))

buildNode = function(params) 
{
  library(spatialSEIR)
  seed = params[[1]]
  outFileName = params[[2]]
  modelComponents = params[[3]]
  
  set.seed(seed)
  
  DataModel = buildDataModel(modelComponents$I_star, 
                             type = "overdispersion", 
                             params = c(10000,10000))
  
  priorBetaIntercept = log(mean(-log(1-(modelComponents$I_star/(modelComponents$N))))) 
  ExposureModel = buildExposureModel(modelComponents$X, modelComponents$Z, 
                                     beta = c(priorBetaIntercept, 
                                              rep(0, ((length(modelComponents$beta_SE))-1))), 
                                     betaPriorPrecision = 0.1)
  ReinfectionModel = buildReinfectionModel("SEIRS", X_prs = modelComponents$X_RS, 
                                           betaPrs = -c(1, rep(0,(length(modelComponents$beta_RS)-1))), 
                                           priorPrecision = 0.1)
  SamplingControl = buildSamplingControl(iterationStride=1000,
                                         sliceWidths = c(0.26,  # S_star
                                                         0.1,  # E_star
                                                         0.15, # I_star
                                                         0.22, # S0
                                                         0.24, # I0
                                                         0.8, # beta
                                                         0.2, # betaPrs
                                                         0.015, # rho
                                                         0.01, # gamma_ei
                                                         0.01, # gamma_ir
                                                         0.01 # phi
                                         ))
  DistanceModel = buildDistanceModel(list(matrix(0)))
  TransitionPriors = buildTransitionPriorsFromProbabilities(1-exp(-modelComponents$gamma_ei), 
                                                            1-exp(-modelComponents$gamma_ir), 
                                                            modelComponents$effectiveTransitionSampleSize,
                                                            modelComponents$effectiveTransitionSampleSize) 
  
  I0 = max(modelComponents$I_star[1], modelComponents$I_star[2], 100)
  E0 = I0
  S0 = modelComponents$N[1] - I0 - E0
  InitContainer = buildInitialValueContainer(modelComponents$I_star, modelComponents$N, 
                                             S0 = S0,I0 = I0, E0 = E0, 
                                             reinfection=TRUE,dataType="I_star")
  res = buildSEIRModel(outFileName,DataModel,ExposureModel,ReinfectionModel,DistanceModel,
                       TransitionPriors, InitContainer, SamplingControl)
  
  res$setRandomSeed(seed)
  res$setTrace(0)
  res$compartmentSamplingMode = 17
  res$useDecorrelation = 50
  res$performHybridStep = 50
  # Burn in tuning parameters
  for (i in 1:(200))
  {
    res$simulate(10)
    res$updateSamplingParameters(0.2, 0.05, 0.01)
  }
  for (i in 1:(20))
  {
    res$simulate(100)
    res$updateSamplingParameters(0.2, 0.05, 0.01)
  }
  
  # Store the model object in the global namespace of the node,
  # we can't pass these between sessions
  localModelObject <<- res
  return(list("model"=res,
              "fileName"=outFileName))    
}




