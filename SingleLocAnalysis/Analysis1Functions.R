
buildParams = function(convergenceSampleSize=20000,
                       convergenceCriterion = 1.05, 
                       extraR0Iterations = 100, 
                       iterationStride = 1000){
  cases = filter(select(read.csv("../Data/DataProcessing/processedData.csv"),
                        YEAR, WEEK, CASES_WNC),
                 YEAR >= 2008 & YEAR <= 2009)
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
                         effectiveTransitionSampleSize=1000000
  )
  
  modelParams = list(list(seed=1901924, 
                          outFileName="./chain_output_influenza_1.txt", 
                          modelComponents=modelComponents,
                          estimateR0=TRUE),
                     list(seed=1555924, 
                          outFileName="./chain_output_influenza_2.txt", 
                          modelComponents=modelComponents,
                          estimateR0=FALSE),
                     list(seed=6612314, 
                          outFileName="./chain_output_influenza_3.txt", 
                          modelComponents=modelComponents,
                          estimateR0=FALSE))
  runParams =list(convergenceSampleSize=convergenceSampleSize,
                  convergenceCriterion=convergenceCriterion,
                  extraR0Iterations=extraR0Iterations, 
                  iterationStride=iterationStride)
  list(modelParams=modelParams,
       runParams=runParams)
}


buildNode = function(x, nodeParams=NA) 
{
  library(spatialSEIR)
  chainNumber <<- x
  if (length(nodeParams) == 1 && is.na(nodeParams)){
    stop("Where my params fool?")
  }
  nodeParams = nodeParams[[x]]
  # Store params and outFileName in global 
  # namespace of node for reference
  params <<- nodeParams
  seed = params[[1]]
  outFileName <<- params[[2]]
  modelComponents <<- params[[3]]
  
  set.seed(seed)
  
  DataModel = buildDataModel(modelComponents$I_star, 
                             type = "overdispersion", 
                             params = c(100000,1000))
  
  priorBetaIntercept = log(mean(-log(1-(modelComponents$I_star/(modelComponents$N))))) 
  ExposureModel = buildExposureModel(modelComponents$X, modelComponents$Z, 
                                     beta = c(priorBetaIntercept, 
                                              rep(0, ((length(modelComponents$beta_SE))-1))), 
                                     betaPriorPrecision = 0.1)
  ReinfectionModel = buildReinfectionModel("SEIRS", X_prs = modelComponents$X_RS, 
                                           betaPrs = -c(1, rep(0,(length(modelComponents$beta_RS)-1))), 
                                           priorMean = c(-2, 3),
                                           priorPrecision = c(100000000, 100000))
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
  
  # Shrink initial compartment guesses towards initial parameter guesses
  res$parameterSamplingMode = 17
  res$compartmentSamplingMode=16
  res$simulate(500)
  
  res$parameterSamplingMode = 8
  res$compartmentSamplingMode=1
  # Burn in tuning parameters
  for (i in 1:200)
  {
    res$simulate(10)
    res$updateSamplingParameters(0.2, 0.05, 0.01)
  }
  for (i in 1:50)
  {
    res$simulate(100)
    res$updateSamplingParameters(0.2, 0.05, 0.01)
  }
  
  
  for (i in 1:100)
  {
    res$simulate(500)
    res$updateSamplingParameters(0.2, 0.05, 0.01)
  }
  res$parameterSamplingMode = 7
  #res$compartmentSamplingMode = 17
  res$useDecorrelation = 10
  res$performHybridStep = 10
  
  # Store the model object in the global namespace of the node,
  # we can't pass these between sessions
  localModelObject <<- res
  return(list("model"=res,
              "fileName"=outFileName))    
}

additionalIterations = function(iterationParams)
{  
  localModelObject$simulate(iterationParams$convergenceSampleSize)
  if (iterationParams$updateSamplingParams)
  {
    localModelObject$updateSamplingParameters(iterationParams$targetAcceptanceRatio, 
                                              0.05,
                                              iterationParams$proportionChange)
  }
}

finishSimulation = function(iterationNumber)
{
  dat = read.csv(outFileName)
  
  ## Do we need to estimate R0 for this chain?
  if (params[["estimateR0"]])
  {  
    R0 = array(0, dim = c(nrow(modelComponents$I_star), ncol(modelComponents$I_star), parameterList$runParams$extraR0Iterations))
    effectiveR0 = array(0, dim = c(nrow(modelComponents$I_star), ncol(modelComponents$I_star), parameterList$runParams$extraR0Iterations))
    empiricalR0 = array(0, dim = c(nrow(modelComponents$I_star), ncol(modelComponents$I_star), parameterList$runParams$extraR0Iterations))
    for (i in 1:parameterList$runParams$extraR0Iterations)
    {
      localModelObject$simulate(1000)
      for (j in 0:(nrow(modelComponents$I_star)-1))
      {
        R0[j+1,,i] = localModelObject$estimateR0(j)
        effectiveR0[j+1,,i] = localModelObject$estimateEffectiveR0(j)
        empiricalR0[j+1,,i] = apply(localModelObject$getIntegratedGenerationMatrix(j), 1, sum)
      }
    }
    
    R0Mean = apply(R0, 1:2, mean)
    R0LB = apply(R0, 1:2, quantile, probs = 0.05)
    R0UB = apply(R0, 1:2, quantile, probs = 0.95)
    effectiveR0Mean = apply(effectiveR0, 1:2, mean)
    effectiveR0LB = apply(effectiveR0, 1:2, quantile, probs = 0.05)
    effectiveR0UB = apply(effectiveR0, 1:2, quantile, probs = 0.95)
    empiricalR0Mean = apply(empiricalR0, 1:2, mean)
    empiricalR0LB = apply(empiricalR0, 1:2, quantile, probs = 0.05)
    empiricalR0UB = apply(empiricalR0, 1:2, quantile, probs = 0.95)
    orig.R0 = R0
    R0 = list("R0" = list("mean"=R0Mean, "LB" = R0LB, "UB" = R0UB),
              "effectiveR0" = list("mean"=effectiveR0Mean, "LB" = effectiveR0LB, 
                                   "UB" = effectiveR0UB),
              "empiricalR0" = list("mean"=empiricalR0Mean, "LB" = empiricalR0LB, 
                                   "UB" = empiricalR0UB))
  } else
  {
    R0 = NULL
    orig.R0 = NULL
  }  
  
  return(list("chainOutput" = dat, "R0" = R0, "rawSamples" = orig.R0))
}
