temp = read.csv("../Data/Temperature/output3.csv", head=FALSE,as.is=TRUE)

is.leap.year = function(years){
  ifelse((years %% 4 != 0), FALSE, 
  ifelse((years %% 100 != 0), TRUE, 
  ifelse((years %% 400 != 0), FALSE, TRUE)))
}

tagStateAbbreviations = function(abbrev){
  NEW = c("CT","RI","MA","NH","VT","ME")
  MAT = c("PA","NY","NJ")
  ENC = c("WI","IL","IN","MI","OH")
  WNC = c("ND","SD","NE","KS","IA", "MO", "MN")
  SAT = c("FL","GA","SC","NC","VA","WV","DC","MD","DE")
  ESC = c("KY","TN","MS","AL")
  WSC = c("TX","OK","AR","LA")
  MTN = c("ID", "MT", "WY", "UT", "CO", "NM","AZ","NV")
  PAC = c("WA", "OR", "CA")
  
  

  as.character(ifelse(abbrev %in% MTN,"MTN",
  ifelse(abbrev %in% WNC,"WNC",
  ifelse(abbrev %in% ENC,"ENC",
  ifelse(abbrev %in% WSC,"WSC",
  ifelse(abbrev %in% ESC,"ESC",
  ifelse(abbrev %in% SAT,"SAT",
  ifelse(abbrev %in% MAT,"MAT",
  ifelse(abbrev %in% NEW,"NEW",
  ifelse(abbrev %in% PAC,"PAC","None"
  ))))))))))
}


categories = as.character(sapply(temp[,1], tagStateAbbreviations))
temp[,1] = categories
names(temp) = c("region", "year", "month", "temp")
monthlyTemp = summarize(group_by(temp, region, year, month), avgTemp = mean(temp))


buildParams = function(convergenceSampleSize=20000,
                       convergenceCriterion = 1.05, 
                       extraR0Iterations = 100, 
                       iterationStride = 5000){
  selectedYears = c(2007, 2008,2009,2010)
  singleLocation = TRUE
  cases = filter(select(read.csv("../Data/DataProcessing/processedDataNoSmooth.csv"),
                         -timeIndex),
                 YEAR %in% selectedYears)
  monthlyTemp.sub = monthlyTemp[monthlyTemp$year %in% selectedYears,]
  includedYears = unique(cases$YEAR)
  leaps = is.leap.year(includedYears)
  
  dayWeights = c()
  for (i in 1:length(includedYears)){
    if (leaps[i]){
      dayWeights = c(dayWeights, 31,29,31,30,31,30,31,31,30,31,30,31)
    }
    else{
      dayWeights = c(dayWeights, 31,28,31,30,31,30,31,31,30,31,30,31)
    }
  }
  
  tmpMat = matrix(monthlyTemp.sub$avgTemp, nrow = (nrow(monthlyTemp.sub)/length(unique( monthlyTemp.sub$region))),
                  ncol = length(unique(monthlyTemp.sub$region)))
  colnames(tmpMat) = unique(monthlyTemp.sub$region)
  
  #avgTemp = monthlyTemp.sub[monthlyTemp.sub$region == "ENC",][1:12,]$avgTemp
  idx=unlist(sapply(1:nrow(tmpMat), function(x){
    as.numeric(rep((1:nrow(tmpMat))[x], dayWeights[x]))}))
  tmpMat = tmpMat[idx,]
  facData = select(summarize(group_by(data.frame(tmpMat, factorA = 
                         cut(1:nrow(tmpMat), 
                             breaks = sum(leaps*53 + (1-leaps)*52), 
                             include.lowest = TRUE)),
                     factorA), 
                     ENC_Temp = mean(ENC),
                               ESC_Temp = mean(ESC),
                               MAT_Temp = mean(MAT),
                               MTN_Temp = mean(MTN),
                               NEW_Temp = mean(NEW),
                               PAC_Temp = mean(PAC),
                               SAT_Temp = mean(SAT),
                               WNC_Temp = mean(WNC),
                               WSC_Temp = mean(WSC)), - factorA)
  
  facData.standard = (facData - mean(as.matrix(facData)))/sd(as.matrix(facData))
    
  N = matrix(100000, nrow=nrow(cases), ncol = ncol(facData.standard))
  
  
  timeIndex = 1:nrow(cases)
  #trigBasis1 = sin((timeIndex + 10)/52*2*pi)
  
  # Create Co-variates
  X = matrix(1, nrow = ncol(facData.standard), ncol = 1)
  Z = matrix(as.numeric(as.matrix(facData.standard)), ncol = 1)
  Z = cbind(Z, rep((timeIndex - mean(timeIndex))/sd(timeIndex), nrow(X)))
  
  H1N1 = exp(-((126):nrow(cases) - 126)/7)
  H1N1 = c(rep(0, nrow(cases) - length(H1N1)), H1N1)

  X_prs = sin((seq(1, nrow(cases)) + 13)/52*2*pi)
  for (i in 2:4){
    X_prs = cbind(X_prs, sin((seq(1, nrow(cases)) + 13*i)/52*2*pi))
  }
  X_prs = cbind(1, H1N1, X_prs)
  
  I_star = as.matrix(cases[,(3:11)[order(names(cases)[3:11])]])
  
  
  DM1 = (1-diag(ncol(I_star)))/ncol(I_star)
  DM2 = matrix(0, ncol = ncol(I_star), nrow = ncol(I_star))
  DM2[1, 6] = 1 # PAC-MTN
  DM2[4, 8] = 1 # MTN-WNC
  DM2[4, 9] = 1 # MTN-WSC
  DM2[8, 9] = 1 # WNC-WSC
  DM2[1, 8] = 1 # ENC-WNC
  DM2[2, 9] = 1 # ESC-WSC
  DM2[1, 2] = 1 # ENC-ESC
  DM2[2, 7] = 1 # ESC-SAT
  DM2[1, 7] = 1 # ENC-SAT
  DM2[1, 3] = 1 # ENC-MAT
  DM2[3, 7] = 1 # SAT-MAT
  DM2[3, 5] = 1 # MAT-NEW
  DM2 = DM2 + t(DM2)
  dmList = list(DM1, DM2)
  
  modelComponents = list(I_star=I_star,
                         N=N,
                         X=X,
                         X_RS=X_prs,
                         Z=Z,
                         dmList=dmList,
                         beta_SE=c(-4, rep(0.1, ncol(Z))),
                         beta_RS=rep(0.1, ncol(X_prs)),
                         gamma_ei=1, # incubation is 1-4 days, so weekly data transition prob is almost 1
                         gamma_ir=0.5, # 5-7 day infectious period
                         singleLocation=singleLocation,
                         effectiveTransitionSampleSize=20000
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
  
  if (!modelComponents$singleLocation){
    stop("Spatial Model Disabled")
  }
  else{    
    DataModel = buildDataModel(modelComponents$I_star[,8], 
                               type = "overdispersion",
			       phi=2)

    priorBetaIntercept = log(mean(-log(1-(modelComponents$I_star[,8]/(modelComponents$N[,8]))))) 
    ExposureModel = buildExposureModel_depricated(modelComponents$X[8,,drop=FALSE], 
                                                  modelComponents$Z[(7*(nrow(modelComponents$I_star))+1):(8*nrow(modelComponents$I_star)),], 
                                                  beta = c(priorBetaIntercept, 
                                                           rep(0.1, ((length(modelComponents$beta_SE))-1))), 
                                                  betaPriorPrecision = rep(1, length(modelComponents$beta_SE)), 
                                                  betaPriorMean = c(priorBetaIntercept, 
                                                                    rep(0.0, (length(modelComponents$beta_SE))-1)),
                                                  offset=rep(7, nrow(modelComponents$I_star)))
    ReinfectionModel = buildReinfectionModel("SEIRS", X_prs = modelComponents$X_RS, 
                                             betaPrs = c(-4,2,rep(0.1,(ncol(modelComponents$X_RS)-2))), 
                                             priorMean =  c(-4,2, rep(0.1,(ncol(modelComponents$X_RS) - 2))),
                                             priorPrecision =  c(4, rep(1,(ncol(modelComponents$X_RS)-1))))
    SamplingControl = buildSamplingControl(iterationStride=5000,
                                           sliceWidths = c(0.2,  # S_star
                                                           0.2,  # E_star
                                                           0.2, # I_star
                                                           0.2, # S0
                                                           0.2, # I0
                                                           0.01, # beta
                                                           0.01, # betaPrs
                                                           0.01, # rho
                                                           0.01, # gamma_ei
                                                           0.01, # gamma_ir
                                                           0.05 # phi
                                           ))
    DistanceModel = buildDistanceModel(list(matrix(0)))
    TransitionPriors = buildTransitionPriorsFromProbabilities(1-exp(-modelComponents$gamma_ei), 
                                                              1-exp(-modelComponents$gamma_ir), 
                                                              modelComponents$effectiveTransitionSampleSize,
                                                              modelComponents$effectiveTransitionSampleSize) 
    
    I0 = max(modelComponents$I_star[1:3,8])
    E0 = I0
    S0 = modelComponents$N[1,8] - I0 - E0
    InitContainer = buildInitialValueContainer(modelComponents$I_star[,8], modelComponents$N[,8], 
                                               S0 = S0, I0 = I0, E0 = E0, 
                                               reinfection=TRUE,dataType="I_star")
    res = buildSEIRModel(outFileName,DataModel,ExposureModel,ReinfectionModel,DistanceModel,
                         TransitionPriors, InitContainer, SamplingControl)
    
    res$setRandomSeed(seed + 1)
    res$setTrace(0)
    
    # Burn in tuning parameters
    #for (i in 1:1000)
    #{
    #  res$simulate(20)
    #  res$updateSamplingParameters(0.1, 0.01, 0.01)
    #}
    #for (i in 1:1000)
    #{
    #  res$simulate(100)
    #  res$updateSamplingParameters(0.05, 0.001, 0.001)
    #}
    
    #res$parameterSamplingMode = 8

    res$compartmentSamplingMode = 17
    res$useDecorrelation = 20
    res$performHybridStep = 100
    res$simulate(10000)
    # Store the model object in the global namespace of the node,
    # we can't pass these between sessions
    localModelObject <<- res
    return(list("model"=res,
                "fileName"=outFileName)) 
  }
}

additionalIterations = function(iterationParams)
{  
  localModelObject$simulate(iterationParams$convergenceSampleSize)
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
