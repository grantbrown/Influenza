library(dplyr)
library(parallel)
library(spatialSEIR)
source("./Analysis1Functions.R")


pred.days = 120
totalSamples = 0
minimumSamples = 1000000
parameterList = buildParams(convergenceSampleSize = 10000,
                            convergenceCriterion = 1.02, 
                            extraR0Iterations = 100, 
                            iterationStride = 1000)


cl = makeCluster(3, outfile = "err.txt")
clusterExport(cl, c("parameterList"))

chains = parLapply(cl, 1:3, buildNode, nodeParams = parameterList$modelParams)

iterationParams = list(convergenceSampleSize=parameterList$runParams$convergenceSampleSize, 
                       targetAcceptanceRatio=0.2,   
                       tolerance=0.05,
                       proportionChange = 0.1,
                       updateSamplingParams = TRUE)


iterationParams = list(iterationParams, iterationParams, iterationParams)

fileNames = c("./chain_output_influenza_1.txt",
              "./chain_output_influenza_2.txt",
              "./chain_output_influenza_3.txt")
conv = FALSE
while (!conv)
{
  cat(paste("Not converged, adding iterations. Total so far: ", totalSamples, 
            "\n", sep =""))
  parLapply(cl, iterationParams, additionalIterations)
  totalSamples = totalSamples + iterationParams[[1]]$convergenceSampleSize
  conv = checkConvergence(fileNames[1], fileNames[2], fileNames[3], 
                          maxVal = parameterList$runParams$convergenceCriterion)
  conv = (conv && (minimumSamples < totalSamples))
}
cat("Chains converged, finishing up...\n")

cleanUpParamsList = list(1,2,3)
chains = parLapply(cl, cleanUpParamsList, finishSimulation)
save("chains", file="./chainOutput.Robj")
stopCluster(cl)

chain1 = chains[[1]]$chainOutput 
chain2 = chains[[2]]$chainOutput 
chain3 = chains[[3]]$chainOutput 

