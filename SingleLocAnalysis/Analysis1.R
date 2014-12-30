library(spatialSEIR)
library(dplyr)

N = as.numeric(filter(read.csv("../Data/DataProcessing/summaryPopulation.csv"), 
                    Region == "WNC")$pop)
cases = filter(select(read.csv("../Data/DataProcessing/processedData.csv"),
                      YEAR, WEEK, CASES_WNC),
               YEAR >= 2007 & YEAR <= 2010)
cases = rename(cases, CASES = CASES_WNC)

