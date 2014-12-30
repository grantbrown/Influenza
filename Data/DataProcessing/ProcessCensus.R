library(dplyr)
dat = read.csv("../Census/PopulationCategorized.csv")

dat_filter = filter(dat, !is.na(Region))
pop2010 = summarize(group_by(dat_filter, Region), pop = sum(Census2010))
write.csv(pop2010, file="summaryPopulation.csv", row.names=FALSE)
