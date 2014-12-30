source("./ProcessCensus.R")
# ILI Data
dat = read.csv("../FluViewPhase2Data/ILINet.csv", skip=1)
# Laboratory Data
labDat = read.csv("../FluViewPhase2Data/WHO_NREVSS.csv")

f = function(x){
  ifelse(x == "East North Central", "ENC",
  ifelse(x == "East South Central", "ESC",
  ifelse(x == "Mid-Atlantic", "MAT",
  ifelse(x == "Mountain", "MTN",        
  ifelse(x == "New England", "NEW",
  ifelse(x == "Pacific", "PAC",       
  ifelse(x == "South Atlantic", "SAT",
  ifelse(x == "West North Central", "WNC",
  ifelse(x == "West South Central", "WSC", "UNK")))))))))
}


dat = filter(group_by(select(mutate(dat, Region = f(REGION), 
                    ILI_PERCENT = X.UNWEIGHTED.ILI), 
                    Region, YEAR, WEEK, ILI_PERCENT), Region), !is.na(ILI_PERCENT))

labDat = filter(group_by(select(mutate(labDat, Region = f(REGION), 
                             INFLZA_PERCENT = PERCENT.POSITIVE), 
                      Region, YEAR, WEEK, INFLZA_PERCENT), Region), !is.na(INFLZA_PERCENT))

smoother = function(x){
  ksmooth(1:length(x), x, kernel = "normal", bandwidth=10)$y
}
uqLoc = unique(labDat$Region)
a = labDat$INFLZA_PERCENT
for (loc in uqLoc){  
  labDat$INFLZA_PERCENT[labDat$Region == loc] = 
    smoother(labDat$INFLZA_PERCENT[labDat$Region == loc])
}

dat = left_join(dat, labDat, by=c("Region", "YEAR", "WEEK"))
# take previous INFLZA_PERCENT if missing
dat$INFLZA_PERCENT[is.na(dat$INFLZA_PERCENT)] = 
  c(0, dat$INFLZA_PERCENT[1:(nrow(dat)-1)])[is.na(dat$INFLZA_PERCENT)]

pop2010 = group_by(pop2010, Region)
dat = filter(
        mutate(
          left_join(dat, pop2010, by = "Region"), 
               CASES = round((ILI_PERCENT/100)*(INFLZA_PERCENT/100)*pop)),
             YEAR > 2002
        )
dat = select(dat, Region, YEAR, WEEK, CASES)
f = function(region, data){
  data[data$Region == region,]
}
uqLoc = unique(dat$Region)
subsets = lapply(uqLoc, f, data = dat)

outDat = subsets[[1]]
outDat[[paste("CASES", uqLoc[1], sep = "_")]] = outDat$CASES
outDat = ungroup(outDat)
outDat = outDat[,which(names(outDat) != "CASES" & names(outDat) != "Region")]

for (i in 2:length(subsets)){
  subsets[[i]][[paste("CASES", uqLoc[i], sep ="_")]] = subsets[[i]]$CASES
  subsets[[i]] = ungroup(subsets[[i]])[,(names(subsets[[i]]) != "CASES" &
                                names(subsets[[i]]) != "Region")]
  outDat = left_join(outDat, subsets[[i]], by=c("YEAR", "WEEK"))
}
write.csv(outDat, "processedData.csv", row.names = FALSE)



