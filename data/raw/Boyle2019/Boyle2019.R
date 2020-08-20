# BOYLE 2019 for meta-analysis
# reprocessing data for a synthesis of plant-associate datasets
# jby 2019.12.16

# starting up ------------------------------------------------------------
# setwd("~/Documents/Academic/Active_projects/plant-associate-lit-search/data/Boyle2019")

require("tidyverse")

source("~/Documents/Academic/Active_projects/shared/Rscripts/base.R") # my special mix of personal functions
source("~/Documents/Academic/Active_projects/shared/Rscripts/base_graphics.R")

require(fields)


#-------------------------------------------------------------------------
# read and prep in data

sites <- read.csv("~/Documents/Academic/Active_projects/plant-associate-lit-search/data/Boyle2019/Boyle2019_geo.txt", h=T, sep="\t")

sites

regions <- sites %>% group_by(Cluster) %>% summarize(mnLon = mean(Lon), mnLat = mean(Lat))

regions

# calculate GS distances between sites
Gdists <- rdist.earth(as.matrix(regions[,c("mnLon","mnLat")]),as.matrix(regions[,c("mnLon","mnLat")]), miles=FALSE) #guhhh
colnames(Gdists) <- rownames(Gdists) <- regions$Cluster
Gdists[Gdists<1e-3] <- 0
Gdists <- as.data.frame(Gdists)%>% gather("To", "Gdist") %>% mutate(From = rep(rownames(Gdists), ncol(Gdists))) %>% filter(!is.na(Gdist)) %>% select(From, To, Gdist) -> Gdists

Gdists

write.table(Gdists, "~/Documents/Academic/Active_projects/plant-associate-lit-search/data/Boyle2019/geo-dists.txt", sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)


