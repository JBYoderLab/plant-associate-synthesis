# GARRIDO et al 2012 for meta-analysis
# reprocessing data for a synthesis of plant-associate datasets
# jby 2019.12.10

# starting up ------------------------------------------------------------
# setwd("~/Documents/Academic/Active_projects/plant-associate-lit-search/data/Yoder2013")

require("tidyverse")

source("~/Documents/Academic/Active_projects/shared/Rscripts/base.R") # my special mix of personal functions
source("~/Documents/Academic/Active_projects/shared/Rscripts/base_graphics.R")

require(fields)


#-------------------------------------------------------------------------
# read and prep in data

sites <- read.csv("~/Documents/Academic/Active_projects/plant-associate-lit-search/data/Garrido2012/Garrido2012.txt", h=T, sep="\t")

# calculate GS distances between sites
Gdists <- rdist.earth(sites[,c("Lon","Lat")],sites[,c("Lon","Lat")], miles=FALSE)
colnames(Gdists) <- rownames(Gdists) <- sites$Site
Gdists[Gdists<1e-3] <- 0
Gdists <- as.data.frame(Gdists)%>% gather("To", "Gdist") %>% mutate(From = rep(rownames(Gdists), ncol(Gdists))) %>% filter(!is.na(Gdist)) %>% select(From, To, Gdist) -> Gdists


write.table(Gdists, "~/Documents/Academic/Active_projects/plant-associate-lit-search/data/Garrido2012/geo-dists.txt", sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)


