# EVANS et al 2013 for meta-analysis
# reprocessing data for a synthesis of plant-associate datasets
# jby 2020.02.07

# starting up ------------------------------------------------------------
# setwd("~/Documents/Academic/Active_projects/plant-associate-lit-search/data/Evans2013")

require("tidyverse")

source("~/Documents/Academic/Active_projects/shared/Rscripts/base.R") # my special mix of personal functions
source("~/Documents/Academic/Active_projects/shared/Rscripts/base_graphics.R")

require("fields")


#-------------------------------------------------------------------------
# read and prep in data

sites = read.csv("~/Documents/Academic/Active_projects/plant-associate-lit-search/data/Evans2013/Evans2013_sites.txt", sep="\t", h=T)


# calculate GS distances between sites
head(sites)

Gdists <- rdist.earth(as.matrix(sites[,c("Lon","Lat")]),as.matrix(sites[,c("Lon","Lat")]), miles=FALSE) #guhhh
colnames(Gdists) <- rownames(Gdists) <- sites$Site
Gdists[Gdists<1e-3] <- 0
Gdists[lower.tri(Gdists, diag=TRUE)] <- NA

Gdists <- as.data.frame(Gdists) %>% gather("To", "Gdist") %>% mutate(From = rep(rownames(Gdists), ncol(Gdists))) %>% filter(!is.na(Gdist)) %>% dplyr::select(From, To, Gdist)

Gdists


write.table(Gdists, "~/Documents/Academic/Active_projects/plant-associate-lit-search/data/Evans2013/Evans2013-geo-dists.txt", sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)




