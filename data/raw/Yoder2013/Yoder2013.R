# YODER et al 2013 for meta-analysis
# reprocessing data for a synthesis of plant-associate datasets
# jby 2019.10.29

# starting up ------------------------------------------------------------
# setwd("~/Documents/Academic/Active_projects/plant-associate-lit-search/data/Yoder2013")

require("tidyverse")

source("~/Documents/Academic/Active_projects/shared/Rscripts/base.R") # my special mix of personal functions
source("~/Documents/Academic/Active_projects/shared/Rscripts/base_graphics.R")

require(fields)


#-------------------------------------------------------------------------
# read and prep in data

sites = read.csv("/Users/jyoder/Documents/Academic/Active_projects/plant-associate-lit-search/data/Yoder2013/Yoder_etal_sites.txt", h=T, sep="\t")

ybr.Fst = read.csv("/Users/jyoder/Documents/Academic/Active_projects/plant-associate-lit-search/data/Yoder2013/Ybrev_Gst.gdv", h=T, sep="\t", row.names="Obs.")
ybr.Fst[lower.tri(ybr.Fst)] <- NA
ybr.Fst %>% gather("To", "Fst") %>% mutate(From = rep(rownames(ybr.Fst), ncol(ybr.Fst))) %>% filter(!is.na(Fst)) %>% select(From, To, Fst) -> ybr.Fst
head(ybr.Fst)
 
tsy.Fst = read.csv("/Users/jyoder/Documents/Academic/Active_projects/plant-associate-lit-search/data/Yoder2013/Tsyn_Gst.gdv", h=T, sep="\t", row.names="Obs.")
tsy.Fst[lower.tri(tsy.Fst)] <- NA
tsy.Fst %>% gather("To", "Fst") %>% mutate(From = rep(rownames(tsy.Fst), ncol(tsy.Fst))) %>% filter(!is.na(Fst)) %>% select(From, To, Fst) -> tsy.Fst
head(tsy.Fst)

tan.Fst = read.csv("/Users/jyoder/Documents/Academic/Active_projects/plant-associate-lit-search/data/Yoder2013/Tant_Gst.gdv", h=T, sep="\t", row.names="Obs.")
tan.Fst[lower.tri(tan.Fst)] <- NA
tan.Fst %>% gather("To", "Fst") %>% mutate(From = rep(rownames(tan.Fst), ncol(tan.Fst))) %>% filter(!is.na(Fst)) %>% select(From, To, Fst) -> tan.Fst
head(tan.Fst)

head(sites)

# calculate GS distances between sites
gsdists <- data.frame(matrix(0,0,3))
names(gsdists) <- c("From", "To", "dist")
gsdists

for(fr in 1:nrow(sites)){

for(to in fr:nrow(sites)){

if(fr == to) out <- data.frame(From=sites[fr,"Population"], To=sites[to,"Population"], dist=0)


if(fr != to) out <- data.frame(From=sites[fr,"Population"], To=sites[to,"Population"], dist=rdist.earth(sites[fr,c("lon_dec","lat_dec")],sites[to,c("lon_dec","lat_dec")], miles=FALSE))

gsdists <- rbind(gsdists, out)

}

}

head(gsdists)
hist(gsdists$dist)

dim(gsdists)
dim(ybr.Fst)
dim(tsy.Fst)

