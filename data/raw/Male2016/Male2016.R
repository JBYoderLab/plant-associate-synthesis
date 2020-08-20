# MALE et al 2016 for meta-analysis
# reprocessing data for a synthesis of plant-associate datasets
# jby 2019.12.17

# starting up ------------------------------------------------------------
# setwd("~/Documents/Academic/Active_projects/plant-associate-lit-search/data/Male2016")

require("tidyverse")

source("~/Documents/Academic/Active_projects/shared/Rscripts/base.R") # my special mix of personal functions
source("~/Documents/Academic/Active_projects/shared/Rscripts/base_graphics.R")

require("fields")


#-------------------------------------------------------------------------
# read and prep in data


plant.Fst = read.csv("~/Documents/Academic/Active_projects/plant-associate-lit-search/data/Male2016/Hirtella_Gst2.gdv", h=T, sep="\t", row.names="Obs.")
plant.Fst[lower.tri(plant.Fst, diag=TRUE)] <- NA
plant.Fst %>% gather("To", "Fst") %>% mutate(From = rep(rownames(plant.Fst), ncol(plant.Fst))) %>% filter(!is.na(Fst)) %>% select(From, To, Fst) -> plant.Fst
head(plant.Fst)
 
assoc.Fst = read.csv("~/Documents/Academic/Active_projects/plant-associate-lit-search/data/Male2016/Allomerus_Gst2.gdv", h=T, sep="\t", row.names="Obs.")
assoc.Fst[lower.tri(assoc.Fst, diag=TRUE)] <- NA
assoc.Fst %>% gather("To", "Fst") %>% mutate(From = rep(rownames(assoc.Fst), ncol(assoc.Fst))) %>% filter(!is.na(Fst)) %>% select(From, To, Fst) -> assoc.Fst
head(assoc.Fst)


# calculate GS distances between sites
head(sites)

Gdists = read.csv("~/Documents/Academic/Active_projects/plant-associate-lit-search/data/Male2016/Male2016_geo.txt", h=T, sep="\t")

Gdists


# bit by bit, putting it together
dim(plant.Fst)
head(plant.Fst)
dim(assoc.Fst)
head(assoc.Fst)
dim(Gdists)
head(Gdists)

intersect(plant.Fst$From, assoc.Fst$From)

alldists <- plant.Fst %>% inner_join(assoc.Fst, by=c("From","To"), suffix=c(".plant",".assoc")) %>% inner_join(Gdists, by=c("From","To"))

dim(alldists)
head(alldists)

write.table(alldists, "~/Documents/Academic/Active_projects/plant-associate-lit-search/data/Male2016/Male2016-all-dists.txt", sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)


