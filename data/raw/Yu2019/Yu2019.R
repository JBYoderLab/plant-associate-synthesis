# YU et al 2019 for meta-analysis
# reprocessing data for a synthesis of plant-associate datasets
# jby 2019.12.16

# starting up ------------------------------------------------------------
# setwd("~/Documents/Academic/Active_projects/plant-associate-lit-search/data/Yu2019")

require("tidyverse")

source("~/Documents/Academic/Active_projects/shared/Rscripts/base.R") # my special mix of personal functions
source("~/Documents/Academic/Active_projects/shared/Rscripts/base_graphics.R")

require("fields")


#-------------------------------------------------------------------------
# read and prep in data

sites = read.csv("~/Documents/Academic/Active_projects/plant-associate-lit-search/data/Yu2019/Yu2019_sites.txt", h=T, sep="\t")

plant.Fst = read.csv("~/Documents/Academic/Active_projects/plant-associate-lit-search/data/Yu2019/Fhirt_Gst.gdv", h=T, sep="\t", row.names="Obs.")
plant.Fst[lower.tri(plant.Fst, diag=TRUE)] <- NA
plant.Fst %>% gather("To", "Fst") %>% mutate(From = rep(rownames(plant.Fst), ncol(plant.Fst))) %>% filter(!is.na(Fst)) %>% select(From, To, Fst) -> plant.Fst
head(plant.Fst)
 
assoc.Fst = read.csv("~/Documents/Academic/Active_projects/plant-associate-lit-search/data/Yu2019/Valisia_Gst.gdv", h=T, sep="\t", row.names="Obs.")
assoc.Fst[lower.tri(assoc.Fst, diag=TRUE)] <- NA
assoc.Fst %>% gather("To", "Fst") %>% mutate(From = rep(rownames(assoc.Fst), ncol(assoc.Fst))) %>% filter(!is.na(Fst)) %>% select(From, To, Fst) -> assoc.Fst
head(assoc.Fst)


# calculate GS distances between sites
head(sites)

Gdists <- rdist.earth(as.matrix(sites[,c("Lon","Lat")]),as.matrix(sites[,c("Lon","Lat")]), miles=FALSE) #guhhh
colnames(Gdists) <- rownames(Gdists) <- sites$Site
Gdists[Gdists<1e-3] <- 0
Gdists[lower.tri(Gdists, diag=TRUE)] <- NA

Gdists <- as.data.frame(Gdists)%>% gather("To", "Gdist") %>% mutate(From = rep(rownames(Gdists), ncol(Gdists))) %>% filter(!is.na(Gdist)) %>% select(From, To, Gdist) -> Gdists

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

write.table(alldists, "~/Documents/Academic/Active_projects/plant-associate-lit-search/data/Yu2019/Yu2019-all-dists.txt", sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)














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

