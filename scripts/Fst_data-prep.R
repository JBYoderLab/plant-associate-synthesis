# Analysis of IBD in plants-associate data
# run locally
# jby 2020.08.10

# starting up ------------------------------------------------------------
# setwd("~/Documents/Academic/Active_projects/plant-associate-synthesis")

require("tidyverse")
require("lme4")
require("brms")
require("fitdistrplus")
require("rsample")

source("~/Documents/Academic/Active_projects/shared/Rscripts/base.R") # my special mix of personal functions
source("~/Documents/Academic/Active_projects/shared/Rscripts/base_graphics.R")

#-------------------------------------------------------------------------
# read and prep in data

dat = read.csv("data/plant-associate-lit-search_Fst.csv", h=T)

head(dat)
dim(dat)

summary(dat)
str(dat) # what have we got here?
table(dat$Author)
table(dat$Assoc.type)
table(dat$DOI,dat$Intxn.type) # LOTS more mutualists huh
table(dat$Assoc.type, dat$Intxn.type)
table(dat$DOI,dat$Assoc.type)

length(unique(dat$DOI)) # hey, 15 studies
length(unique(paste(dat$DOI, dat$Plant.sp, dat$Assoc.sp))) # to do: break up Joshua trees
unique(dat$Author)

hist(dat$Plant.Fst)
hist(dat$Assoc.Fst)
hist(dat$Geog.dist.km)

# some data-munging --------------

dat$Plant.Fst[dat$Plant.Fst<0] <- 0
dat$Assoc.Fst[dat$Assoc.Fst<0] <- 0
dat <- filter(dat, Geog.dist.km > 0)

dat$Pair <- paste(dat$Plant.sp, dat$Assoc.sp, sep="/")

hist(dat$Plant.Fst/(1-dat$Plant.Fst))
hist(dat$Assoc.Fst/(1-dat$Assoc.Fst))
hist(log10(dat$Geog.dist.km))

dat <- dat %>% mutate(pfst.x = Plant.Fst/(1-Plant.Fst), afst.x = Assoc.Fst/(1-Assoc.Fst), geo.x = log10(Geog.dist.km))

#---------------------------------

hist(dat$pfst.x) # transformed (whew)
hist(dat$afst.x) # transformed (whew)
hist(dat$geo.x) # transformed (whew)

length(which(is.na(dat$pfst.x)))
length(which(is.nan(dat$pfst.x)))
length(which(dat$pfst.x==Inf))

length(which(is.na(dat$afst.x)))
length(which(is.nan(dat$afst.x)))
length(which(dat$afst.x==Inf))

length(which(is.na(dat$geo.x)))
length(which(is.nan(dat$geo.x)))
length(which(dat$geo.x==Inf))

dat <- filter(dat, afst.x!=Inf, pfst.x!=Inf)

table(dat$Pair)
keepers <- names(which(table(dat$Pair)>5)) 
keepers # 17 retained okay

dat.working <- dat %>% filter(Pair %in% keepers) %>% dplyr::select(DOI, Plant.sp, Assoc.sp, Pair, Intxn.type, Assoc.type, Pop1, Pop2, Plant.Fst, Assoc.Fst, Geog.dist.km, pfst.x, afst.x, geo.x, Entry.by)

write.table(dat.working, "data/Fst_working.txt", sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE) # and this is what we're using going forward





