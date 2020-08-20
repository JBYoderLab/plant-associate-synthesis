# Analysis of IBD in plants-associate data
# run locally
# jby 2020.08.18

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

dat = read.csv("data/Fst_working.txt", sep="\t", h=T) # cleaned, transformed data

str(dat)
unique(dat$Pair) # 17 because cleaned data

#-------------------------------------------------------------------------
# model fitting

descdist(dat$Assoc.Fst) # prolly not actually beta? (oh but maybe)
descdist(dat$afst.x) # suuuuper weird

# LOL this is a messssss
mlm.base <- lmer(Assoc.Fst~Plant.Fst+geo.x+(1|Pair), data=filter(dat, Pair %in% keepers)) # maybe?

summary(mlm.base) # not clear plant effect is nonzero; geography significant

# okay this is the basics ...
mod.base <- brm(afst.x~pfst.x+geo.x+(1|Pair), data=dat, family=hurdle_gamma) # yikes? this works??

save(mod.base, file="output/models/brms_mod.base.rdata")
load("output/models/brms_mod.base.rdata")

plot(mod.base) # diagnostics look good
summary(mod.base) # significantly nonzero plant effect, though much smaller than geography


# more hypothesis-testy ...
mod.types <- brm(afst.x~pfst.x+geo.x+Intxn.type+pfst.x:Intxn.type+(1|Pair), data=dat, family=hurdle_gamma) # yikes? this works??

save(mod.types, file="output/models/brms_mod.type.rdata")
load("output/models/brms_mod.type.rdata")

plot(mod.types) # diagnostics look good
summary(mod.types) # MUCH bigger plant effect; but changed by an interaction with interaction type --- 

# Population-Level Effects: 
#                            Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS
# Intercept                     -2.93      0.66    -4.26    -1.66 1.00     1588
# pfst.x                         3.54      0.41     2.73     4.32 1.00     3050
# geo.x                          0.30      0.04     0.22     0.38 1.00     4425
# Intxn.typeMutualism            0.70      0.77    -0.81     2.22 1.00     1018
# pfst.x:Intxn.typeMutualism    -3.50      0.41    -4.28    -2.68 1.00     3043

mod.types.red <- brm(afst.x~pfst.x+geo.x+Intxn.type+(1|Pair), data=dat, family=hurdle_gamma) # yikes? this works??

save(mod.types.red, file="output/models/brms_mod.type-red.rdata")
load("output/models/brms_mod.type-red.rdata")

plot(mod.types.red) # diagnostics look good
summary(mod.types.red) # MUCH bigger plant effect; but changed by an interaction with 




mod.assoc <- brm(afst.x~pfst.x+geo.x+Assoc.type+pfst.x:Assoc.type+(1|Pair), data=dat, family=hurdle_gamma) # yikes? this works??

save(mod.assoc, file="output/models/brms_mod.assoc.rdata")
load("output/models/brms_mod.assoc.rdata")

plot(mod.assoc) # diagnostics look good
summary(mod.assoc) # MUCH bigger plant effect; but changed by an interaction with 


mod.assoc.red <- brm(afst.x~pfst.x+geo.x+Assoc.type+(1|Pair), data=dat, family=hurdle_gamma) # yikes? this works??

save(mod.assoc.red, file="output/models/brms_mod.assoc.red.rdata")
load("output/models/brms_mod.assoc.red.rdata")

plot(mod.assoc.red) # diagnostics look good
summary(mod.assoc.red) # MUCH bigger plant effect; but changed by an interaction with 



LOO_all_mods <- LOO(mod.base, mod.types, mod.types.red, mod.assoc, mod.assoc.red, reloo=TRUE) # am I overfitting??
save(LOO_all_mods, file="output/models/brms_LOO_all.rdata")

LOO_all_mods

# Model comparisons:
#          elpd_diff se_diff
# mod.types   0.0       0.0 < !! 
# mod.all   -42.1      17.7





