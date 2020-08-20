# Analysis of IBD in plants-associate data
# raw correlations, bootstrapping
# run locally
# jby 2020.08.18

# starting up ------------------------------------------------------------
# setwd("~/Documents/Academic/Active_projects/plant-associate-synthesis")

require("tidyverse")
require("lme4")
require("brms")
require("fitdistrplus")
require("rsample")
require("meta")

source("~/Documents/Academic/Active_projects/shared/Rscripts/base.R") # my special mix of personal functions
source("~/Documents/Academic/Active_projects/shared/Rscripts/base_graphics.R")

#-------------------------------------------------------------------------
# read and prep in data

dat = read.csv("data/Fst_working.txt", h=T, sep="\t") # cleaned and transformed data

table(dat$Pair)
str(dat)

#-------------------------------------------------------------------------
# simple cor tests

cor.test(~pfst.x+geo.x, data=filter(dat, Pair == "Roridula gorgonius Pameridea roridulae"), method="spearman")

# who shows IBD? (sign. cor between Fst and geodist)
pibd <- dat %>% group_by(Pair) %>% do(broom::tidy(cor.test(~pfst.x+geo.x, data=., method="spearman"))) %>% mutate(test = "IBD, plant")
pibd %>% filter(p.value < 0.05) # 10/17 with IBD

pibd <- dat %>% group_by(Pair) %>% summarize(n=length(which(!is.na(pfst.x)))) %>% left_join(pibd) # I'll want ns later ...

aibd <- dat %>% group_by(Pair) %>% do(broom::tidy(cor.test(~afst.x+geo.x, data=., method="spearman"))) %>% mutate(test = "IBD, assoc.")
aibd %>% filter(p.value < 0.05) # 7/17 with IBD

aibd <- dat %>% group_by(Pair) %>% summarize(n=length(which(!is.na(afst.x)))) %>% left_join(aibd) # I'll want ns later ...

# who shows host-assoc cor?
avp <- dat %>% group_by(Pair) %>% do(broom::tidy(cor.test(~afst.x+pfst.x, data=., method="spearman"))) %>% mutate(test = "Plant v. assoc.")
avp %>% filter(p.value < 0.05) # 9/18 with host-assoc cor ...

avp <- dat %>% group_by(Pair) %>% summarize(n=length(which(!is.na(pfst.x)))) %>% left_join(avp) # I'll want ns later ...

avp %>% filter(estimate < 0) # LOL, Triponez ... and it's in the original paper! Weeeiird.


allcors <- rbind(pibd, aibd, avp)

allcors

write.table(allcors, "output/spearman-cors.txt", sep="\t", col.names=TRUE, row.names=FALSE) # save for later

allcors <- read.csv("output/spearman-cors.txt", sep="\t", h=TRUE) # later

#-------------------------------------------------------------------------
# bootstrapping instead 
nboot <- 1000 # tweak this to change N

# helper functions
pibd_split <- function(split) cor.test(~pfst.x+geo.x, analysis(split), method="spearman")
aibd_split <- function(split) cor.test(~afst.x+geo.x, analysis(split), method="spearman")
pva_split <- function(split) cor.test(~afst.x+pfst.x, analysis(split), method="spearman")

dat.cors.bs <- data.frame(matrix(0,0, 6)) # data frame looking for data woo
colnames(dat.cors.bs) <- c("Pair", "test", "id", "estimate", "statistic", "p.value")

for(pr in unique(dat$Pair)){

# pr <- keepers[1]

dat.pr.bs <- dat %>% filter(Pair==pr) %>% bootstraps(nboot)
pibd.bs <- dat.pr.bs %>% mutate(pibd = map(splits, pibd_split), pibd_info=map(pibd, tidy)) %>% unnest(pibd_info) %>% dplyr::select(id, estimate, statistic, p.value)
aibd.bs <- dat.pr.bs %>% mutate(aibd = map(splits, aibd_split), aibd_info=map(aibd, tidy)) %>% unnest(aibd_info) %>% dplyr::select(id, estimate, statistic, p.value)
pva.bs <- dat.pr.bs %>% mutate(pva = map(splits, pva_split), pva_info=map(pva, tidy)) %>% unnest(pva_info) %>% dplyr::select(id, estimate, statistic, p.value)

out <- data.frame(Pair=pr, test=rep(c("IBD, plant","IBD, assoc.","Plant v. assoc."), each=nboot), rbind(pibd.bs, aibd.bs, pva.bs))

dat.cors.bs <- rbind(dat.cors.bs, out)

}

# add annotations
dat.cors.bs <- dat %>% filter(!duplicated(Pair)) %>% dplyr::select(DOI, Plant.sp, Assoc.sp, Assoc.type, Intxn.type, Pair) %>% right_join(dat.cors.bs)

dim(dat.cors.bs)
head(dat.cors.bs)

# save for later
write.table(dat.cors.bs, "output/correlations_bootstrap.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

dat.cors.bs <- read.csv("output/correlations_bootstrap.txt", sep="\t", h=TRUE) # later



#-------------------------------------------------------------------------
# assembling a dataset for display

table(dat$Pair,dat$Intxn.type) # going to need to organize these

pairorda <- c("Datura stramonium/Lema trilineata", "Pinus banksiana/Arceuthobium americanum", "Pinus contorta/Arceuthobium americanum", "Populus angustifolia/Aceria parapopuli", "Rhus chinensis/Schlectendalia chinensis")
pairordm <- c("Ficus hirta/Valisia spp", "Ficus pumila/Wiebesia sp. 1", "Hirtella physophora/Allmoerus decemarticulatus", "Lysimachia vulgaris/Macropis europaea", "Lysimachia vulgaris/Macropis fulvipes", "Medicago lupulina/Ensifer medicae, E. meliloti", "Roridula gorgonius/Pameridea roridulae", "Silene latifolia/Hadena bicruris", "Vachellia drepanolobium/Crematogaster nigriceps", "Vachellia drepanolobium/Tetraponera penzigi", "Yucca brevifolia/Tegeticula synthetica", "Yucca jaegeriana/Tegeticula antithetica")

pairord <- c(pairorda, pairordm)

setdiff(pairord, dat$Pair)

cors.summ <- dat.cors.bs %>% mutate(Pair = factor(Pair, pairord), test = factor(test, c("IBD, plant","IBD, assoc.","Plant v. assoc."))) %>% group_by(Pair, test, Intxn.type, Assoc.type) %>% summarize(bs.med = median(estimate, na.rm=TRUE), lo.95 = quantile(estimate, 0.025, na.rm=TRUE), up.95 = quantile(estimate, 0.975, na.rm=TRUE), lo.50 = quantile(estimate, 0.25, na.rm=TRUE), up.50 = quantile(estimate, 0.75, na.rm=TRUE)) 

# create a vector of x-axis positions, for convenience (because I have no other good approach LOL)
cors.summ$plotpos <- rep(1:17,each=3)

cors.summ

# add in actual observed values!
cors.disp <- allcors %>% ungroup() %>% mutate(Pair = factor(Pair, pairord), test = factor(test, c("IBD, plant","IBD, assoc.","Plant v. assoc.")), Meta="Individual correlations") %>% dplyr::select(Meta, Pair, test, n, estimate) %>% rename(obs.cor=estimate) %>% right_join(cors.summ) %>% dplyr::select(-bs.med)

cors.disp

# nb, observed and bootstrap median are essentially identical

head(cors.disp)
dim(cors.disp)

table(cors.disp$Assoc.type) 

#-------------------------------------------------------------------------
# meta-analytic pooling within interaction types

head(cors.disp)

pibd.meta <- metacor(cor=cors.disp$obs.cor, n=cors.disp$n, subset=(cors.disp$test=="IBD, plant"))
aibd.meta <- metacor(cor=cors.disp$obs.cor, n=cors.disp$n, subset=(cors.disp$test=="IBD, assoc."))
pva.meta <- metacor(cor=cors.disp$obs.cor, n=cors.disp$n, subset=(cors.disp$test=="Plant v. assoc."))

# by Intxn.type
pibd.mut.meta <- metacor(cor=cors.disp$obs.cor, n=cors.disp$n, subset=(cors.disp$test=="IBD, plant"&cors.disp$Intxn.type=="Mutualism"))
aibd.mut.meta <- metacor(cor=cors.disp$obs.cor, n=cors.disp$n, subset=(cors.disp$test=="IBD, assoc."&cors.disp$Intxn.type=="Mutualism"))
pva.mut.meta <- metacor(cor=cors.disp$obs.cor, n=cors.disp$n, subset=(cors.disp$test=="Plant v. assoc."&cors.disp$Intxn.type=="Mutualism"))

pibd.ant.meta <- metacor(cor=cors.disp$obs.cor, n=cors.disp$n, subset=(cors.disp$test=="IBD, plant"&cors.disp$Intxn.type=="Antagonism"))
aibd.ant.meta <- metacor(cor=cors.disp$obs.cor, n=cors.disp$n, subset=(cors.disp$test=="IBD, assoc."&cors.disp$Intxn.type=="Antagonism"))
pva.ant.meta <- metacor(cor=cors.disp$obs.cor, n=cors.disp$n, subset=(cors.disp$test=="Plant v. assoc."&cors.disp$Intxn.type=="Antagonism"))

# by Assoc.type
pibd.ins.meta <- metacor(cor=cors.disp$obs.cor, n=cors.disp$n, subset=(cors.disp$test=="IBD, plant"&cors.disp$Assoc.type=="Insect"))
aibd.ins.meta <- metacor(cor=cors.disp$obs.cor, n=cors.disp$n, subset=(cors.disp$test=="IBD, assoc."&cors.disp$Assoc.type=="Insect"))
pva.ins.meta <- metacor(cor=cors.disp$obs.cor, n=cors.disp$n, subset=(cors.disp$test=="Plant v. assoc."&cors.disp$Assoc.type=="Insect"))

pibd.mic.meta <- metacor(cor=cors.disp$obs.cor, n=cors.disp$n, subset=(cors.disp$test=="IBD, plant"&cors.disp$Assoc.type=="Microbe"))
aibd.mic.meta <- metacor(cor=cors.disp$obs.cor, n=cors.disp$n, subset=(cors.disp$test=="IBD, assoc."&cors.disp$Assoc.type=="Microbe"))
pva.mic.meta <- metacor(cor=cors.disp$obs.cor, n=cors.disp$n, subset=(cors.disp$test=="Plant v. assoc."&cors.disp$Assoc.type=="Microbe"))


# faaaaaascinating
# for mutualism, strong meta-IBD in both players, weaker but significant pva
# for antagonism, weak meta-IBD in both players, strong pva

# extracting summary stats for comparison ....
# test for heterogeneity with p < 0.0001 in all cases but aibd.ant and pva.ant; going to use random effects methods across the board

# back-transform per https://www.meta-analysis.com/downloads/Meta-analysis%20Effect%20sizes%20based%20on%20correlations.pdf

rex <- function(maz){
return((exp(2*maz)-1)/(exp(2*maz)+1))
}

names(pibd.mut.meta)
pibd.mut.meta$TE.random
rex(pibd.mut.meta$TE.random) # pooled correlation
rex(pibd.mut.meta$lower.random) # lower CI
rex(pibd.mut.meta$upper.random) # upper CI

pooled <- data.frame(Meta="Pooled", Pair=rep(c("Mutualisms", "Antagonisms", "Insects", "Microbes"), each=3), n=c(12,12,12,5,5,5,14,14,14,3,3,3), plotpos=rep(2:1, each=3), test=factor(rep(c("IBD, plant", "IBD, assoc.", "Plant v. assoc."),4), c("IBD, plant","IBD, assoc.","Plant v. assoc.")), Intxn.type=c(rep(c("Mutualism", "Antagonism"), each=3), rep(NA, 6)), Assoc.type=c(rep("Both", 6), rep(c("Insect", "Microbe"), each=3)), obs.cor=c(rex(pibd.mut.meta$TE.random), rex(aibd.mut.meta$TE.random), rex(pva.mut.meta$TE.random), rex(pibd.ant.meta$TE.random), rex(aibd.ant.meta$TE.random), rex(pva.ant.meta$TE.random), rex(pibd.ins.meta$TE.random), rex(aibd.ins.meta$TE.random), rex(pva.ins.meta$TE.random), rex(pibd.mic.meta$TE.random), rex(aibd.mic.meta$TE.random), rex(pva.mic.meta$TE.random)), lo.95=c(rex(pibd.mut.meta$lower.random), rex(aibd.mut.meta$lower.random), rex(pva.mut.meta$lower.random), rex(pibd.ant.meta$lower.random), rex(aibd.ant.meta$lower.random), rex(pva.ant.meta$lower.random), rex(pibd.ins.meta$lower.random), rex(aibd.ins.meta$lower.random), rex(pva.ins.meta$lower.random), rex(pibd.mic.meta$lower.random), rex(aibd.mic.meta$lower.random), rex(pva.mic.meta$lower.random)), up.95=c(rex(pibd.mut.meta$upper.random), rex(aibd.mut.meta$upper.random), rex(pva.mut.meta$upper.random), rex(pibd.ant.meta$upper.random), rex(aibd.ant.meta$upper.random), rex(pva.ant.meta$upper.random), rex(pibd.ins.meta$upper.random), rex(aibd.ins.meta$upper.random), rex(pva.ins.meta$upper.random), rex(pibd.mic.meta$upper.random), rex(aibd.mic.meta$upper.random), rex(pva.mic.meta$upper.random)), lo.50=NA, up.50=NA)

pooled # finally

ggplot() + geom_hline(yintercept=0, lty=2) + geom_linerange(data=pooled, aes(x=plotpos, ymax=up.95, ymin=lo.95, color=Intxn.type), size=0.75) + geom_point(data=pooled, aes(x=plotpos, y=obs.cor, shape=test)) + scale_color_manual(values=colors[c(2,5)], name="Type") + scale_shape_manual(values=c(21,25,23), name="Correlation", labels=c("IBD, plant", "IBD, assoc.", "Plant v. assoc.")) + scale_x_continuous(breaks=c(6,2), labels=c("Mutualism", "Antagonism")) + coord_flip() + labs(x = "Plant/Associate pair", y = "Spearman rank correlation") + theme_ipsum(axis_title_size=13) + theme(legend.position="bottom", legend.title=element_blank())


global <- data.frame(Meta="", Pair="Global", n=17, plotpos=1.5, test=factor(c("IBD, plant", "IBD, assoc.", "Plant v. assoc."), c("IBD, plant","IBD, assoc.","Plant v. assoc.")), Intxn.type=NA, Assoc.type="Both", obs.cor=c(rex(pibd.meta$TE.random), rex(aibd.meta$TE.random), rex(pva.meta$TE.random)), lo.95=c(rex(pibd.meta$lower.random), rex(aibd.meta$lower.random), rex(pva.meta$lower.random)), up.95=c(rex(pibd.meta$upper.random), rex(aibd.meta$upper.random), rex(pva.meta$upper.random)), lo.50=NA, up.50=NA)


cors.meta <- rbind(cors.disp,pooled,global) %>% mutate(Meta = factor(Meta, c("Individual correlations", "Pooled", "")))

write.table(cors.meta, "output/cors_figure_data.txt", sep="\t", col.names=TRUE, row.names=FALSE) # save for later

cors.meta <- read.csv("output/cors_figure_data.txt", sep="\t", h=TRUE) %>% mutate(test=factor(test, c("IBD, plant", "IBD, assoc.", "Plant v. assoc.")), Pair=factor(Pair, c("Datura stramonium/Lema trilineata", "Pinus banksiana/Arceuthobium americanum", "Pinus contorta/Arceuthobium americanum", "Populus angustifolia/Aceria parapopuli", "Rhus chinensis/Schlectendalia chinensis", "Ficus hirta/Valisia spp", "Ficus pumila/Wiebesia sp. 1", "Hirtella physophora/Allmoerus decemarticulatus", "Lysimachia vulgaris/Macropis europaea", "Lysimachia vulgaris/Macropis fulvipes", "Medicago lupulina/Ensifer medicae, E. meliloti", "Roridula gorgonius/Pameridea roridulae", "Silene latifolia/Hadena bicruris", "Vachellia drepanolobium/Crematogaster nigriceps", "Vachellia drepanolobium/Tetraponera penzigi", "Yucca brevifolia/Tegeticula synthetica", "Yucca jaegeriana/Tegeticula antithetica", "Mutualism, pooled", "Antagonism, pooled", "Insects, pooled", "Microbes, pooled", "Global"))) # later



#-------------------------------------------------------------------------
# building figures

colors <- park_palette("BlueRidgePkwy")

# here's individual correlations PLUS meta-analysis pooled and global cors

{cairo_pdf("output/figures/spearman_cors_meta.pdf", height=7, width=8.5)

ggplot() + geom_hline(yintercept=0, lty=2) + 
geom_linerange(data=cors.meta, aes(x=Pair, ymax=up.95, ymin=lo.95, color=Intxn.type), lwd=0.5) + 
geom_linerange(data=cors.meta, aes(x=Pair, ymax=up.50, ymin=lo.50, color=Intxn.type), lwd=1.25) + 
geom_point(data=cors.meta, aes(x=Pair, y=obs.cor, shape=Assoc.type, fill=Intxn.type, size=Meta)) + 

scale_color_manual(values=colors[c(2,5)], na.translate=TRUE, na.value=colors[4], name="Type") + 
scale_fill_manual(values=colors[c(2,5)], na.translate=TRUE, na.value=colors[4], name="Type") + 
scale_size_manual(values=c(1.5,2.5,2.5), guide=NULL) + 
scale_shape_manual(values=c(21,25,23), name="Correlation", labels=c("IBD, plant", "IBD, assoc.", "Plant v. assoc.")) + 

coord_flip() + facet_grid(Meta~test, scale="free_y", space="free_y") + labs(x = "Plant/Associate pair", y = "Spearman rank correlation") + theme_ipsum(axis_title_size=13) + theme(legend.position="none", legend.title=element_blank(), panel.spacing.y=unit(0.5,"line"))

}
dev.off()

# just the individual study correlations
indcors <- filter(cors.meta, Meta=="Individual correlations")

{cairo_pdf("output/figures/spearman_cors.pdf", height=6, width=8.5)

ggplot() + geom_hline(yintercept=0, color="gray65") + 
geom_linerange(data=indcors, aes(x=Pair, ymax=up.95, ymin=lo.95, color=Intxn.type), lwd=0.75) + 
geom_linerange(data=indcors, aes(x=Pair, ymax=up.50, ymin=lo.50, color=Intxn.type), lwd=1.5) + 
geom_point(data=indcors, aes(x=Pair, y=obs.cor, shape=Assoc.type), color="white", size=4) +
geom_point(data=indcors, aes(x=Pair, y=obs.cor, shape=Assoc.type, color=Intxn.type), size=3) + 

scale_color_manual(values=colors[c(2,5)]) + 
scale_fill_manual(values=colors[c(2,5)]) +  
scale_shape_manual(values=c(18,20), name="Correlation", labels=c("Insect", "Microbe")) +  

coord_flip() + facet_wrap("test") + labs(x = "Plant/Associate pair", y = "Spearman rank correlation (± 50 and 95% CI)") + theme_ipsum(axis_title_size=13) + theme(legend.position="bottom", legend.title=element_blank(), panel.spacing.x=unit(1,"line"))

}
dev.off()


# just the pooled correlations

globcor <- dplyr::filter(cors.meta, Meta=="") %>% dplyr::select(-Meta)
pooled.disp <- dplyr::filter(cors.meta, Meta%in%c("Pooled", ""))
pooled.disp$set <- c(rep("Intxn", 6), rep("Assoc", 6), rep("Global", 3))

{cairo_pdf("output/figures/spearman_cors_meta_pooled.pdf", height=4, width=8)

ggplot() + geom_hline(yintercept=0, color="gray65") + geom_hline(data=globcor, aes(yintercept=obs.cor), lty=2, lwd=0.5, color=colors[6]) +
geom_linerange(data=pooled.disp, aes(x=Pair, ymax=up.95, ymin=lo.95, color=set), lwd=0.75) + 
geom_point(data=pooled.disp, aes(x=Pair, y=obs.cor, fill=set, size=n), pch=21, color="white") + 

scale_color_manual(values=colors[c(1,6,4)]) + 
scale_fill_manual(values=colors[c(1,6,4)]) + 

coord_flip() + facet_wrap("test") + labs(x = NULL, y = "Spearman rank correlation (± 95% CI)") + theme_ipsum(axis_title_size=13) + theme(legend.position="none", legend.title=element_blank(), panel.spacing.y=unit(0.5,"line"))

}
dev.off() 

