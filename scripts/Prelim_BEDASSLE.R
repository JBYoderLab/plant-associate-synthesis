# Attempting BEDASSLE with compiled plants-associate data
# run locally
# jby 2020.02.21

# starting up ------------------------------------------------------------
# setwd("~/Documents/Academic/Active_projects/plant-associate-lit-search")

require("tidyverse")
require("lme4")
require("brms")
require("fitdistrplus")
require("rsample")
require("bedassle")
require("dils")

source("~/Documents/Academic/Active_projects/shared/Rscripts/base.R") # my special mix of personal functions
source("~/Documents/Academic/Active_projects/shared/Rscripts/base_graphics.R")

#-------------------------------------------------------------------------
# read and prep in data

dat = read.csv("data/plant-associate-lit-search_Fst.csv", h=T)

head(dat)
dim(dat)

dat$Pair = paste(dat$Author, dat$Plant.sp)
dat <- filter(dat, Geog.dist.km > 0)

#-------------------------------------------------------------------------
# data munging for a test run

table(dat$Pair)

# function to create bedassle input matrices from compiled data format
long2matrix <- function(dat, pr){

	# pr <- "Yoder Yucca brevifolia" 

	seg <- filter(dat, Pair == pr) %>% dplyr::select(Pop1, Pop2, Plant.Fst, Assoc.Fst, Geog.dist.km)	
	sites <- sort(unique(c(seg$Pop1, seg$Pop2)))
	
	adist <- seg %>% dplyr::select(Pop1, Pop2, Assoc.Fst) %>% pivot_wider(names_from = Pop2, values_from = Assoc.Fst)
	adist.mat <- as.matrix(adist[,-1])
	rownames(adist.mat) <- adist$Pop1
	adist.mat <- adist.mat[sites, sites]
	adist.mat[is.na(adist.mat)] <- 0
	
	pdist <- seg.c %>% dplyr::select(Pop1, Pop2, Plant.Fst) %>% pivot_wider(names_from = Pop2, values_from = Plant.Fst)
	pdist.mat <- as.matrix(pdist[,-1])
	rownames(pdist.mat) <- pdist$Pop1
	pdist.mat <- pdist.mat[sites, sites]
	pdist.mat[is.na(pdist.mat)] <- 0
	
	gdist <- seg.c %>% dplyr::select(Pop1, Pop2, Geog.dist.km) %>% pivot_wider(names_from = Pop2, values_from = Geog.dist.km)
	gdist.mat <- as.matrix(gdist[,-1])
	rownames(gdist.mat) <- gdist$Pop1
	gdist.mat <- gdist.mat[sites, sites]
	gdist.mat[is.na(gdist.mat)] <- 0
	
	return(list(assoc=adist.mat, plant=pdist.mat, geo=gdist.mat))

}

test.data <- long2matrix(dat, "Yoder Yucca brevifolia") # trial run
test.data

#-------------------------------------------------------------------------
# run BEDASSLE

test.bd <- run.bedassle(
	genDist = test.data$assoc,
	geoDist = test.data$geo,
	envDist = test.data$plant,
	nLoci = 8,
	prefix = "Magalhaes",
	n.chains = 4,
	n.iter = 2000,
	make.figs = TRUE,
	save.stanfit = TRUE
) # hmmm

pwp2allelicCov(test.data$assoc)















summary(dat)
str(dat) # what have we got here?
table(dat$DOI)
table(dat$Assoc.type) # LOTS more mutualists huh
table(dat$DOI,dat$Assoc.type)
length(unique(dat$DOI)) # hey, 15 studies
length(unique(paste(dat$DOI, dat$Plant.sp, dat$Assoc.sp))) # hey, 19 study-plant-assoc combos
unique(dat$Author)

hist(dat$Plant.Fst)
hist(dat$Assoc.Fst)
hist(dat$Geog.dist.km)

# some data-munging --------------

dat$Plant.Fst[dat$Plant.Fst<0] <- 0
dat$Assoc.Fst[dat$Assoc.Fst<0] <- 0

dat$Pair <- paste(dat$Plant.sp, dat$Assoc.sp)

dat <- dat %>% mutate(pfst.x = Plant.Fst/(1-Plant.Fst), afst.x = Assoc.Fst/(1-Assoc.Fst), geo.x = log10(Geog.dist.km)) %>% group_by(DOI, Plant.sp) %>% mutate(pfst.x.st = (pfst.x-mean(pfst.x, na.rm=T))/sd(pfst.x, na.rm=T), afst.x.st = (afst.x-mean(afst.x, na.rm=T))/sd(afst.x, na.rm=T), geo.x.st = (geo.x-mean(geo.x, na.rm=T))/sd(geo.x, na.rm=T))

#---------------------------------

hist(dat$pfst.x) # transformed (whew)
hist(dat$pfst.x.st) # standardized within species/studies
hist(dat$afst.x) # transformed (whew)
hist(dat$afst.x.st) # standardized within species/studies
hist(dat$geo.x) # transformed (whew)
hist(dat$geo.x.st) # standardized within species/studies

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
keepers # 16 retained okay

#-------------------------------------------------------------------------
# Statistical stuff

pibd <- dat %>% group_by(Pair, Intxn.type) %>% filter(length(pfst.x)>2) %>% do(broom::tidy(cor.test(~pfst.x+geo.x, data=., method="spearman")))
pibd %>% filter(p.value < 0.05) # only 9/18 (one has too few observations)

aibd <- dat %>% group_by(Pair, Intxn.type) %>% filter(length(pfst.x)>2) %>% do(broom::tidy(cor.test(~afst.x+geo.x, data=., method="spearman")))
aibd %>% filter(p.value < 0.05) # only 7/18

avp <- dat %>% group_by(Pair, Intxn.type) %>% filter(length(pfst.x)>2) %>% do(broom::tidy(cor.test(~afst.x+pfst.x, data=., method="spearman")))
avp %>% filter(p.value < 0.05) # only 7/18
avp %>% filter(estimate < 0) # it's Triponez ... and it's in the original paper! Weeeiird.

# bootstrapping instead ---------------
nboot <- 1000 # tweak this to change N

# helper functions
pibd_split <- function(split) cor.test(~pfst.x+geo.x, analysis(split), method="spearman")
aibd_split <- function(split) cor.test(~afst.x+geo.x, analysis(split), method="spearman")
pva_split <- function(split) cor.test(~afst.x+pfst.x, analysis(split), method="spearman")

dat.cors.bs <- data.frame(matrix(0,0, 6)) # data frame looking for data woo
colnames(dat.cors.bs) <- c("Pair", "test", "id", "estimate", "statistic", "p.value")

for(pr in keepers){

# pr <- keepers[1]

dat.pr.bs <- dat %>% filter(Pair==pr) %>% bootstraps(nboot)
pibd.bs <- dat.pr.bs %>% mutate(pibd = map(splits, pibd_split), pibd_info=map(pibd, tidy)) %>% unnest(pibd_info) %>% dplyr::select(id, estimate, statistic, p.value)
aibd.bs <- dat.pr.bs %>% mutate(aibd = map(splits, aibd_split), aibd_info=map(aibd, tidy)) %>% unnest(aibd_info) %>% dplyr::select(id, estimate, statistic, p.value)
pva.bs <- dat.pr.bs %>% mutate(pva = map(splits, pva_split), pva_info=map(pva, tidy)) %>% unnest(pva_info) %>% dplyr::select(id, estimate, statistic, p.value)

out <- data.frame(Pair=pr, test=rep(c("pibd", "aibd", "pva"), each=nboot), rbind(pibd.bs, aibd.bs, pva.bs))

dat.cors.bs <- rbind(dat.cors.bs, out)

}

# add annotations
dat.cors.bs <- dat %>% filter(!duplicated(Pair)) %>% dplyr::select(Author, Year, DOI, Plant.sp, Assoc.sp, Assoc.type, Intxn.type, Pair) %>% right_join(dat.cors.bs)

dim(dat.cors.bs)
head(dat.cors.bs)

# save for later
write.table(dat.cors.bs, "output/correlations_bootstrap.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")


#--------------------------------------

descdist(dat$Assoc.Fst) # prolly not actually beta? (oh but maybe)
descdist(dat$afst.x) # suuuuper weird
descdist(dat$afst.x.st) # lognormal, nice

# LOL this is a messssss
mlm.all <- lmer(Assoc.Fst~Plant.Fst+geo.x+Intxn.type+(1|Pair), data=dat) # maybe? nope

summary(mlm.all)

anova(mlm.all)


# get me residuals of Assoc.Fst~Geog.dist.km

ag <- lm(Assoc.Fst~Geog.dist.km, data=dat)

ag$residuals

dat$Assoc.resid[!is.na(dat$Geog.dist.km)] <- ag$residuals

#-------------------------------------------------------------------------
# Visualization

ggplot(dat, aes(x=Plant.Fst, y=Assoc.Fst, color=interaction(Author, Plant.sp))) + geom_point(alpha=0.3) + geom_smooth(method="lm", se=FALSE)


# per-study correlations, visualized
ggplot() + geom_smooth(data=filter(dat, Pair %in% filter(pibd, p.value < 0.05)$Pair), aes(x=geo.x, y=pfst.x), method="lm", color="white") + geom_point(data=dat, aes(x=geo.x, y=pfst.x, color=Intxn.type), alpha=0.5) + scale_color_manual(values=c("darkorange", "forestgreen")) + facet_wrap("Pair", scale="free") + labs(y=expression(F[ST]/(1-F[ST])), x=expression(log[10]("geographic distance")), title="Isolation by distance, plants") + theme_ipsum(axis_title_just="lt") + theme(legend.position=c(0.9, 0.1))

ggplot() + geom_smooth(data=filter(dat, Pair %in% filter(aibd, p.value < 0.05)$Pair), aes(x=geo.x, y=afst.x), method="lm", color="white") + geom_point(data=dat, aes(x=geo.x, y=afst.x, color=Intxn.type), alpha=0.5) + scale_color_manual(values=c("darkorange", "forestgreen")) + facet_wrap("Pair", scale="free") + labs(y=expression(F[ST]/(1-F[ST])), x=expression(log[10]("geographic distance")), title="Isolation by distance, associates") + theme_ipsum(axis_title_just="lt") + theme(legend.position=c(0.9, 0.1))

ggplot() + geom_smooth(data=filter(dat, Pair %in% filter(avp, p.value < 0.05)$Pair), aes(x=pfst.x, y=afst.x), method="lm", color="white") + geom_point(data=dat, aes(x=pfst.x, y=afst.x, color=Intxn.type), alpha=0.5) + scale_color_manual(values=c("darkorange", "forestgreen")) + facet_wrap("Pair", scale="free") + labs(y=expression(F[ST]/(1-F[ST])), x=expression(F[ST]/(1-F[ST])), title="Plant vs associate genetic distance") + theme_ipsum(axis_title_just="lt") + theme(legend.position=c(0.9, 0.1))



# IBD and PvA correlations, bootstrapped
cors.sum <- dat.cors.bs %>% group_by(Pair, test, Assoc.type, Intxn.type) %>% summarize(mncor = mean(estimate, na.rm=TRUE), mdcor=median(estimate, na.rm=TRUE), lo=quantile(estimate, 0.025, na.rm=TRUE), hi=quantile(estimate, 0.975, na.rm=TRUE))

head(cors.sum)

ggplot(cors.sum, aes(x=interaction(test, Pair), y=mdcor, ymin=lo, ymax=hi, color=Intxn.type, shape=test)) + geom_abline(slope=0, intercept=0, color="gray70") + geom_pointrange() + coord_flip() + theme(legend.position="none") # hmmm




