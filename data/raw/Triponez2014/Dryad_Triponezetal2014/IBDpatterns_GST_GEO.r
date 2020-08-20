####### Mod 2019.12.11 jby

require(ade4)
require(maps)
require(geosphere)
require(vegan)

rm(list=ls())

setwd("~/Documents/Academic/Active_projects/plant-associate-lit-search/data/Triponez2014/Dryad_Triponezetal2014")

### Functions
NicePlot = function(X, Y, name = "Scatterplot_IBD_MF.pdf"){
  pdf(file = name)
  tline = data.frame(geo = X, insect = Y)
  tline = tline[ rowSums(tline) > 0,]
  plot(tline, xlab = "Dist. Geo between populations (Km)", ylab = "GST", bty = 'n', col = "grey", pch = 16, cex = .5)
  test = lm(insect~geo, data = tline)
  abline(test)
  test = summary(test)
  coeff = test$coefficients
  Rsq = test$r.squared
  text(min(tline[,1], na.rm = T) + 0.75 * (max(tline[,1], na.rm = T) - min(tline[,1], na.rm = T)), min(tline[,2], na.rm = T) + 0.05 * (max(tline[,2], na.rm = T) - min(tline[,2], na.rm = T)), labels = paste("Effect = ", signif(coeff[2, 1], 2), "\nRsq =", round(Rsq, 2), "\nPval = ", signif(coeff[2, 4], 2)), cex = 1.2)
  dev.off()
  }

COSINE = function(X, Y){
  geo = cbind(X, Y)
  ix = 1:length(X)
  geo.d = outer(ix[-1], ix[-length(ix)], function(ivec, jvec, ...) sapply(seq(along = ivec),
    function(k) {
      i = ivec[k]
      j = jvec[k]
      if (i > j)
      distCosine(geo[i,1:2], geo[j, 1:2]) / 1000
      else NA
      }))
  as.matrix(as.dist(rbind(NA, cbind(geo.d, NA))))
  }

GST_pairwise = function(matm, pops){
  # prepare data
  pops = as.character(pops)
  freqsP = aggregate(matm, by=list(pops), mean)[, -1] #get band frequencies per pops
  logs = freqsP * log2(freqsP)
  logs[is.na(logs)] = 0
  poplevs = levels(as.factor(pops))
  level.names = poplevs

  # compute pairwise GSTs
  GST_fast=function(i, j, ...){
	  matm1 = matm[pops==i | pops==j, ]
	  freqsT = colMeans(matm1)
	  logsT = freqsT * log2(freqsT)
	  logsT[is.na(logsT)] = 0
	  shanT = -mean(logsT)
	  logsP1 = logs[poplevs==i | poplevs==j, ]
	  shanP = -mean(rowMeans(logsP1))
	  (shanT - shanP) / shanT
	  }

  ix = seq(along = level.names)
  names(ix) = level.names
  pp = outer(ix[-1], ix[-length(ix)], function(ivec, jvec, ...) sapply(seq(along = ivec),
    function(k) {
      i = ivec[k]
      j = jvec[k]
      if (i > j)
      GST_fast(names(ix)[i], names(ix)[j], matm, logs, poplevs)
      else NA
      }))
  pp2 = rbind(pp[, 1], pp)
  pp3 = as.matrix(as.dist(cbind(pp2, NA)))
  pp3[is.na(pp3)] = 0
  rownames(pp3) = colnames(pp3) = level.names
  pp3
  }


################
### Prepare Lysi
plant = read.delim("Lysi439ind497loc.txt", header = F, row.names = 1)
plant = plant[ -1, -1]
plant = plant[, -ncol(plant)]
pops = substr(rownames(plant), 0, 3)
plant = cbind(pops, plant)

### Prepare insect 1 (change ME to MF accordingly, pay attention to ALSO change name in exports)
insect1 = read.delim("ME159ind119loc.txt", header = F, row.names = 1)
insect1 = insect1[ -1, -1]
insect1 = insect1[, -ncol(insect1)]
pops = substr(rownames(insect1), 0, 3)
insect1 = cbind(pops, insect1)
# filter out single guys
counts = table(as.character(pops))
targ = names(counts[counts < 2])
insect1 = insect1[ -match(targ, pops),]

### Insect 2
insect2 = read.delim("MF74ind83loc.txt", header = F, row.names = 1)
insect2 = insect2[ -1, -1]
insect2 = insect2[, -ncol(insect2)]
pops = substr(rownames(insect2), 0, 3)
insect2 = cbind(pops, insect2)
# filter out single guys
counts = table(as.character(pops))
targ = names(counts[counts < 2])
insect2 = insect2[ -match(targ, pops),]



### Prepare geo
# get metadata
info = read.delim("LysiMacro.info", header = T, row.names = 1)
insect1 = insect1[ !is.na(match(insect1$pops, info$Pop)), ]
insect2 = insect2[ !is.na(match(insect2$pops, info$Pop)), ]


### IBD among pops
# GST
plant.mat = plant[, -1]
plant.mat = plant.mat[, colSums(plant.mat) > 0]
plant.pops = plant[, 1]
plant.d = GST_pairwise(plant.mat, plant.pops)

insect1.mat = insect1[, -1]
insect1.mat = insect1.mat[, colSums(insect1.mat) > 0]
insect1.pops = insect1[, 1]
insect1.d = GST_pairwise(insect1.mat, insect1.pops)

insect2.mat = insect2[, -1]
insect2.mat = insect2.mat[, colSums(insect2.mat) > 0]
insect2.pops = insect2[, 1]
insect2.d = GST_pairwise(insect2.mat, insect2.pops)


# Dist geo
geo.plant = info[match(rownames(plant.d), info$Pop), ]

geo.insect1 = info[match(rownames(insect1.d), info$Pop), ]
geo.insect2 = info[match(rownames(insect2.d), info$Pop), ]

geo.insect1.d = COSINE(geo.insect1[, 3], geo.insect1[, 2])
geo.insect2.d = COSINE(geo.insect2[, 3], geo.insect2[, 2])

geo.plant.d = COSINE(geo.plant[, 3], geo.plant[, 2])

NicePlot(as.vector(geo.insect1.d), as.vector(insect1.d), name = "Scatterplot_IBDpops_GST_Geo_ME.pdf")
NicePlot(as.vector(geo.insect2.d), as.vector(insect2.d), name = "Scatterplot_IBDpops_GST_Geo_MF.pdf")
NicePlot(as.vector(geo.plant.d), as.vector(plant.d), name = "Scatterplot_IBDpops_GST_Geo_Lysi.pdf")

################
### new for aggregation

require(tidyverse)

length(intersect(colnames(insect1.d), colnames(plant.d))) # cool

both1 <- intersect(colnames(insect1.d), colnames(plant.d))
both2 <- intersect(colnames(insect2.d), colnames(plant.d))

head(geo.plant.d)

plant.d[lower.tri(plant.d, diag=TRUE)] <- NA
insect1.d[lower.tri(insect1.d, diag=TRUE)] <- NA
insect2.d[lower.tri(insect2.d, diag=TRUE)] <- NA

plant1.d <- plant.d[both1, both1]
insect1.d <- insect1.d[both1, both1]
plant2.d <- plant.d[both2, both2]
insect2.d <- insect2.d[both2, both2]


all1.d <- data.frame(From=rownames(plant1.d),plant1.d) %>% gather("To","Gst",-From) %>% filter(!is.na(Gst))
head(all1.d) # nice

all2.d <- data.frame(From=rownames(plant2.d),plant2.d) %>% gather("To","Gst",-From) %>% filter(!is.na(Gst))
head(all2.d) # nice

all1.d <- data.frame(From=rownames(insect1.d),insect1.d) %>% gather("To","Gst",-From) %>% filter(!is.na(Gst)) %>% left_join(all1.d, by=c("From", "To"), suffix=c(".plant",".insect"))
head(all1.d) # nice

all2.d <- data.frame(From=rownames(insect2.d),insect2.d) %>% gather("To","Gst",-From) %>% filter(!is.na(Gst)) %>% left_join(all2.d, by=c("From", "To"), suffix=c(".plant",".insect"))
head(all2.d) # nice


rownames(info) <- info$Pop

geo1.d <- COSINE(info[both1,3],info[both1,2])
geo2.d <- COSINE(info[both2,3],info[both2,2])


dim(geo1.d)
colnames(geo1.d) <- rownames(geo1.d) <- both1
geo1.d[lower.tri(geo1.d, diag=TRUE)] <- NA

dim(geo2.d)
colnames(geo2.d) <- rownames(geo2.d) <- both2
geo2.d[lower.tri(geo2.d, diag=TRUE)] <- NA

head(geo1.d)
head(geo2.d)

all1.d <- data.frame(From=rownames(geo1.d),geo1.d) %>% gather("To","geo",-From) %>% filter(!is.na(geo)) %>% left_join(all1.d, by=c("From", "To")) %>% select(From, To, Gst.plant, Gst.insect, geo)

all2.d <- data.frame(From=rownames(geo2.d),geo2.d) %>% gather("To","geo",-From) %>% filter(!is.na(geo)) %>% left_join(all2.d, by=c("From", "To")) %>% select(From, To, Gst.plant, Gst.insect, geo)

head(all1.d) # nice
dim(all1.d) # nice
head(all2.d) # nice
dim(all2.d) # nice

cor.test(~Gst.plant+Gst.insect, data=all1.d, method="sp") # welp, it starts here: significant NEGATIVE cor between plant and insect Gst?
cor.test(~Gst.plant+geo, data=all1.d, method="sp") # positive IBD cor
cor.test(~Gst.insect+geo, data=all1.d, method="sp") # n.s. IBD cor

cor.test(~Gst.plant+Gst.insect, data=all2.d, method="sp") # n.s. plant-insect for other sp
cor.test(~Gst.plant+geo, data=all2.d, method="sp") # STRONG positive IBD cor
cor.test(~Gst.insect+geo, data=all2.d, method="sp") # n.s. IBD cor


write.table(all1.d, "Triponez_dists_1.txt", sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE) # gotta look at this again, re-run for second associate

write.table(all2.d, "Triponez_dists_2.txt", sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE) # gotta look at this again, re-run for second associate

