require(ade4)
require(maps)
require(geosphere)
require(vegan)

rm(list=ls())

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

### Prepare insect (change ME to MF accordingly, pay attention to ALSO change name in exports)
insect = read.delim("ME159ind119loc.txt", header = F, row.names = 1)
# insect = read.delim("MF74ind83loc.txt", header = F, row.names = 1)
insect = insect[ -1, -1]
insect = insect[, -ncol(insect)]
pops = substr(rownames(insect), 0, 3)
insect = cbind(pops, insect)
# filter out single guys
counts = table(as.character(pops))
targ = names(counts[counts < 2])
insect = insect[ -match(targ, pops),]


### Prepare geo
# get metadata
info = read.delim("LysiMacro.info", header = T, row.names = 1)
insect = insect[ !is.na(match(insect$pops, info$Pop)), ]


### IBD among pops
# GST
plant.mat = plant[, -1]
plant.mat = plant.mat[, colSums(plant.mat) > 0]
plant.pops = plant[, 1]
plant.d = GST_pairwise(plant.mat, plant.pops)

insect.mat = insect[, -1]
insect.mat = insect.mat[, colSums(insect.mat) > 0]
insect.pops = insect[, 1]
insect.d = GST_pairwise(insect.mat, insect.pops)

# Dist geo
geo.plant = info[match(rownames(plant.d), info$Pop), ]
geo.insect = info[match(rownames(insect.d), info$Pop), ]
geo.insect.d = COSINE(geo.insect[, 3], geo.insect[, 2])
geo.plant.d = COSINE(geo.plant[, 3], geo.plant[, 2])

NicePlot(as.vector(geo.insect.d), as.vector(insect.d), name = "Scatterplot_IBDpops_GST_Geo_ME.pdf")
NicePlot(as.vector(geo.plant.d), as.vector(plant.d), name = "Scatterplot_IBDpops_GST_Geo_Lysi.pdf")



