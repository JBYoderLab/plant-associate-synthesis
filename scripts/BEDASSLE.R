# BEDASSLE analysis of plant-associate data, prototyping
# Noninteractive; converting from hulk to MAJEL environment
# jby 2020.02.07

# starting up ------------------------------------------------------------
# R --vanilla --args Spr.Assn2 data all 10000 no

require("MASS")
require("BEDASSLE",lib.loc="../shared/Rscripts")

source("../shared/Rscripts/base.R") # my special mix of personal functions

# parse the command line ------
ftag = commandArgs(trailingOnly=T)[1] # what's the fileset?
datdir = commandArgs(trailingOnly=T)[2] # where are they? 
e = commandArgs(trailingOnly=T)[3] # what env?
runfor = as.numeric(commandArgs(trailingOnly=T)[4]) # how long to run?
con = commandArgs(trailingOnly=T)[5] == "continue" # extend prior run from exisiting output files?
# ftag = "Spr.Assn2"
# ftag = "Pine.Assn2"

BPbase = paste(datdir,"/",ftag,".test", sep="") 
samfile = paste(datdir,"/",ftag,".samples", sep="") 
envfile = paste(datdir,"/",ftag,".baypass.env", sep="") 
envnames = paste(datdir,"/", ftag,".baypass.env-names", sep="") # think this will always work

#-------------------------------------------------------------------------
# read in data 
# population and environmental data
samples = read.csv(samfile, h=F, sep="")$V1
pops = read.csv(paste(datdir,"/",ftag,".baypass.pops",sep=""), h=F, sep="\t", col.names = c("id","pop"))[,c(2,1)]
envs = data.frame(pop=as.numeric(unique(pops$pop)), t(read.csv(envfile, h=F, s="", row.names=read.csv(envnames, h=F, s="")$V1)))
# allele count data
bp = data.frame(pop=rep(unique(pops$pop), each=2), al=rep(c(1,2),length(unique(pops$pop))), t(read.csv(paste(BPbase,".baypass",sep=""), h=F, sep="\t")))
bp[1:5,1:10]
dim(bp)

# if you want, reduce data for prototyping
# bp = bp[1:40,1:502]

#-------------------------------------------------------------------------
# format for BEDASSLE
# allele counts
cts = as.matrix(subset(bp, al==1)[,3:ncol(bp)])
dimnames(cts) = list(unique(bp$pop), read.csv(paste(BPbase,".baypass.SNPinfo",sep=""),h=T,sep="")$probeset_id)

# sample sizes
ss = as.matrix(subset(bp, al==1)[,3:ncol(bp)]) + as.matrix(subset(bp, al==2)[,3:ncol(bp)])
dimnames(ss) = list(unique(bp$pop), read.csv(paste(BPbase,".baypass.SNPinfo",sep=""),h=T,sep="")$probeset_id)

# toss sites missing lat-long coordinates
nolalo = (is.na(envs$LAT) | is.na(envs$LONG))

if(length(which(nolalo))>0){
cat(paste("DROPPING data from",length(which(nolalo)),"sites missing geographic data.\n\n"))
cts = cts[!nolalo, ]
ss = ss[!nolalo, ]
envs = envs[!nolalo, ]
}

# geo distance
lalo = as.matrix(envs[,c("LAT","LONG")]); rownames(lalo) = envs$pop
Gdist = sapply(as.character(envs$pop), function(a) sapply(as.character(envs$pop), function(b) if(a==b) 0 else gcd.slc(lalo[a,2],lalo[a,1],lalo[b,2],lalo[b,1]))) # boom

# env distance
# (which env?)

if(e!="all"){

if(!(e%in%colnames(envs))) stop(paste("Error: Environment code does not match any ENV in data. Specify one of:",paste(colnames(envs)[-1], collapse=", ")))

en = envs[,e]; names(en) = envs$pop
en = round((en - mean(en))/sd(en),6) # standardize to mean=0, sd=1
Edist = sapply(as.character(unique(bp$pop)), function(a) sapply(as.character(unique(bp$pop)), function(b) if(a==b) 0 else abs(en[a]-en[b])))
}

if(e=="all"){

envsStd = apply(envs[,-1], 2, function(x) round((x - mean(x))/sd(x),6))

envsList = lapply(colnames(envsStd), function(x) sapply(envsStd[,x], function(a) abs(a-envsStd[,x])) )
names(envsList) = colnames(envsStd)

}


#-------------------------------------------------------------------------
# run BEDASSLE
# see help entry for `MCMC()`

# recover end-state parameters from prior run, if continuing
if(con){

confile = sort(list.files("output", pattern=paste(ftag,"_",e,"_\\d+_MCMC_output1.Robj",sep=""),full.names=T),d=T)[1]
conit = as.numeric(gsub(paste("output/",ftag,"_",e,"_(\\d+)_MCMC_output1.Robj",sep=""), "\\1", confile))+1
cpfile = gsub(paste("output/(",ftag,"_",e,"_\\d+_MCMC_output1).Robj",sep=""), "\\1_cp.Robj", confile)

make.continuing.params(confile, paste("output/",cpfile,sep=""))

cat(paste("CONTINUING from", confile, "; new iteration is index", conit, "\n\n"))

}

# The value of delta may set off warnings, 
#so temporarily disable warnings.
op <- options("warn")
options(warn = -1)

# Call the Markov chain Monte Carlo for the standard model ... BB because baseline breaks!
# allowing for different tuning for each species ...
if(ftag=="Spr.Assn2"){
MCMC_BB(
	counts = cts,
	sample_sizes = ss,
	D = round(Gdist/1000,5),
	E = if(e=="all") envsList else round(Edist,5),
	k = nrow(cts),
	loci = ncol(cts),
	delta = 0.01,
	aD_stp = 0.07,
	aE_stp = 0.002,
	a2_stp = 0.0018,
	phi_stp = 2,
	thetas_stp = 0.075,
	mu_stp = 0.2,
	ngen = runfor,
	printfreq = 2,
	savefreq = 100,
	samplefreq = 10,
	directory = "output",
	prefix = if(con) paste(ftag,"_",e,"_",conit,"_" ,sep="") else paste(ftag,"_",e,"_1_",sep=""), # output will have `_MCMC_output1.Robj` appended
	continue = con,
	continuing.params = if(con) cpfile else NULL
)
}
if(ftag=="Pine.Assn2"){
MCMC_BB(
	counts = cts,
	sample_sizes = ss,
	D = round(Gdist/1000,5),
	E = if(e=="all") envsList else round(Edist,5),
	k = nrow(cts),
	loci = ncol(cts),
	delta = 0.01,
	aD_stp = 0.05,
	aE_stp = 0.005,
	a2_stp = 0.001,
	phi_stp = 1.9,
	thetas_stp = 0.075,
	mu_stp = 0.1,
	ngen = runfor,
	printfreq = 2,
	savefreq = 100,
	samplefreq = 10,
	directory = "output",
	prefix = if(con) paste(ftag,"_",e,"_",conit,"_" ,sep="") else paste(ftag,"_",e,"_1_",sep=""), # output will have `_MCMC_output1.Robj` appended
	continue = con,
	continuing.params = if(con) cpfile else NULL
)
}

# Re-enable warnings
options(op)
# Reset wd
setwd("..")

q()