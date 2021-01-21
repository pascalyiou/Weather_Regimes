## Classification d'une simulation du modele de l'IPSL
## Pascal Yiou (LSCE), Jan. 2021
## Peut se lancer en batch par:
## R CMD BATCH "--args fname freg fout" classif_IPSL.R
library(Weather_Regimes)

## Parametres de l'analyse
args=(commandArgs(TRUE))
print(args)
i=1
if(length(args)>0){
    fname=args[i];i=i+1 ## Fichier d'entree
    varname=args[i];i=i+1 ## Variable d'analyse
    freg=args[i];i=i+1 ## Fichier des regimes de reference
    fout=args[i];i=i+1 ## Fichier de sortie
}else{
    fname="inst/extdata/SLP_IPSLCM5MR_19500101_19991231_daily.nc"
    varname="slp"
    freg="myfileWR.Rdata"
    fout="myfileCl.txt"
}

## Definitions des saisons
env_seas = new.env(parent = emptyenv())
source(system.file("inst", "scripts", "def_seasons.R", package="Weather_Regime"), local = env_seas)
l.seas = env_seas$l.seas

## Lecture des regimes sur la simulation de controle
# save(file=fout,dat.class,pc.dat,nreg,fname,seas,dat.MOD.time,ISEAS,
env_inputs = new.env(parent = emptyenv())
load(freg, envir = env_inputs)
dat.class = env_inputs$dat.class
seas = env_inputs$seas
nreg = env_inputs$nreg

## Lecture des données de pression/Z500 sur les autres simulations
dat.IPSL = readnc(fname=fname, varname=varname)

## "readipslnc" <- function(varname="t2m",fname,yr.range,ical=360)

## Soustraction du cycle saisonnier & calcul d'anomalies saisonnieres
dat.IPSL.dum=sousseasmean(dat.IPSL$dat,dat.IPSL$time)

dat.IPSL.time = dat.IPSL$time

dat.IPSL$anom=dat.IPSL.dum$anom
dat.IPSL$seascyc=dat.IPSL.dum$seascyc

I.seas=which(dat.IPSL.time$month %in% l.seas[[seas]])

## Calcul des distances a chaque WR identifie (nreg)
Xdiff=matrix(NA, nrow=length(I.seas), ncol=nreg)
for(i in 1:nreg){
  dum=t(t(dat.IPSL$anom[I.seas,])-dat.class$reg.var[,i])
  dum=dum^2
  Xdiff[, i]=apply(dum,1,sum)/ncol(dum)
}
## Determination du regime le plus proche: classification
class.Xdiff=apply(Xdiff,1,which.min)
dist.reg=sqrt(apply(Xdiff,1,min)/nrow(dat.class$reg.var))


Xout = cbind(dat.IPSL.time[I.seas,], "class"  = class.Xdiff, "dist" = dist.reg)
write.table(file=fout,Xout,row.names=FALSE, quote=FALSE,col.names=TRUE)


