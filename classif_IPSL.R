## Classification d'une simulation du modele de l'IPSL
## Pascal Yiou (LSCE), Jan. 2021
## Peut se lancer en batch par:
## R CMD BATCH "--args fname freg fout" classif_IPSL.R

library(ncdf4)

SI=Sys.info()
if(SI[[1]] == "Darwin"){
  Rsource="/Users/yiou/programs/RStat/WREGIMES/"
}
if(SI[[1]] == "Linux"){
  Rsource="/home/users/yiou/RStat/WREGIMES/"
}

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
    fname="nom_de_simu.nc"
    varname="slp"
    freg="nom_de_regimes.Rdata"
    fout="nom_de_resultat.Rdata"
}

## Initialisation des fonctions
source(paste(Rsource,"read_ncfiles.R",sep="")) ## Lecture des fichiers ncdf
source(paste(Rsource,"preproc_WR.R",sep="")) ## Calcul des cycles saisonniers

## Lecture des regimes sur la simulation de controle
load(freg)

## Lecture des données de pression/Z500 sur les autres simulations
yr1=1996
yr2=2005
dat.IPSL = readipslnc(fname=fname, varname=varname,ical=360,
                      yr.range=c(yr1,yr2))

## "readipslnc" <- function(varname="t2m",fname,yr.range,ical=360)

## Calcul d'une ponderation par la racine du cos de la latitude
pond.z500=1/sqrt(cos(dat.IPSL$lat*pi/180))
scale.z500=rep(pond.z500,length(dat.IPSL$lon))

## Soustraction du cycle saisonnier & calcul d'anomalies saisonnieres
dat.IPSL.dum=sousseasmean(dat.IPSL$dat,dat.IPSL$conv.time)

dat.IPSL.time = dat.IPSL$time

dat.IPSL$anom=dat.IPSL.dum$anom
dat.IPSL$seascyc=dat.IPSL.dum$seascyc

I.seas=which(dat.IPSL.time$month %in% l.seas[[seas]])

## Calcul des distances a chaque WR identifie (nreg)
Xdiff=c()
for(i in 1:nreg){
  dum=t(t(dat.IPSL$anom[I.seas,])-dat.class$reg.var[,i])
  dum=dum^2
  sdum=apply(dum,1,sum)/ncol(dum)
  Xdiff=cbind(Xdiff,sqrt(sdum))
}
## Determination du regime le plus proche: classification
class.Xdiff=apply(Xdiff,1,which.min)
dist.reg=sqrt(apply(Xdiff,1,min)/nrow(dat.class$reg.var))


Xout = cbind(dat.IPSL.time[I.seas], class.Xdiff, dist.reg)
write.table(file=fout,Xout,row.names=FALSE, quote=FALSE,col.names=FALSE)

q("no")

