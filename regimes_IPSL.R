## Classification d'une simulation du modele de l'IPSL
## Pascal Yiou (LSCE), Jan. 2021
## Le script peut etre lance en "sh" par
## R CMD "--args fname seas nreg fout" regimes_IPSL.R
## Cette commande peut etre incluse dans un script sh
## Il faut avoir installe les librairies: ncdf4, Mclust, maps
## > install.packages(c("ncdf4","mclust","maps"), dependencies=TRUE)
library(ncdf4)
Rsource="./"

## Parametres de l'analyse
args=(commandArgs(TRUE))
print(args)
i=1
if(length(args)>0){
    fname=args[i];i=i+1 ## Fichier d'entree
    seas=args[i];i=i+1 ## Saison d'analyse
    varname=args[i];i=i+1 ## Variable a analyser
    nreg=as.integer(args[i]);i=i+1 ## Nombre de regimes
    fout=args[i] ## Fichier de sortie
}else{
    fname="inst/extdata/SLP_IPSLCM5MR_19500101_19991231_daily.nc"
    seas="DJF"
    varname="slp"
    nreg=4
    fout="myfileWR.Rdata"
}

## Initialisation des fonctions
source(paste(Rsource,"imagecont_WR.R",sep="")) ## Pour tracer les resultats
source(paste(Rsource,"read_ncfiles.R",sep="")) ## Lecture des fichiers ncdf
source(paste(Rsource,"preproc_WR.R",sep="")) ## Calcul des cycles saisonniers
source(paste(Rsource,"compu_WR.R",sep="")) ## Calcul des regimes

## Definitions des saisons
env_seas = new.env()
source("def_seasons.R", local = env_seas)
l.seas = env_seas$l.seas


## Lecture des donnees a classer dans un fichier ncdf
datMOD = readipslnc(varname=varname,fname,yr.range=c(1950,1999),ical=365)
dat.MOD.time = datMOD$time
## Soustraction du cycle saisonnier & calcul d'anomalies saisonnieres
# debug(sousseasmean)
dat.MOD.dum=sousseasmean(datMOD$dat,datMOD$time)


datMOD$anom=dat.MOD.dum$anom
datMOD$seascyc=dat.MOD.dum$seascyc

## Calcul des poids sur la latitude pour le calcul des PC/EOF
pond.z500=1/sqrt(cos(datMOD$lat*pi/180))
scale.z500=rep(pond.z500,length(datMOD$lon))

## Selection des jours correspondant a la saison seas
ISEAS=which(datMOD$time$month %in% l.seas[[seas]])
dat.m=datMOD$anom[ISEAS,]

## Calcul des PCs
pc.dat=prcomp(dat.m,scale.=scale.z500)

## Calcul des regimes et classification
#debug(classnorm)
dat.class=classnorm(pc.dat,nreg=nreg,lon=datMOD$lon,lat=datMOD$lat)

## Sauvegarde dans f.out au format Rdat
save(file=fout,dat.class,pc.dat,nreg,fname,seas,dat.MOD.time,ISEAS,
     l.seas,varname)


