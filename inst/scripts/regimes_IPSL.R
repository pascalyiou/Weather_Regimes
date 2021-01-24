## Classification d'une simulation du modele de l'IPSL
## Pascal Yiou (LSCE), Jan. 2021
## Le script peut etre lance en "sh" par
## R CMD "--args fname seas nreg fout" regimes_IPSL.R
## Cette commande peut etre incluse dans un script sh
## Il faut avoir installe les librairies: ncdf4, Mclust, maps
## > install.packages(c("ncdf4","mclust","maps"), dependencies=TRUE)

library(WeatherRegimes)

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
    if(any(is.na(args[seq.int(i)]))) stop("at least one argument is missing")
}else{
    fname=system.file("extdata", "SLP_IPSLCM5MR_19500101_19991231_daily.nc", package="WeatherRegimes")
    seas="DJF"
    varname="slp"
    nreg=4
    fout="myfileWR.Rdata"
}

## Definitions des saisons
env_seas = new.env()
rfile=system.file("scripts", "def_seasons.R", package="WeatherRegimes")
source(rfile, local = env_seas)
l.seas = env_seas$l.seas


## Lecture des donnees a classer dans un fichier ncdf
dat = readnc(varname=varname,fname)
dat.time = dat$time
## Soustraction du cycle saisonnier & calcul d'anomalies saisonnieres
# debug(sousseasmean)
dat.dum=sousseasmean(dat$dat,dat$time)


dat$anom=dat.dum$anom
dat$seascyc=dat.dum$seascyc

## Calcul des poids sur la latitude pour le calcul des PC/EOF
pond.z500=1/sqrt(cos(dat$lat*pi/180))
scale.z500=rep(pond.z500, each=length(dat$lon))

## Pour faciliter la sauvegarde des longitudes et latitudes
lon=dat$lon
lat=dat$lat

## Selection des jours correspondant a la saison seas
iseas=which(dat$time$month %in% l.seas[[seas]])
if(length(iseas) == -1) stop("no data found for the selected season")
dat.m=dat$anom[iseas,]

## Calcul des PCs
pc.dat=prcomp(dat.m,scale.=scale.z500)

## Calcul des regimes et classification
#debug(classnorm)
dat.class=classnorm(pc.dat,nreg=nreg)

## Sauvegarde dans f.out au format Rdat
save(file=fout,dat.class,pc.dat,nreg,fname,seas,varname,lon,lat,dat.time)


