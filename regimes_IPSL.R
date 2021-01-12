## Classification d'une simulation du modele de l'IPSL
## Pascal Yiou (LSCE), Jan. 2021
## Le script peut etre lance en "sh" par
## R CMD "--args fname seas nreg fout" regimes_IPSL.R
## Cette commande peut etre incluse dans un script sh
## Il faut avoir installe les librairies: ncdf4, Mclust, maps
## > install.packages(c("ncdf4","mclust","maps"), dependencies=TRUE)
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
    seas=args[i];i=i+1 ## Saison d'analyse
    varname=args[i];i=i+1 ## Variable a analyser
    nreg=as.integer(args[i]);i=i+1 ## Nombre de regimes
    fout=args[i];i=iBonjour +1 ## Fichier de sortie
}else{
    fname="nom_de_simu.nc"
    seas="DJF"
    varname="slp"
    nreg=4
    fout="nom_de_resultat.Rdata"
}

## Initialisation des fonctions
source(paste(Rsource,"imagecont_WR.R",sep="")) ## Pour tracer les resultats
source(paste(Rsource,"read_ncfiles.R",sep="")) ## Lecture des fichiers ncdf
source(paste(Rsource,"preproc_WR.R",sep="")) ## Calcul des cycles saisonniers
source(paste(Rsource,"compu_WR.R",sep="")) ## Calcul des regimes

## Definitions des saisons
l.seas=list(JJA=6:8,SON=9:11,DJF=c(12,1,2),MAM=c(3,4,5),
            SONDJF=c(9:12,1,2),JJAS=6:9,NDJFM=c(11,12,1,2,3))

## Lecture des donnees a classer dans un fichier ncdf
datMOD = readipslnc(varname=varname,fname,yr.range=c(1950,2000),ical=360))

## Soustraction du cycle saisonnier & calcul d'anomalies saisonnieres
dat.MOD.dum=sousseasmean(datMOD$dat,datMOD$conv.time)

dat.MOD.time = datMOD$time

datMOD$anom=dat.MOD.dum$anom
datMOD$seascyc=dat.MOD.dum$seascyc

## Calcul des poids sur la latitude pour le calcul des PC/EOF
pond.z500=1/sqrt(cos(datMOD$lat*pi/180))
scale.z500=rep(pond.z500,length(datMOD$lon))

## Selection des jours correspondant a la saison seas
ISEAS=which(datMOD$time$month %in% l.seas[[seas]])
dat.m=datMODP$anom[ISEAS,]

## Calcul des PCs
pc.dat=prcomp(dat.m,scale.=scale.z500)

## Calcul des regimes et classification
dat.class=classnorm(pc.dat,nreg=nreg,lon=lon,lat=lat)

## Sauvegarde dans f.out au format Rdat
save(file=fout,dat.class,pc.dat,nreg,fname,seas,dat.MOD.time,ISEAS,
     l.seas,varname)

q("no")

