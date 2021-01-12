## Routines de classification par kmeans avec Monte-Carlo et
## classification des classifications
## Par Pascal Yiou (LSCE). Janvier 2021 (update de la version de 2014).

## Fonction generale de determination des regimes et classification
## Demande d'avoir calcule les PC/EOF du champ a classer
## La procedure est decrite dans l'article:
## Yiou, P.; Goubanova, K.; Li, Z. X. & Nogaj, M. Weather regime dependence of extreme value statistics for summer temperature and precipitation, Nonlinear Processes in Geophysics, 2008, 15, 365-378
## https://hal.archives-ouvertes.fr/hal-00331119/

## Il faut avoir installe le package 'mclust' en amont:
## install.packages("mclust",dependencies=TRUE)

"classnorm"=function(pc.dat,nreg=4,npc=10,nsim=200,simuname="",
                     varname="z500",yr1="",yr2="",lon="",lat="")
{
## Classification par kmeans
## On effectue nsim=100 classifications et on recupere les centroides
    require(mclust)
    if(class(pc.dat)=="princomp"){
        kmeans.dat=kmeans(pc.dat$scores[,1:npc],nreg)
    }else{
        kmeans.dat=kmeans(pc.dat$x[,1:npc],nreg)
    }
    ndat=ifelse(class(pc.dat)=="princomp",nrow(pc.dat$scores),nrow(pc.dat$x))
## Initialisations aleatoires du kmeans      
    dum=kmeans.dat$centers
    for(i in 2:nsim){
        kmeans.dat=kmeans(pc.dat$x[sample(1:ndat,ndat),1:npc],nreg)
        dum=rbind(dum,kmeans.dat$centers)
    }
## Classification des nreg*nsim centroides par mixture modeling
    dum.mclust=Mclust(dum) ## Mclust est dans le package mclust
## On met la classif par groupes de nreg
    dum.mcl.cla=t(matrix(dum.mclust$classification,nreg,nsim))
## On ordonne la classification des classifications
    dum.str=c()
    for(i in 1:nsim){
      s.a=paste(sort(dum.mcl.cla[i,]),sep="",collapse="")
      dum.str=c(dum.str,s.a)
    }
    dum.levels=levels(factor(dum.str))
## On determine la classe de classif la plus probable
    dum.class=c()
    for(lev in dum.levels) dum.class=c(dum.class,length(which(dum.str==lev)))
    class.max=which.max(dum.class)
## Determination des manips qui conduisent a cette classif
    I.max=which(dum.str==dum.levels[class.max])
    II=(rep(I.max,each=nreg)-1)*nreg+c(1:nreg)
    dum.II=dum[II,]
    kmeans.dat=kmeans(pc.dat$x[,1:npc],dum[II[1:nreg],])

## Calcul des regimes dans l'espace physique
    reg.var.kmeans=pc.dat$rotation[,1:npc] %*% t(kmeans.dat$centers[1:nreg,])
## Frequence de chaque regime
    perc.r=c()
    for(i in 1:nreg) perc.r=c(perc.r,length(which(kmeans.dat$cluster==i))/
                                     length(kmeans.dat$cluster)*100)
    classif.out=list(kmeans=kmeans.dat,reg.var=reg.var.kmeans,
                     perc.r=perc.r,lon,lat)
    detach(package:mclust)
    invisible(classif.out)
}#end function
