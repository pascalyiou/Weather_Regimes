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

"classnorm"=function(pc.dat, nreg=4, npc=10, nsim=200)
{
    require(mclust)
    ## Classification par kmeans
    ## On effectue nsim=100 classifications et on recupere les centroides
    if(class(pc.dat)=="princomp"){
        pc = pc.dat$scores[,1:npc]
        rot = pc.dat$loadings[,1:npc]
    }else{
        pc = pc.dat$x[,1:npc]
        rot = pc.dat$rotation[,1:npc]
    }
    ndat=nrow(pc)
    ## Initialisations aleatoires du kmeans      
    dum=matrix(NA, nrow = nsim * nreg, ncol = npc)
    kmeans.dat=kmeans(pc,nreg)
    dum[1:nreg, ]=kmeans.dat$centers
    for(i in 2:nsim){
        kmeans.dat=kmeans(pc[sample(1:ndat,ndat),],nreg)
        dum[(i-1) * nreg + 1:nreg, ]=kmeans.dat$centers
    }
    ## Classification des nreg*nsim centroides par mixture modeling
    dum.mclust=mclust::Mclust(dum) ## Mclust est dans le package mclust
    ## On met la classif par groupes de nreg
    dum.mcl.cla=t(matrix(dum.mclust$classification,nreg,nsim))
    ## On ordonne la classification des classifications
    dum.str=character(nsim)
    for(i in 1:nsim){
      dum.str[i]=paste(sort(dum.mcl.cla[i,]),sep="",collapse="")
    }
    ## On determine la classe de classif la plus probable
    dum.class=table(dum.str)
    dum.levels=names(dum.class)
    class.max=which.max(dum.class)
    ## Determination des manips qui conduisent a cette classif
    I.max=which(dum.str==dum.levels[class.max])
    II=(rep(I.max,each=nreg)-1)*nreg + (1:nreg)
    kmeans.dat=kmeans(pc,dum[II[1:nreg],])

    ## Calcul des regimes dans l'espace physique
    reg.var.kmeans=rot %*% t(kmeans.dat$centers)
    ## Frequence de chaque regime
    perc.r=numeric(nreg)
    for(i in 1:nreg){
      perc.r=c(perc.r,length(which(kmeans.dat$cluster==i))/
                                   length(kmeans.dat$cluster)*100)
    }
    classif.out=list(kmeans=kmeans.dat,reg.var=reg.var.kmeans,
                     perc.r=perc.r)
    invisible(classif.out)
}
