## Routines de petits calculs qui aident a la classification
## Pascal Yiou (LSCE), Janvier 2021 (tire d'une version de 2014)

## Soustraction de la moyenne saisonniere jour a jour pour une matrice
## Le cycle saisonnier peut etre estime sur un sous ensemble d'annees
"sousseasmean"<-function(dat,conv.time,l.mon=1:12,
                         l.year=unique(conv.time$year),rprint=FALSE)
{
    dat.m=dat
    seas.cyc=c()
    time.cyc=c()
## Calcul du cycle saisonnier pour une liste d'annees dans l.year
    if(rprint) print("Seasonal cycle")
    for(mon in l.mon){
        r.day=range(conv.time$day[conv.time$month==mon])
        for(day in r.day[1]:r.day[2]){
            m.day=apply(dat[conv.time$day==day & conv.time$month==mon &
                            conv.time$year %in% l.year,],2,mean,na.rm=TRUE)
            seas.cyc=rbind(seas.cyc,m.day)
            time.cyc$month=c(time.cyc$month,mon)
            time.cyc$day=c(time.cyc$day,day)
        }
    }
## Lissage par spline du cycle saisonnier avec preservation des derivees aux bords
    if(rprint) print("Seasonal cycle smoothing")
    seas.cyc.spl=rbind(seas.cyc,seas.cyc,seas.cyc)
    for(i in 1:ncol(seas.cyc)){
        seas.cyc.spl[,i]=smooth.spline(seas.cyc.spl[,i],spar=0.8)$y
    }
    seas.cyc.spl=seas.cyc.spl[(nrow(seas.cyc)+1):(2*nrow(seas.cyc)),]

## Soustraction du cycle saisonnier
    if(rprint) print("Seasonal anomalies")
    for(t in 1:nrow(dat)){
        mon=conv.time$month[t]
        day=conv.time$day[t]
        ii= which(time.cyc$month %in% mon & time.cyc$day %in% day)
        dat.m[t,]=dat[t,]-seas.cyc.spl[ii,]
    }
    datsub=list(anom=dat.m,seascyc=list(seascyc=seas.cyc.spl,timecyc=time.cyc))
    invisible(datsub)
}# fin de la definition de fonction

##---------------------------------------------------------------------------
## Les fonctions qui suivent sont moins utiles
## Classification des regimes par ordre decroissance de frequence
"sort.reg" <- function(Xclust,Xcenter=NULL)
{
## Calcul des frequences de regimes
    freq=c()
    for(i in 1:max(Xclust)) freq=c(freq,length(which(Xclust==i)))
    sort.clust=sort(freq,decreasing=TRUE,index=TRUE)
    sort.clust$ix
}


## Determine l'indice de cluster le plus proche (au sens des moindres carres)
## pour un champ X. reg est issu de kmeans
"which.regime"<-function(X,reg)
{
    regi=list()
    diff=c()
    for(i in c(1:nrow(reg))){
        diffX=t(t(X)-reg[i,])^2
#  diffX=(X-t(matrix(rep(reg[i,],nrow(X)),ncol(X),nrow(X))))^2
        diffXsq=apply(diffX,1,sum)
        diff=cbind(diff,diffXsq)
    }
    regi$cluster=apply(diff,1,which.min)
    regi$dist=apply(diff,1,min)
    return(regi)
} # End function


"classdist" <- function(X,centers)
{
nd=dim(X)
dist=0
for(i in 1:ncol(centers)){
  diffX=t(X)-centers[,i]
  diffsq=apply(diffX^2,2,sum)
  dist=dist+sum(diffsq)
}
return(dist)
}

# Pour un champ, calcule les regimes de jours consecutifs
# qui realisent un critere predefini Ifilt.data
"reg.extr"<-function(Ifilt.data,tt,t.h,clust)
{
    regextr=c()
    for(i in c(1:ncol(Ifilt.data))){
        II=which(Ifilt.data[,i])
        I5=II
        I5=which(clust[II]==clust[II+1] &
                 clust[II]==clust[II+2] &
                 clust[II]==clust[II+3] &
                 clust[II]==clust[II+4] &
                 clust[II]==clust[II+5] &
                 t.h[tt[II+5]]-t.h[tt[II]]==120)
        h=hist(clust[II[I5]],breaks=c(0:4),plot=FALSE)
        regextr=c(regextr,which.max(h$counts))
    }
    return(regextr)
}

# Calcul du cycle saisonnier jour a jour
"seasmean"<-function(dat,conv.time,l.mon=1:12,l.year=unique(conv.time$year),
                     rprint=FALSE)
  {
    seas.cyc=c()
    time.cyc=c()
# Calcul du cycle saisonnier pour une liste d'annees dans l.year
    if(rprint) print("Seasonal cycle")
    for(mon in l.mon){
      r.day=range(conv.time$day[conv.time$month==mon])
      for(day in r.day[1]:r.day[2]){
        m.day=apply(dat[conv.time$day==day & conv.time$month==mon &
          conv.time$year %in% l.year,],2,mean,na.rm=TRUE)
        seas.cyc=rbind(seas.cyc,m.day)
        time.cyc$month=c(time.cyc$month,mon)
        time.cyc$day=c(time.cyc$day,day)
      }
    }
# Lissage par spline du cycle saisonnier
    if(rprint) print("Seasonal cycle smoothing")
    seas.cyc.spl=rbind(seas.cyc,seas.cyc,seas.cyc)
    for(i in 1:ncol(seas.cyc)){
      seas.cyc.spl[,i]=smooth.spline(seas.cyc.spl[,i],spar=0.8)$y
    }
    seas.cyc.spl=seas.cyc.spl[(nrow(seas.cyc)+1):(2*nrow(seas.cyc)),]
    datsub=list(seascyc=seas.cyc.spl,timecyc=time.cyc)
    invisible(datsub)
  }#end function




## Trace des regimes par contours
"regplot.c" <-function(lon,lat,pc.dat,kmeans.dat,title="",lclust=NULL,offset=0)
{
    marg=c(3,3,4,2)
    if(is.null(lclust)) lclust=1:nrow(kmeans.dat$centers)
    nreg=length(lclust)
    reg.var.kmeans=pc.dat$rotation[,1:npc] %*% t(kmeans.dat$centers[lclust,])
    reg.var.kmeans=reg.var.kmeans+offset
    layout(matrix(1:(2*ceiling(nreg/2)),ceiling(nreg/2),2))
    for(i in 1:nreg)
        image.cont.c(lon,lat,reg.var.kmeans[,i],titre=paste("(",letters[i],")",title,"Reg.",lclust[i]),mar=marg,xlab="",ylab="")
}

