## Routines de traces de cartes de champs 2D
## Par Pascal Yiou (LSCE), Jan 2021

## Trace une carte "longitude-latitude" d'un champ, avec la ligne de continents
## qui correspond
"image.cont" <-
    function(lonF,latF,champ,titre="",Ichange=numeric(0),
             zlev=seq(min(champ,na.rm=TRUE),max(champ,na.rm=TRUE),length=11),
             transpose=TRUE,mar=c(5,5,5,6),legend=TRUE,xlab="Longitude",
             ylab="Latitude",satur=FALSE,paquet="fields",add=FALSE,
             lonrange=range(lonF),latrange=range(latF))
{
    col10=rainbow(length(zlev)-1,start=0,end=2/6)
    par( mar=mar)
    if(satur){
        champ[champ<=min(zlev)]=min(zlev)
        champ[champ>=max(zlev)]=max(zlev)
    }
    if(length(Ichange)>0) champ[Ichange]=1
    if (transpose)
        dum=t(matrix(champ,length(latF),length(lonF)))
    else
        dum=matrix(champ,length(lonF),length(latF))
    latF.sort=sort(latF,index.return=TRUE)
    lonF.sort=sort(lonF,index.return=TRUE)
    plot(lonrange,latrange,type="n",xlab=xlab,ylab=ylab,xlim=lonrange,
         ylim=latrange)
    image(sort(lonF),sort(latF),dum[lonF.sort$ix,latF.sort$ix],
          col=col10[length(col10):1],
          xlab=xlab,ylab=ylab,main=titre,breaks=zlev,add=TRUE)
    if(paquet=="fields"){
        library(fields)
        world(add=TRUE)
    }
    if(paquet=="maps"){
        library(maps)
        map(add=TRUE)
    }
    if(legend) image.plot(dum[,length(latF):1],col=col10[length(col10):1],
                          legend.only=TRUE,zlim=range(zlev))
}

##Trace une carte "longitude-latitude" de contours d'un champ, avec la ligne 
##de continents qui correspond
"image.cont.c" <-
    function(lonF,latF,champ,titre="",Ichange=numeric(0),
             zlev=pretty(champ,10),
             transpose=TRUE,mar=c(5,5,5,6),xlab="Longitude",
             ylab="Latitude",col="blue",add=FALSE,paquet="fields",lty=1)
{
    col10=rainbow(length(zlev)-1,start=0,end=2/6)
    par( mar=mar)
    if(length(Ichange)>0) champ[Ichange]=1
    if (transpose)
        dum=t(matrix(champ,length(latF),length(lonF)))
    else
        dum=matrix(champ,length(lonF),length(latF))
    latF.sort=sort(latF,index.return=TRUE)
    nlev=length(zlev)
    contour(lonF,sort(latF),dum[,latF.sort$ix],
            xlab=xlab,ylab=ylab,main=titre,col=col,add=add,nlevels=nlev,
            levels=zlev,lty=lty)
    if(paquet=="fields"){
        library(fields)
        world(xlim=range(lonF),ylim=range(latF),add=TRUE)}
    if(paquet=="maps"){
        library(maps)
        map(add=TRUE)
    }
}

###Trace une carte "longitude-latitude" de contours d'un champ d'anomalies,
## (donc centre en 0) avec la ligne 
## de continents qui correspond
## Les isolignes <0 sont en pointilles et le 0 en gras
## P. Yiou, Oct. 2019.
"image.cont.c.ano" <-
    function(lonF,latF,champ,titre="",Ichange=numeric(0),
             zlev=pretty(champ,10),
             transpose=TRUE,mar=c(5,5,5,6),xlab="Longitude",
             ylab="Latitude",colo="blue",add=FALSE,paquet="fields",lty=1)
{
    par( mar=mar)
    if(length(Ichange)>0) champ[Ichange]=1
    if (transpose){
        dum=t(matrix(champ,length(latF),length(lonF)))
    }else{
        dum=matrix(champ,length(lonF),length(latF))
    }
    latF.sort=sort(latF,index.return=TRUE)
    nlev=length(zlev)
    llty=rep(1,nlev)
    llty[zlev<0]=2
    llwd=rep(1,nlev)
    llwd[zlev==0]=2
    contour(lonF,sort(latF),dum[,latF.sort$ix],
            xlab=xlab,ylab=ylab,main=titre,col=colo,add=add,nlevels=nlev,
            levels=zlev,lty=llty,lwd=llwd)
    if(paquet=="fields"){
        library(fields)
        world(xlim=range(lonF),ylim=range(latF),add=TRUE)}
    if(paquet=="maps"){
        library(maps)
        map(add=TRUE)
    }
}




