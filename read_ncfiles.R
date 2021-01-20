## Routines de lecture de fichiers ncdf
# Pascal Yiou (LSCE)
# Lecture d'un fichier 2D netcdf3
## Modification Juillet 2017 pour lire les fichiers netcdf4
## Mise a jour Janvier 2021
## demande d'avoir fait: library(ncdf4)

## Lecture des donnees, longitude-latitude, temps, scaling
"lirevarnc" <- function(nc,varname,n.lon="lon",n.lat="lat",n.time="time",
  scaling=c(1,0),y.ref=1800)
{
    varnc=nc$var[[varname]]
    varsize=varnc$varsize
    ndims=varnc$ndims
## Extraction des longitudes, latitudes et du temps
    nclon=nc$dim[[n.lon]]
    lon.dat=nclon$vals
    nclat=nc$dim[[n.lat]]
    lat.dat=nclat$vals
## Traitement du temps
    nctime = nc$dim[[n.time]]
    time=nctime$vals
    conv.time.dat=caldat(time/24+julday(1,1,y.ref))
## Conversion en yyyymmdd    
    time.dat=conv.time.dat$year*10000+conv.time.dat$month*100+
      conv.time.dat$day
## Lecture du fichier netcdf
    dat.m=extractnc(nc,varnc)
## Scaling des donnees    
    dat.m=dat.m/scaling[1]-scaling[2]
    datout=list(dat=dat.m,lon=lon.dat,lat=lat.dat,
      time=time.dat,conv.time=conv.time.dat,timenc=time)
    rm(dat.m)
    invisible(datout)
}

## Extraction d'une matrice a partir d'un fichier netcdf
"extractnc" <- function(nc, varnc)
{
## Lecture des donnees de la saison ISEAS pour l'annee year
    dat=ncvar_get(nc,varnc)
    dat.dim=dim(dat)
    if(length(dat.dim) != 3) stop(paste0("varnc, ", varnc, " must have 3 dimensions: lon, lat, time"))
    if(is.null(ISEAS)){
      ISEAS=1:dat.dim[3]
    }
    nt=length(ISEAS)
    nx=dat.dim[1]
    ny=dat.dim[2]
    dat=dat[,,ISEAS]
    ## Remise dans l'ordre temps-lon-lat
    dat=aperm(dat, c(3, 1, 2))
    ## On prefere les tableaux a deux dimensions temps-position
    dim(dat)=c(nt,nx*ny)
    invisible(dat)
}

## Routine de lecture d'un fichier ncdf provenant du modele de l'IPSL
## creation d'un calendrier de 360 jours
"readipslnc" <- function(varname="t2m",fname,yr.range,ical=360)
{
    nc = nc_open(fname) #open.ncdf(fname)
    varnc=nc$var[[varname]]
    varsize=varnc$varsize
    ndims=varnc$ndims
## Extraction des longitudes, latitudes et du temps
    nclon=nc$dim[['lon']]
    lon=nclon$vals
    nclat=nc$dim[['lat']]
    lat=nclat$vals
##Traitement du temps: creation d'un calendrier a 360 ou 365 jours
    itime=nc$dim[["time_counter"]]$vals
    dat=extractnc(nc,varnc,NULL)
    nc_close(nc)
    if(ical==360){
        ndyear=360
        nmyear=12
        lmo=1:12

        nyear=yr.range[2]-yr.range[1]+1
        day=rep(1:30,times=nmyear*nyear)
        month=rep(rep(lmo,each=30),times=nyear)
        year=rep(yr.range[1]:yr.range[2],each=ndyear)

    }
    if(ical==365){
        year=month=day=integer(nrow(dat))
        monthdum=c(rep(1,31),rep(2,28),rep(3,31),rep(4,30),rep(5,31),rep(6,30),
                   rep(7,31),rep(8,31),rep(9,30),rep(10,31),rep(11,30),
                   rep(12,31))
        daydum=c((1:31),(1:28),(1:31),(1:30),(1:31),(1:30),(1:31),(1:31),
        (1:30),(1:31),(1:30),(1:31))
        ndyear=length(daydum)

        nyear=yr.range[2]-yr.range[1]+1
        year=rep(yr.range[1]:yr.range[2],each=ndyear)
        month=rep(monthdum,times=nyear)
        day=rep(daydum,times=nyear)
    }
    conv.time=data.frame(year=year, month=month, day=day)

    data.nc=list(lon=lon,lat=lat,dat=dat,time=conv.time)
    invisible(data.nc)
}
## end function




## The function computes month, day, and year from Julian days. 
"caldat" <- function (julian) 
{
    igreg = 2299161
    julian <- trunc(julian)
    jalpha <- julian * 0
    ja <- julian * 0
    im <- (julian >= igreg)
    if (sum(im) > 0) {
        jalpha[im] <- trunc(((julian - 1867216) - 0.25)/36524.25)
        ja[im] <- julian + 1 + jalpha - trunc(0.25 * jalpha)
    }
    im <- (julian < igreg)
    if (sum(im) > 0) 
        ja[im] <- julian[im]
    jb <- ja + 1524
    jc <- trunc(6680 + ((jb - 2439870) - 122.1)/365.25)
    jd <- 365 * jc + trunc(0.25 * jc)
    je <- trunc((jb - jd)/30.6001)
    id <- jb - jd - trunc(30.6001 * je)
    mm <- je - 1
    im <- (mm > 12)
    if (sum(im) > 0) 
        mm[im] <- mm[im] - 12
    iyyy <- jc - 4715
    im <- (mm > 2)
    if (sum(im) > 0) 
        iyyy[im] <- iyyy[im] - 1
    im <- (iyyy <= 0)
    if (sum(im) > 0) 
        iyyy <- iyyy - 1
    caldat <- list(month = mm, day = id, year = iyyy)
    invisible(caldat)
}

## The function computes Julian days from month, day, and year. 
"julday" <- function (mm, id, iyyy)
{
    igreg <- 588829
    mm <- trunc(mm)
    id <- trunc(id)
    iyyy <- trunc(iyyy)
    im <- (iyyy == 0)
    if (sum(im, na.rm = TRUE) > 0)
        return("There is no year zero!")
    if ((length(mm) != length(id)) | (length(mm) != length(iyyy)) |
        (length(iyyy) != length(id)))
        return("The vectors must have same length!")
    im <- (iyyy < 0)
    if (sum(im) > 0)
        iyyy[im] <- iyyy[im] + 1
    jy <- mm * 0
    jm <- mm * 0
    ja <- mm * 0
    im <- (mm > 2)
    if (sum(im) > 0) {
        jy[im] <- iyyy[im]
        jm[im] <- mm[im] + 1
    }
    im <- (mm <= 2)
    if (sum(im) > 0) {
        jy[im] <- iyyy[im] - 1
        jm[im] <- mm[im] + 13
    }
    jul <- trunc(365.25 * jy) + trunc(30.6001 * jm) + id + 1720995
    im <- (id + 31 * (mm + 12 * iyyy) >= igreg)
    if (sum(im) > 0) {
        ja[im] <- trunc(0.01 * jy)
        jul[im] <- jul + 2 - ja[im] + trunc(0.25 * ja[im])
    }
    julday <- jul
    invisible(julday)
} 
