## Routines de lecture de fichiers ncdf
# Pascal Yiou (LSCE)
# Lecture d'un fichier 2D netcdf3
## Modification Juillet 2017 pour lire les fichiers netcdf4
## Mise a jour Janvier 2021
## demande d'avoir fait: library(ncdf4)

## Extraction d'une matrice a partir d'un fichier netcdf
"extractnc" <- function(nc, varnc)
{
    dat=ncvar_get(nc,varnc)
    dat.dim=dim(dat)
    if(length(dat.dim) != 3) stop(paste0("varnc, ", varnc, " must have 3 dimensions: lon, lat, time"))
    nt=dat.dim[3]
    nx=dat.dim[1]
    ny=dat.dim[2]
    ## Remise dans l'ordre temps-lon-lat
    dat=aperm(dat, c(3, 1, 2))
    ## On prefere les tableaux a deux dimensions temps-position
    dim(dat)=c(nt,nx*ny)
    invisible(dat)
}

## Routine de lecture d'un fichier ncdf provenant du modele de l'IPSL
## creation d'un calendrier de 360 jours
"readnc" <- function(varname, fname)
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
    dat=extractnc(nc,varnc)
    ts=nc.get.time.series(nc)
    cal=attr(ts, "cal")
    nc = nc_close(nc) #open.ncdf(fname)
    yyyymmdd=as.character(ts, format="%Y%m%d")
    if(cal == "360"){
      yyyymmdd=cal360_30dbm(yyyymmdd)
    }
    time=data.frame(
      year=as.numeric(substr(yyyymmdd, 1, 4)),
      month=as.numeric(substr(yyyymmdd, 5, 6)),
      day=as.numeric(substr(yyyymmdd, 7, 8))
    )
    data.nc=list(lon=lon,lat=lat,dat=dat,time=time)
    invisible(data.nc)
}
## end function

# TODO: need to be tested.
"cal360_30dbm" <- function(yyyymmdd_PCICt){
  day360=rep(1:30,times=12)
  month360=rep(1:12, each=30)
  dayPCICt=c(1:30, 1:28, 1:31, 1:30, 1:30, 1:30, 1:30, 1:31, 1:30, 1:30, 1:30, 1:30)
  monthPCICt=rep(1:12,  c(30, 28, 31, 30, 30,30, 30, 31, 30, 30, 30, 30))
  convtable=data.frame(
   "PCICt" = sprintf("%02d%02d", monthPCICt, dayPCICt),
   "cal360" = sprintf("%02d%02d", month360, day360)
  )
  mmdd_PCICt=substr(yyyymmdd_PCICt, 5, 8)
  yyyy_PCICt=substr(yyyymmdd_PCICt, 1, 4)
  imatch = match(mmdd_PCICt, convtable$PCICt)
  if(any(is.na(imatch))) stop("at least one date with no correspondance date found")
  return(paste0(yyyy_PCICt, convtable$cal360[imatch]))
}

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
