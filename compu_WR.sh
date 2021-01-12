#!/bin/sh -l
## Calcul de regimes de temps et classification sur une sortie de
## modele de l'IPSL
## Pascal Yiou (LSCE) Jan. 2021
## Se lance par
## ./compu_WR.sh

start_date=`date +"%m/%d/%Y (%H:%M)"`
echo -e "\n\nStarting script at: ${start_date}\n"

module load R/4.0.3 ## Optional

## Regimes de temps
varname="slp"
fname=/mypath/myfileref.nc
seas=DJF
nreg=4
fout=/mypath/myfileWR.Rdata

R CMD "--args ${fname} ${seas} ${varname} ${nreg} ${fout}" regimes_IPSL.R

## Classification sur les regimes de temps
fname=/mypath/myfile.nc
freg=/mypath/myfileWR.Rdata
fout=/mypath/myfile_classif.Rdata

R CMD BATCH "--args ${fname} ${varname} ${freg} ${fout}" classif_IPSL.R

end_date=`date +"%m/%d/%Y (%H:%M)"`
echo -e "\nScript finished at: ${end_date}\n"
