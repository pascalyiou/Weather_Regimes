#!/bin/sh -l
## Calcul de regimes de temps et classification sur une sortie de
## modele de l'IPSL
## Pascal Yiou (LSCE) Jan. 2021
## Se lance par
## ./compu_WR.sh

set -o errexit
set -o nounset
set -o pipefail

start_date=`date +"%m/%d/%Y (%H:%M)"`
echo -e "\n\nStarting script at: ${start_date}\n"

# module load R/4.0.3 ## Optional

## Regimes de temps
varname="slp"
fname="inst/extdata/SLP_IPSLCM5MR_19500101_19991231_daily.nc"
seas="DJF"
nreg=4
fout="myfileWR.Rdata"

R CMD BATCH $"--args ${fname} ${seas} ${varname} ${nreg} ${fout}" regimes_IPSL.R

## Classification sur les regimes de temps
fname="inst/extdata/SLP_IPSLCM5MR_19500101_19991231_daily.nc"
freg="myfileWR.Rdata"
fout="myfileCl.txt"

R CMD BATCH "--args ${fname} ${varname} ${freg} ${fout}" classif_IPSL.R

end_date=`date +"%m/%d/%Y (%H:%M)"`
echo -e "\nScript finished at: ${end_date}\n"
