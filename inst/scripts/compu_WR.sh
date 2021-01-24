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

# Retrieve location of the installed R package WeatherRegimes
wr_folder=$(Rscript -e 'system.file(package="WeatherRegimes")' | sed 's/.*"\(.*\)".*/\1/')
echo $wr_folder


## Calcul des representants des regimes de temps
varname="slp"
fname="$wr_folder/extdata/SLP_IPSLCM5MR_19500101_19991231_daily.nc"
seas="DJF"
nreg=4
fout="myfileWR.Rdata"

Rscript --vanilla $wr_folder/scripts/regimes_IPSL.R ${fname} ${seas} ${varname} ${nreg} ${fout} 

## Classification en regimes de temps
varname="slp"
fname="$wr_folder/extdata/SLP_IPSLCM5MR_19500101_19991231_daily.nc"
freg="myfileWR.Rdata"
fout="myfileCl.txt"

Rscript --vanilla $wr_folder/scripts/classif_IPSL.R ${fname} ${varname} ${freg} ${fout} 

end_date=`date +"%m/%d/%Y (%H:%M)"`
echo -e "\nScript finished at: ${end_date}\n"
