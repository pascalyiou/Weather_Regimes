# WeatherRegimes

Computation of weather regimes and classification
Pascal Yiou (LSCE), Jan. 2021

This R package contains functions to compute weather regimes from netcdf file. The exported functions are:
- readnc: reads data from a netcdf file
- sousseasmean: computes the seasonal cycle and the anomalies from a dataset
- classnorm: defines the weather regimes from a principal components analysis

## Installation

It can be installed with the package `remotes` by typing in an R console:

`remotes::install_github("thaos/WeatherRegimes", ref = "main")`

The the package can then be loaded with:

`library(WeatherRegimes)`

The package `remotes` can be installed with:

`install.packages("remotes")`

## Examples

Scripts that serve as examples are provided with this package. The directory contaning  those scripts can be found by typing in a R console the following commands:

`system.file("scripts", package="WeatherRegimes")`

### regimes_IPSL.R 
The script is used to define the weather regimes. This is done with multiple kmeans classifications, on which an additional classification is done with a mixture model. This identifies the most probable classification. It is preferrably applied on a long dataset (e.g. control or historical simulation).

The script takes inputs in the following order:
- The path to a netcdf files used to define the weather regimes.
- The season for which the weather regimes are defined. Available choices are: JJA, SON, DJF, MAM, SONDJF, JJAS or NDJFAM.
- The name of the variable to read from the netcdf. The variable should be a 2d field.
- The number of weather regimes to define.
- The output file name.

The script returns a .Rdata file containing the following variables:
- `dat.class`, results of the function classnorm.
- `pc.dat`, results of the principal component analysis.
- `nreg`, number of regimes.
- `fname`, name of the input netcdf file.
- `seas`, the season.
- `varname`, the variable name.


### classif_IPSL.R 

Classification of a dataset on identified weather regimes. The distance to each weather regime is computed, then the closest weather regime is determined. This can done on potentially short datasets.

The script takes inputs in the following order:
- The path to a netcdf files for which the associated weather regimes are to be computed.
- The name of the variable to read from the netcdf. The variable should be a 2d field.
- The path to an .Rdata file which is the output of `regimes_IPSL.R`.
- The output file name.

The script returns a text file the following columns:
- year
- month
- day
- class, the numbre indentifing the associated weather regime.
- dist, the distance to the associated weather regime.


### compu_WR.sh

A sample bash script to illustrate how to call the scripts `regimes_IPSL.R` 
and `classif_IPSL.R` .

This script can be adapted to your need by changing the arguments of the scripts `regimes_IPSL.R` and `classif_IPSL.R`.


## Non-exported functions

Some functions are present in the R directory of the package but not exported (i.e. hidden from users and undocumented).

They are organized in the following way:

`read_ncfiles.R`: read netcdf files, extract what is needed, etc.

`preproc_WR.R`: computation of seasonal cycle and other functions

`compu_WR.R`: computation of weather regimes with classification

`imagecont_WR.R`: plotting (optional)


## Disclaimer
The program is provided "as is", without warranty. It uses the procedure described by:
Yiou, P.; Goubanova, K.; Li, Z. X. & Nogaj, M. Weather regime dependence of extreme value statistics for summer temperature and precipitation, Nonlinear Processes in Geophysics, 2008, 15, 365-378

The scripts can be distributed freely (CECILL license). For a commercial use, please contact Pascal Yiou (pascal.yiou "at" lsce.ipsl.fr).
