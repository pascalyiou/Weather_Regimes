# Weather_Regimes

Computation of weather regimes and classification
Pascal Yiou (LSCE), Jan. 2021

This archive contains R scripts to compute weather regimes in a netcdf file (the example is given for a simulation of the IPSL climate model).

There are two steps:
Computation of weather regimes. This is done with multiple kmeans classifications, on which an additional classification is done with a mixture model. This identifies the most probable classification. This is done with the "regimes_IPSL.R" script, preferrably on a long dataset (e.g. control or historical simulation).
The script can be modified to work on reanalyses or other climate model output.

Classification of a dataset on identified weather regimes. The distance to each weather regime is computed, then the closest weather regime is determined. This is done with "classif_IPSL.R", on potentially short datasets.

Those two scripts require an initialization with:
read_ncfiles.R: read netcdf files, extract what is needed, etc.
preproc_WR.R: computation of seasonal cycle and other functions
compu_WR.R: computation of weather regimes with classification
imagecont_WR.R: plotting (optional)

Some paths need to be adapted in the script headers.

A sample sh script is provided to illustrate how it works:
compu_WR.sh

The program is provided "as is", without warranty. It uses the procedure described by:
Yiou, P.; Goubanova, K.; Li, Z. X. & Nogaj, M. Weather regime dependence of extreme value statistics for summer temperature and precipitation, Nonlinear Processes in Geophysics, 2008, 15, 365-378

The scripts can be distributed freely (CECILL license). For a commercial use, please contact Pascal Yiou (pascal.yiou "at" lsce.ipsl.fr).
