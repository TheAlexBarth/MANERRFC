! Warning This is an on-going project with incomplete analyses

# MANERRFC

Processing of Mission Aransas NERR flowcam data

# Getting started

To work through these scripts, follow the numeric flow of the files for the R scripts. Many of these are data generating and need not be ran multiple times. I try to write code to function as units - scripts should be able to be run from terminal.

Data products correspond to the script which produced them.

# R/

-   00: import ecotaxa data and inital formatting
-   01: import environmental
-   02: merge environmental to scales matching etx data
    -   This is really just to keep code clean down the line.

# Stan/

A generalized stan script to support the mixture model.

# Data/

All data are stored in the data folder. Note there are some preexisting (non-numbered) data pieces. SWMP data came from pre-processed set of data exported from CDMO. PDSI data was downloaded from {INSERT?}