#!/bin/bash

# Bash script to compile the navigation correction code

rm aster2txt cns_noaa-p3_rv netcdf2text 

NETCDF="/usr/local"

FLIB="-L$NETCDF/lib -lnetcdf -lnetcdff -lcurl -lhdf5 -lhdf5_hl"

FINC="-I$NETCDF/include -I/usr/include/geotiff"

CLIB="-L/usr/lib -lgeotiff -ltiff -L/usr/local/lib"

CINC="-I/usr/include/geotiff"

gfortran -g -o cns_noaa-p3_rv cns_noaa-p3_rv.f90 chol_inv.f $FLIB $FINC

gfortran -g -o netcdf2text netcdf2text.f90 $FLIB $FINC

gcc -o aster2txt aster2txt.c $CLIB $CINC

