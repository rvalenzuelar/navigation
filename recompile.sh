#!/bin/bash

# Bash script to compile the navigation correction code

rm cns_eldo_cai

NETCDF="/usr/local"

FLIB="-L$NETCDF/lib -lnetcdf -lnetcdff -lcurl -lhdf5 -lhdf5_hl"

FINC="-I$NETCDF/include -I/usr/include/geotiff"

gfortran -g -o cns_eldo_cai cns_eldo_cai.f chol_inv.f $FLIB $FINC
