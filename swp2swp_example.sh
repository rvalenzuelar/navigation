#!/bin/bash

#	MBell script for replacing GPS coordinates with GPS-INS
#	corrected coordinates

# make sure there aren't any line continuation characters (ie: \)
export DORADE_DIR='/home/rvalenzuela/P3/dorade/case04_gpsins_corrected'
export AC_NETCDF_FILES='010125I.nc'
export AC_NETCDF_ALIASES='ONLY LATC < lat LONC < lon '
export INPUT_FORMAT='SWEEP_FILES'
export OUTPUT_FLAGS='SWEEP_FILES'

# Optional parameters
#export NCP_THRESHOLD_VAL "NCP 0.25"
#export ALTITUDE_LIMITS "-5. < 22.5"
#export COMPRESSION_SCHEME "HRD_COMPRESSION"
#export FIRST_GOOD_GATE  2
#export DERIVED_FIELDS "DEFAULT(VR > VG  VG > VT)"
#export OPTIONS "DESC_SEARCH"

export BATCH_MODE=' '
export SELECT_RADARS='TA43P3'
xltrsii
export SELECT_RADARS='TF43P3'
xltrsii
