#!/bin/csh

#	MBell script for replacing GPS coordinates with GPS-INS
#	corrected coordinates

# make sure there aren't any line continuation characters (ie: \)
setenv DORADE_DIR "./"
setenv AC_NETCDF_FILES "052403_n42.nc"
setenv AC_NETCDF_ALIASES "LATC < LATC LONC < LONC PALT < PALT HGME < HGM232 GALT < GGALT PITCH < NOT ROLL < NOT VEWC < VEWC VNSC < VNSC THDG < NOT  UIC < UIC VIC < VIC WIC < WIC"
setenv INPUT_FORMAT "SWEEP_FILES"
setenv OUTPUT_FLAGS "SWEEP_FILES"

# Optional parameters
#setenv NCP_THRESHOLD_VAL "NCP 0.25"
#setenv ALTITUDE_LIMITS "-5. < 22.5"
#setenv COMPRESSION_SCHEME "HRD_COMPRESSION"
#setenv FIRST_GOOD_GATE  2
#setenv DERIVED_FIELDS "DEFAULT(VR > VG  VG > VT)"
#setenv OPTIONS "DESC_SEARCH"

setenv BATCH_MODE ""
setenv SELECT_RADARS "TA42P3"
xltrsii
setenv SELECT_RADARS "TF42P3"
xltrsii
