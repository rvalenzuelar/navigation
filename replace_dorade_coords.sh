#!/bin/bash


INDIR="/home/rvalenzuela/P3/dorade/case04"
OUTDIR="/home/rvalenzuela/P3/dorade/case04_coords_cor"
STDTAPE="/home/rvalenzuela/Github/navigation/010125I.nc"
PYFUN="/home/rvalenzuela/Github/navigation/replace_cfradial_coords.py"

echo
echo "Changing to $INDIR"
echo
cd $INDIR
echo 'Running RadxConvert'
RadxConvert -f swp* -cfradial -outdir .

RDXOUT="$(ls -d */)"
echo
echo "Changing to $RDXOUT'"
echo
cd $RDXOUT
echo 'Running replace_cfradial_coords.py'
echo
python $PYFUN $STDTAPE
echo 
echo 'Done'
echo










