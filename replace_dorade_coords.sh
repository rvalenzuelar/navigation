#!/bin/bash


# I/O directories
# INDIR="$HOME/P3/dorade/case04"
# # OUTDIR="$HOME/P3/dorade/case04_coords_cor2"
# OUTDIR="$HOME/P3/dorade/dummy"

# INDIR="$HOME/P3/dorade/case03/leg01"
# OUTDIR="$HOME/P3/dorade/case03_coords_cor"

INDIR="$HOME/P3/dorade/case03/leg03"
OUTDIR="$HOME/P3/dorade/case03_coords_cor/leg03"

# standard tape file
# STDTAPE="$HOME/Github/navigation/010125I.nc"
STDTAPE="$HOME/Github/navigation/010123I.nc"

# python function
PYFUN="$HOME/Github/navigation/replace_cfradial_coords.py"

echo
echo "Changing to input directory: $INDIR"
echo
cd $INDIR
echo 'Running RadxConvert'
RadxConvert -f swp* -cfradial -outdir .

RDXOUT="$(ls -d */)"
echo
echo "Changing to RadxConvert directory: $RDXOUT'"
echo
cd $RDXOUT
echo 'Running replace_cfradial_coords.py'
echo
python $PYFUN $STDTAPE
echo 
echo 'Coordinates replaced'
echo
echo "Cleaning and moving files to $OUTDIR"
mkdir $OUTDIR/cfrad 
mv cfrad.* $OUTDIR/cfrad 
cd $INDIR
rm -rf $RDXOUT
cd $OUTDIR/cfrad
RadxConvert -f cfrad* -dorade -outdir $OUTDIR
cd $OUTDIR
cd $RDXOUT
mv swp* $OUTDIR
cd $OUTDIR
rm -rf $RDXOUT
echo 'Done'
echo









