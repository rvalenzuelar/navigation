#!/bin/bash
#
# Bash script that uses Radx application to convert files
# from dorade (sweep) format to CfRadial format
#
#  Raul Valenzuela
#  July, 2104

#inputfiles=/home/rvalenzuela/P3/dorade/case03/leg01
#outputdir=/home/rvalenzuela/P3/CfRadial/case03_new

# inputfiles=/home/rvalenzuela/P3/dorade/case04
# outputdir=/home/rvalenzuela/P3/CfRadial/case04

inputfiles=/home/rvalenzuela/P3/dorade/case04_gpsins_corrected
outputdir=/home/rvalenzuela/P3/dorade/case04_gpsins_corrected

echo 
echo "Converting dorade files from: "  
echo $inputfiles 
echo "CfRadial Output directory: "
echo $outputdir
echo

response=" "

read -r -p "Do you want to continue? [y/N] " response
case $response in
    [yY][eE][sS]|[yY]) 
	echo "Using RadxConvert ..."
        	RadxConvert -params convert_params -f "$inputfiles/swp"* -outdir $outputdir
        	echo
        	;;
    *)
	echo
        ;;
esac

exit



