#!/bin/bash

## Requires mcverter, command line tool provided by MRIConvert, 
## http://lcni.uoregon.edu/~jolinda/MRIConvert/

command -v mcverter >/dev/null 2>&1 || { echo >&2 "I require mcverter but it's not installed.  Aborting."; exit 1; }

dataDirectory=../data-imaging
convertedSubdirectory="Mcverter_Dicom_conversion"

## find all grandchildren subdirectories of the data directory, like 
##   ./data-imaging/Lin002/20131022
dirs=$(find -L ${dataDirectory%%/} -maxdepth 2 -mindepth 2 -type d)

for dir in $dirs 
do 
    echo "Processing $dir"
    thisDirConverted=${dir%%/}/$convertedSubdirectory
    if [ -d "$thisDirConverted" ]; then 
        echo "  Converted subdirectory already exists"
    else
        numBoldDirs=$(find ${dir%%/} -type d -regex ".*ep2d_bold.*" | wc -w)
        echo "  Creating subdirecotry $thisDirConverted"
        mkdir $thisDirConverted
        echo "  Running mcverter..."
        mcverter --output $thisDirConverted --nii --format nifti --split_dir ${dir%%/}/*/*.dcm
    fi
done
