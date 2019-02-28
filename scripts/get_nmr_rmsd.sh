#!/bin/bash

# script uses convpdb.pl from MMTSB to calculate RMSDs of structure in the combined
# obsolete and updated ensemble
# Uses info from tests_info figure out which reference structure to use.

for i in {1..5}
do
    pdb=`awk -v i=${i} '{if (NR==i) print $1}' data/tests_info`
    ref=`awk -v i=${i} '{if (NR==i) print $4+1}' data/tests_info`
    refpdb="coors/${pdb}_${ref}.pdb"
    for j in {1..50}
    do
        comppdb="coors/${pdb}_${j}.pdb"
        if [[  -f ${comppdb} ]]
        then
            rmsd=`rms.pl -fitsel heavy  ${refpdb} ${comppdb} | awk '{print $1}'`
            echo "${j} ${comppdb} ${rmsd}"
        fi
    done
done