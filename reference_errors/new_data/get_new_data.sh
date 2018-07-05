#!/bin/bash

# script to add new data to reference error
rnas="5KQE 5IEM 5LWJ 2RVO 5UZT 5UF3 5WQ1 5V17 5V16 5KH8 5LSN 6EZ0 5N5C"
for rna in $rnas
do
    rm -rf ${rna}/
    mkdir ${rna}
    cd ${rna}
    
    # get coordinate files
    split_pdb.bash ${rna} ${rna} 0
    mv ${rna}_1.pdb reference.pdb
    rm ${rna}*pdb
    
    # get shifts
    downloadSTR ${rna} > chemical_shifts.txt
    gzip -9 chemical_shifts.txt
    
    # predict shifts
    pdb2shifts reference.pdb larmord 0 ${rna} > chemical_shifts.larmord.txt
    gzip -9 chemical_shifts.larmord.txt
    cd ..
done