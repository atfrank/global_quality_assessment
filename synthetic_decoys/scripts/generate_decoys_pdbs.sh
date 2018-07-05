#!/bin/bash
if [[ $# -ne 1 ]]
then
    echo "usage: $0 <ndecoys>"
else
    source ~/.bashrc
    ndecoys=$1
    ndecoys=`seq 1 ${ndecoys}`

    # script to convert decoys PDB files to DCD
    rnas="1XHP 1YSV 1Z2J 1ZC5 2FDT 2JWV 2K66 2KOC 2L1V 2L3E 2LBJ 2LBL 2LDL 2LDT 2LHP 2LI4 2LK3 2LP9 2LPA 2LPS 2LQZ 2LU0 2LUB 2LUN 2LV0 2M4W 2M5U 2M8K 2M12 2M21 2M22 2M24 2MEQ 2MFD 2MHI 2MNC 2MXL 2N2O 2N2P 2N4L 2N6S 2N6T 2N6W 2N6X 2N7X 2NBY 2NBZ 2NC0 2NCI 2QH2 2QH4 2RVO 2Y95 4A4S 4A4T 4A4U 5A17 5A18 5IEM 5KH8 5KMZ 5KQE 5LSN 5LWJ 5N5C 5UF3 5UZT 5V16 5V17 5WQ1 6EZ0 28SR"
    export vmd="/Applications/VMD.1.9.2.app/Contents/vmd/vmd_MACOSXX86"
    # loop over rnas and get data
    
    mkdir decoy_files
    for rna in ${rnas}
    do
        #catdcd -o decoy_files/${rna}_nlb_decoys.dcd -otype dcd -s pdb_files/${rna}.pdb -stype pdb -pdb decoy_files/${rna}_nlb_decoys.pdb
        convpdb.pl -model 1 decoy_files/${rna}_nlb_decoys.pdb > decoy_files/${rna}_nlb_decoys_reference.pdb
        for i in ${ndecoys}
        do
            catdcd -first ${i} -last ${i} -o tmp.pdb -otype pdb -s decoy_files/${rna}_nlb_decoys_reference.pdb -stype pdb -pdb decoy_files/${rna}_nlb_decoys.pdb
            $vmd -dispdev text -e scripts/build.tcl
            mv tmp_autopsf.pdb decoy_files/${rna}_nlb_decoy_${i}.pdb
            rm -f tmp*pdb tmp*psf tmp*log
        done
        rm -f decoy_files/${rna}_nlb_decoys_reference.pdb
    done
fi