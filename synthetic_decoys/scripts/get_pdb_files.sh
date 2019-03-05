#!/bin/bash
# script to loop over PDB entries and get restraint files
rnas="1XHP 1YSV 1Z2J 1ZC5 2FDT 2JWV 2K66 2KOC 2L3E 2LBJ 2LBL 2LDL 2LDT 2LHP 2LI4 
      2LK3 2LP9 2LPA 2LQZ 2LU0 2LUB 2LUN 2LV0 2M4W 2M5U 2M8K 2M12 2M24 2MEQ 2MFD 
      2MHI 2MNC 2MXL 2N2O 2N2P 2N4L 2N6S 2N6T 2N6W 2N6X 2N7X 2NBY 2NBZ 2NC0 2NCI 
      2QH2 2QH4 2RVO 2Y95 4A4S 4A4T 4A4U 5A17 5A18 5IEM 5KH8 5KQE 5LSN 5LWJ 5N5C 
      5UF3 5UZT 5V16 5V17 5WQ1 6EZ0 28SR"
# loop over rnas and get data
mkdir pdb_files
for rna in ${rnas}
do
	split_pdb.bash ${rna} ${rna} 0
	mv ${rna}_1.pdb pdb_files/${rna}.pdb
	rm -v ${rna}*pdb	
done