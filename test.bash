#!/bin/bash

# compile files decoy charmm trajectories
source ~/.bashrc
cd /Users/atfrank/GitSoftware/global_quality_assessment
pdbs="1KKA 1L1W 1LC6 1LDZ 1NC0 1OW9 1PJY 1R7W 1R7Z 1SCL 1UUU 1XHP 1YSV 1Z2J 1ZC5 2FDT 2JWV 2K66 2KOC 2L3E 2LBJ 2LBL 2LDL 2LDT 2LHP 2LI4 2LK3 2LP9 2LPA 2LQZ 2LU0 2LUB 2LUN 2LV0 2M12 2M21 2M22 2M24 2M4W 2M5U 2M8K 2MEQ 2MFD 2MHI 2MIS 2MNC 2MXL 2N2O 2N2P 2N4L 2N6S 2N6T 2N6W 2N6X 2N7X 2NBY 2NBZ 2NC0 2NCI 2QH2 2QH4 2Y95 4A4S 4A4T 4A4U 5A17 5A18 5KQE"
for pdb in ${pdbs}
do
    mkdir -p reference_errors/${pdb}/            
    cp ~/GitSoftware/assessing_sensivity/decoys/${pdb}/structure_info/reference.pdb reference_errors/${pdb}/.
    cp ~/GitSoftware/assessing_sensivity/decoys/${pdb}/chemical_shifts/chemical_shifts.txt.gz reference_errors/${pdb}/.
    pdb2shifts reference_errors/${pdb}/reference.pdb larmord 0 ${pdb} > reference_errors/${pdb}/chemical_shifts.larmord.txt
    pdb2shifts reference_errors/${pdb}/reference.pdb ramsey 0 ${pdb} > reference_errors/${pdb}/chemical_shifts.ramsey.txt
    gzip -9 reference_errors/${pdb}/chemical_shifts.larmord.txt
    gzip -9 reference_errors/${pdb}/chemical_shifts.ramsey.txt    
done

