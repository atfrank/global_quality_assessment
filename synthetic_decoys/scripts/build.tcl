package require autopsf
mol new tmp.pdb waitfor all
set sel [atomselect top "all"]
autopsf -mol 0 -nucleic $sel
exit
