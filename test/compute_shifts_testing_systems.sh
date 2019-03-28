#!/bin/bash
rnas="1R2P 2LPS 2N6Q 5KMZ 2H2X 2M21 2FRL 2M22 2KFC 2L1V"
# download and split pdb 
for rna in $rnas
do
  split_pdb.bash ${rna} ${rna}_model 0
done

# combine old and new structures
declare -a rnas_a=("1R2P" "2N6Q" "2H2X" "2FRL" "2KFC")
declare -a rnas_b=("2LPS" "5KMZ" "2M21" "2M22" "2L1V")
for i in {0..4}
do
  x=${rnas_a[i]}
  y=${rnas_b[i]}
  n_old=`ls ${x}_*.pdb | wc -l`
  n_new=`ls ${y}_*.pdb | wc -l`
  for j in `seq 1 $n_new`
  do
    newIdx=$((j+n_old))
    #echo $newIdx
    mv ${y}_model_${j}.pdb ${y}_model_${newIdx}.pdb
  done
  for j in `seq 1 $n_old`
  do
    mv ${x}_model_${j}.pdb ${y}_model_${j}.pdb
  done
done

# predict chemical shifts using larmord 
rnas="2LPS 5KMZ 2M21 2M22 2L1V"
for rna in $rnas
do
  ndecoys=`ls ${rna}* | wc -l`
  for i in `seq 1 $ndecoys`
  do  
    pdb2shifts ${rna}_model_${i}.pdb larmord ${i} ${rna} | tee -a ${rna}.txt
  done
done