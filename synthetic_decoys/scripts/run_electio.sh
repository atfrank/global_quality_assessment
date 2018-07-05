#!/bin/bash

# script to run Electio for paper
source ~/.bashrc
rnas="5KH8 1XHP 1YSV 1Z2J 1ZC5 2FDT 2JWV 2K66 2KOC 2L1V 2L3E 2LBJ 2LBL 2LDL 2LDT 2LHP 2LI4 2LK3 2LP9 2LPA 2LPS 2LQZ 2LU0 2LUB 2LUN 2LV0 2M4W 2M5U 2M8K 2M12 2M21 2M22 2M24 2MEQ 2MFD 2MHI 2MNC 2MXL 2N2O 2N2P 2N4L 2N6S 2N6T 2N6W 2N6X 2N7X 2NBY 2NBZ 2NC0 2NCI 2QH2 2QH4 2RVO 2Y95 4A4S 4A4T 4A4U 5A17 5A18 5IEM 5KH8 5KMZ 5KQE 5LSN 5LWJ 5N5C 5UF3 5UZT 5V16 5V17 5WQ1 6EZ0 28SR"
for rna in $rnas
do
    # make data ready for running electio
    cat chemical_shift_files/${rna}.txt | awk '{print $6, $2, $3, $4, $5, $1}' > chemical_shift_files/tmp.txt

    electio_prep \
    --accuracyfile ../data/larmord_accuracy_resname_nucleus.txt \
    --observed_vector_filename electio_selections/observed_electio_vector_${rna}.txt \
    --predicted_matrix_filename electio_selections/predicted_matrix_${rna}.txt \
    --weights_vector_filename electio_selections/weights_vector_${rna}.txt \
    ../reference_errors/${rna}/chemical_shifts_corrected_larmord.txt \
    chemical_shift_files/tmp.txt

    # run electio with weights
    # in this mode, supplied weights are used to compute fitness during selection
    # binary mode
    #electio \
    #--ga_binary \
    #--ga_alpha 1.0 \
    #--ga_populations 20 \
    #--ga_cycles 100 \
    #--ga_seed ${RANDOM} \
    #--ga_monitor \
    #--input rmsd_files/${rna}.txt \
    #--ga_weights electio_selections/weights_vector_${rna}.txt \
    #--output electio_selections/${rna}_binary.txt \
    #electio_selections/observed_electio_vector_${rna}.txt \
    #electio_selections/predicted_matrix_${rna}.txt
    
    #electio \
    #--ga_alpha 1.0 \
    #--ga_populations 20 \
    #--ga_cycles 100 \
    #--ga_seed ${RANDOM} \
    #--ga_monitor \
    #--input rmsd_files/${rna}.txt \
    #--ga_weights electio_selections/weights_vector_${rna}.txt \
    #--output electio_selections/${rna}_real.txt \
    #--ga_mask electio_selections/${rna}_binary.txt \
    #--ga_mcolumn 4 \
    #electio_selections/observed_electio_vector_${rna}.txt \
    #electio_selections/predicted_matrix_${rna}.txt

    electio \
    --ga_alpha 1000.0 \
    --ga_populations 100 \
    --ga_cycles 5000 \
    --ga_seed ${RANDOM} \
    --ga_monitor \
    --input rmsd_files/${rna}.txt \
    --ga_weights electio_selections/weights_vector_${rna}.txt \
    --output electio_selections/${rna}_real.txt \
    electio_selections/observed_electio_vector_${rna}.txt \
    electio_selections/predicted_matrix_${rna}.txt

done
exit
