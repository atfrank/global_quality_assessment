#!/bin/bash

source activate mdanalysis
export LARMORPYTHON=/Users/atfrank/GitSoftware/RingCurrentModels/python
export PYTHONPATH=$LARMORPYTHON:$PYTHONPATH
cd /Users/atfrank/GitSoftware/global_quality_assessment/spath

python ${LARMORPYTHON}/larmorExtractor.py --output=features_1R2P_2LPS.txt --id=1R2P_2LPS --nproc=8 --chemshift=scratch/measured_shifts_2LPS.dat --verbose coors/2lps_ave_min.psf scratch/1r2p_2lps_trajectory.dcd
python ${LARMORPYTHON}/larmorExtractor.py --output=features_2FRL_2M22.txt --id=2FRL_2M22 --nproc=8 --chemshift=scratch/measured_shifts_2M22.dat --verbose coors/2m22_ave_min.psf scratch/2frl_2m22_trajectory.dcd
python ${LARMORPYTHON}/larmorExtractor.py --output=features_2H2X_2M21.txt --id=2H2X_2M21 --nproc=8 --chemshift=scratch/measured_shifts_2M21.dat --verbose coors/2m21_ave_min.psf scratch/2h2x_2m21_trajectory.dcd
python ${LARMORPYTHON}/larmorExtractor.py --output=features_2KFC_2L1V.txt --id=2KFC_2L1V --nproc=8 --chemshift=scratch/measured_shifts_2L1V.dat --verbose coors/2l1v_ave_min.psf scratch/2kfc_2l1v_trajectory.dcd
python ${LARMORPYTHON}/larmorExtractor.py --output=features_2N6Q_5KMZ.txt --id=2N6Q_5KMZ --nproc=8 --chemshift=scratch/measured_shifts_5KMZ.dat --verbose coors/5kmz_ave_min.psf scratch/2n6q_5kmz_trajectory.dcd

source deactivate

