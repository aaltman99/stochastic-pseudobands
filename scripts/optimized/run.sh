#!/bin/bash -l
#SBATCH -p debug
#SBATCH -t 00:30:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -C knl
##SBATCH -S 4

module load python
export HDF5_USE_FILE_LOCKING=FALSE

python3 /pseudobands_opt.py --fname_in WFN.h5 --fname_in_q WFNq.h5 --fname_out WFN_SPB.h5 --fname_out_q WFN_SPB_q.h5 --nv 5 --nc 100 --nslice_v 5 --uniform_width 0.0 --nslice_c 20 --max_freq 0.0 --nspbps_v 1 --nspbps_c 1 --verbosity 2 &
wait
