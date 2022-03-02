#!/bin/bash -l
#SBATCH -p debug
#SBATCH -t 00:15:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -C knl

module load python3
export HDF5_USE_FILE_LOCKING=FALSE

# --uniform_width 0.00734 is in Ry. 0.00734 Ry == 100 meV

python3 $PATH-TO-SCRIPT/pseudobands.py --fname_in WFN_in.h5 --fname_in_q WFNq.h5 --fname_out WFN_SPB.h5 --fname_out_q WFN_SPB_q.h5 --nv 50 --nc 50 --efrac_v 0.01 --uniform_width 0.00734 --efrac_c 0.02 --max_freq 1.0 --nspbps_v 1 --nspbps_c 1 --verbosity 1 &
wait
