#!/bin/bash -l
#SBATCH -p debug
#SBATCH -t 00:30:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -C knl
##SBATCH -S 4

module load python
export HDF5_USE_FILE_LOCKING=FALSE


python3 ${PATH_TO_SCRIPT}/pseudobands_opt.py --fname_in WFN.h5 --fname_in_q WFNq.h5 --fname_out WFN_SPB.h5 --fname_out_q WFN_SPB_q.h5 --nv 1 --nc 1 --nslice_v 5 --nslice_c 175 --max_freq 0.0 --nspbps_v 5 --nspbps_c 5 --verbosity 2 


# you can also just run this on the login node, since its a serial script
# if you do, remove the HDF5_USE_FILE_LOCKING variable
