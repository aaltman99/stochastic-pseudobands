#!/bin/bash -l

module load python


python3 ${PATH_TO_SCRIPT}/pseudobands-convergence.py --fname_in WFN_SPB.h5 --do_val false --do_cond true
