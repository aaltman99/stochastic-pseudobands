# STOCHASTIC PSEUDOBANDS (SPBs) 
## (**optimized version**)

***Aaron R. Altman (02/2022)***

**Originally written by Felipe H. da Jornada (2015)**

old/pseudobands_FHJ.py is Felipe's original code

**pseudobands_opt.py is the current working version.** 

## **Purpose**

The **WFN** file output by these scripts *exponentially* reduces the number of bands required to perform the GW calculation. There is currently **no mean-field speedup**, but eventually a spectral slicing algorithm will be implemented in ParaBands that can cirvumvent the current requirement to generate deterministic unoccupied states. **The GW calculation scales quasi-quadratically rather than quarticly with this approach.**


## **How It Works**

Given a set of input parameters (explained below) and input **WFN.h5** file containing many bands, *pseudobands.py* constructs a set of energy slices for both valence and conduction bands, and creates stochastic linear combinations of the bands in each energy slice. The stochastic pseudobands in a given slice are all degenerate, with energy equal to the mean energy of the slice. There are also a number of **protected bands**, which are simply copied from the input file, and not included in the slicing process. These should at least contain the bands for which the self energy is computed. The output is a **WFN_SPB.h5** file which has the same format as **WFN.h5** but contains the stochastic pseudobands, and a **phases.h5** file that contains the random coefficients used to construct the valence SPBs. Different random coefficients are used for different k-points. **WFN.h5 and WFNq.h5 are treated simultaneously, and must both be provided as input.**

### **Construction Of Slices**
Slices are constructed according an optimization procedure designed to minimize the expected error of the Green's function given the input parameters. The slices are given by an exponential distribution, i.e. the energy-width of each slice is given by:
        **w_j = alpha * exp(beta * j)**
for some parameters **alpha, beta** that are obtained from the optimization. The parameters going into the optiization are **E0, Emax, nslice,** and **nspbps,** where **E0** and **Emax** are the energy of the first unprotected band and the energy of the last band in the input **WFN.h5** file, respectively.

Note that the **nslice** parameters are directly related to the computational cost of the GW calculation, with epsilon scaling as **nslice_v * nspbps_v * nslice_c * nspbps_c**. Therefore, these parameters should be determined by the resources you are willing to expend on the GW calculation.

### **Convergence Testing**
Convergence testing is very straighforward in this scheme. Once **WFN_SPB.h5** and **WFN_SPB_q.h5** are obtained, the script *pseudobands-convergence.py* can be run on each WFN_SPB separately to remove one SPB per slice in the WFN, thereby reducing the total number of bands, and allowing convergence testing with respect to **nspbps**. The script can be run recursively if you want to remove multiple bands per slice. *Note*: this script creates a separate file and leaves the original untouched, and can handle both the fortran and python versions of pseudobands. The user can also specify whether to do this procedure for only valence pseudobands, conduction pseudobands, or both. See *run-conv.sh* for an example run script.

**Note:** you must have at least 2 bands per slice for this procedure to make sense. Since converged results only appear after 2 bands per slice, it is recommended to use at least **nspbps = 3** if you plan on doing convergence testing.

If you wish to perform convergence testing for epsilon, then **WFN_SPB.h5** and **WFN_SPB_q.h5** must be treated identially. If you wish to test sigma, then only run on **WFN_SPB.h5** for conduction pseudobands only (see below).

## ***Warning*: Important Use Cases**
The approximations made by this algorithm have vanishing stochastic error only when the computed quantity has a 1/energy term. This means epsilon can handle both conduction and valence pseudobands but <span style="color:red">**sigma should not be used with valence pseudobands**</span> due to the exchange term that weights all valence bands equally!! 


## **Input Files**
- **WFN.h5**: WFN file in *h5* format. Should contain many unoccupied bands from either ParaBands or nscf calculation from a DFT code that has been run through pw2bgw.x and then wfn2hdf.x
- **WFNq.h5**: q-shifted WFN file in *h5* format. Should contain all occupied and at least one unoccupied band.

## **Output Files**
- **WFN_SPB.h5**: WFN file containing stochastic pseudobands
- **WFN_SPB.log**: Log file containing all input parameters, filenames, and errors, if any.
- **phases.h5**: h5 file containing the random coefficients used to construct the *valence* SPBs. Used for constructing SPBs for **WFNq**. 

## **Input paramaters**
***Required***
- **fname_in (str)**: Input WFN.h5, in HDF5 format
- **fname_out (str)**: Output WFN.h5 with pseudobands, in HDF5 format
- **fname_in_q (str)**: Input WFNq.h5, in HDF5 format
- **fname_out_q (str)**: Output WFNq.h5 with pseudobands, in HDF5 format

***Recommended***
- **nv (int >= -1)**: Number of protected valence bands counting from VBM. If nv == -1 then all valence bands are copied (i.e. no valence SPBs), which is preferable if there are less than ~100 valence states. Default == -1.
- **nc (int >= -1)**: Number of protected conduction bands counting from CBM. If nc == -1 then all conduction bands are copied (i.e. no conduction SPBs). Default == 100
- **nslice_v (int >= 0)**:Number of subspaces spanning the total energy range of the valence bands. Default == 10
- **nslice_c (int >= 0)**: Number of subspaces spanning the total energy range of the conduction bands. Default == 100
- **nspbps_v (int >= 2)**: Number of stochastic pseudobands constructed per valence slice. Typically set this higher than nspbps_c, or do not use valence SPBs for sigma. **You must set this value to be at least 2!!** Default == 2.
- **nspbps_c (int >= 2)**: Number of stochastic pseudobands constructed per conduction slice. **You must set this value to be at least 2!!** Default == 2.

***Optional***
- **NNS (int=0,1)**: If using a separate WFNq.h5 for NNS, set **NNS = 1**. The NNS WFNq fname_in/fname_out flags will be ignored without this flag. Default == 0.
- **fname_in_NNS (str)**: Input NNS WFNq.h5, in HDF5 format. Default == None
- **fname_out_NNS (str)**: Output NNS WFNq.h5 with pseudobands, in HDF5 format. Default == None.
- **max_freq (float)** (optional): Maximum energy (Ry) for usage of uniform_freq, for full frequency calculations. Should be greater than the largest energy in the full-frequency calculation in epsilon. If **uniform_width** is not set then this flag does nothing. Default == 1.0.
- **uniform_width (float)** (optional): energy width of constant slice slices. Required for use of **max_frequency**. Should be half of the spacing of the freqency grid. Default == None.
- **copydirectly (bool)** (optional): Direct copying for protected bands. If False, then copying is done in chunks to limit memory usage. Set to False is you have a large number of protected bands. Default == True.
- **fname_phases (str)** (optional): Phases.h5 file, containing random coefficients for SPB construction for valence states. Should be consistent with all other parameters. Intended for use with WFNq calculation. Default == None.
- **verbosity (int 0-3)** (optional): Verbosity of output, written to the log file. Default == 0.


### **Example Run Command**
python3 ${PATH-TO-SCRIPT}/pseudobands.py --fname_in WFN_in.h5 --fname_in_q WFNq.h5 --fname_out WFN_SPB.h5 --fname_out_q WFN_SPB_q.h5 --nv 50 --nc 50 --efrac_v 0.01 --efrac_c 0.02 --uniform_width 0.0035 --max_freq 1.0 --nspbps_v 2 --nspbps_c 1

See also *run.sh* for an example job submission script.

### **Best Practices**
- If you run several calculations you will need to keep track of the **WFN_SPB.h5** output files. The SPB parameters are NOT logged in the WFN_SPB.h5 file. This is because it must maintain the standard **WFN** file format. **Set the parameters used as part of the WFN_SPB.h5 filename**.
- If running multiple calculations, use the workflow scripts in the *workflow* folder to stay organized and save time.
