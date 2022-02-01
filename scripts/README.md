# STOCHASTIC PSEUDOBANDS (SPBs)

***Aaron R. Altman (01/2022)***

**Originally written by Felipe H. da Jornada (2015)**

old/pseudobands_FHJ.py is Felipe's original code

**pseudobands.py is the current working version.** 

## **Purpose**

The **WFN** file output by these scripts *exponentially* reduces the number of bands required to perform the GW calculation. There is currently **no mean-field speedup**, but eventually a spectral slicing algorithm will be implemented in ParaBands that can cirvumvent the current requirement to generate deterministic unoccupied states. **The GW calculation becomes almost independent of system size with this approach.**


## **How It Works**

Given a set of input parameters (explained below) and input **WFN.h5** file containing many bands, *pseudobands.py* constructs a set of energy slices for both valence and conduction bands, and creates stochastic linear combinations of the bands in each energy slice. The stochastic pseudobands in a given slice are all degenerate, with energy equal to the mean energy of the slice. There are also a number of **protected bands**, which are simply copied from the input file, and not included in the slicing process. These should at least contain the bands for which the self energy is computed. The output is a **WFN_SPB.h5** file which has the same format as **WFN.h5** but contains the stochastic pseudobands, and a **phases.h5** file that contains the random coefficients used to construct the valence SPBs. Different random coefficients are used for different k-points. **WFN and WFNq are treated simultaneously, and must both be provided as input.**

### **Construction Of Slices**
Slices are constructed according to the input parameters **efrac_v**. Each slice is given as the bands that fit in the energy interval 
		**[E_first, E_last] = [E, E(1 + efrac)]**, 
where 
		**E_first**(interval [*i*]) = **E_last**(interval [*i - 1*]), 
and **E_first**(interval [*1*]) is set by **nv**/**nc**.

## **Input Files**
- **WFN.h5**: WFN file in *h5* format. Should contain many unoccupied bands from either ParaBands or nscf calculation from a DFT code that has been run through pw2bgw.x and then wfn2hdf.x
- **WFNq.h5**: q-shifted WFN file in *h5* format. Should contain all occupied and at least one unoccupied band.

## **Output Files**
- **WFN_SPB.h5**: WFN file containing stochastic pseudobands
- **WFN_SPB.log**: Log file containing all input parameters, filenames, and errors, if any.
- **phases.h5**: h5 file containing the random coefficients used to construct the *valence* SPBs. Used for constructing SPBs in for **WFNq**. 

## **Input paramaters**
***Required***
- **fname_in (str)**: Input WFN.h5, in HDF5 format
- **fname_out (str)**: Output WFN.h5 with pseudobands, in HDF5 format
- **fname_in_q (str)**: Input WFNq.h5, in HDF5 format
- **fname_out_q (str)**: Output WFNq.h5 with pseudobands, in HDF5 format

***Recommended***
- **nv (int >= -1)**: Number of protected valence bands counting from VBM. If nv == -1 then all valence bands are copied (i.e. no valence SPBs), which is preferable if there are less than ~100 electrons. Default == -1.
- **nc (int >= 0)**: Number of protected conduction bands counting from CBM. Must be nc >= 0. Default == 100
- **efrac_v (float < 1)**: Accumulation window for valence slices, as a fraction of the energy of the band in each subspace. Default == 0.01.
- **efrac_c (float < 1)**: Accumulation window for conduction slices, as a fraction of the energy of the band in each subspace. Default == 0.01.
- **nspbps_v (int >= 1)**: Number of stochastic pseudobands constructed per valence slice. Typically set this higher than nspbps_c, or do not use valence SPBs for sigma. Default == 1.
- **nspbps_c (int >= 1)**: Number of stochastic pseudobands constructed per conduction slice. Default == 1.

***Optional***
- **efrac_c_fine (float < 1)** (optional): Accumulation window for conduction slices with energies <= max_freq, as a fraction of the energy of the band in each subspace. Intended for use with full-frequency in epsilon, otherwise leave unset. Default: **efrac_c**.
- **max_freq (float)** (optional): Maximum energy (Ry) for usage of efrac_c_fine. Should be greater than the largest energy in the full-frequency calculation in epsilon. If **efrac_c_fine** is not set then this flag does nothing. Default == 1.0.
- **copydirectly (bool)** (optional): Direct copying for protected bands. If False, then copying is done in chunks to limit memory usage. Set to False is you have a large number of protected bands. Default == True.
- **fname_phases (str)** (optional): Phases.h5 file, containing random coefficients for SPB construction for valence states. Should be consistent with all other parameters. Intended for use with WFNq calculation. Default == None.
- **verbosity (int 0-3)** (optional): Verbosity of output, written to the log file. Default == 0.


### **Example Run Command**
python3 $PATH-TO-SCRIPT/pseudobands.py --fname_in WFN_in.h5 --fname_in_q WFNq.h5 --fname_out WFN_SPB.h5 --fname_out_q WFN_SPB_q.h5 --nv 50 --nc 50 --efrac_v 0.01 --efrac_c_fine 0.01 --efrac_c 0.02 --max_freq 1.0 --nspbps_v 2 --nspbps_c 1

### **Best Practices**
- If you run several calculations you will need to keep track of the **WFN_SPB.h5** output files. The SPB parameters are NOT logged in the WFN_SPB.h5 file. This is because it must maintain the standard **WFN** file format. **Set the parameters used as part of the WFN_SPB.h5 filename**.
- Setting the output filename is not part of the code because you may want to run pseudobands multiple times to estimate the variance on the GW quantities. In this case you should append a run number to your output filename, like **WFN_SPB...(parameters)..._run1.h5**.
- If running multiple calculations, use the workflow scripts in the *workflow* folder to stay organized and save time.
