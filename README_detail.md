# STOCHASTIC PSEUDOBANDS

*Altman A. R., Kundu S., da Jornada F. H., "Mixed Stochastic-Deterministic Approach For Many-Body Perturbation Theory Calculations" (2023).*

## **Purpose**

This script massively compresses the mean-field wavefunciton to speed up GW calculations based on sum-over-bands formalisms. It reduced the scaling of the GW code, speedup up the calculation of the dieleectric matrix and self-energy often by 1-3 orders of magnitude or more, and eliminates traditional band-truncation convergence parameters. 

## **How It Works**

Given a set of input parameters (see below) and input **WFN.h5** file in BerkeleyGW format (http://manual.berkeleygw.org/3.0/wfn_h5_spec/) containing many bands, *pseudobands.py* constructs a set of energy subspaces/slices for conduction bands and optionally valence bands, and creates stochastic linear combinations of the mean-field bands in each energy slice. The stochastic pseudobands in a given slice are all degenerate, with energy equal to the mean energy of the slice. There are also a number of **protected bands**, which are simply copied from the input file, and not included in the slicing process. These should at least contain the bands for which the self-energy is computed. The output is a **WFN_SPB.h5** file which has the same format as **WFN.h5** but contains the stochastic pseudobands, and a **phases.h5** file that contains the random coefficients used to construct the valence SPBs. Different random coefficients are used for different k-points. The **q**-shifted wavefunction **WFNq.h5** used for computing the dielectric function as **q --> 0** has to be treated simultaneously, so **WFN.h5 and WFNq.h5 must both be provided as input.**

### **Construction Of Slices**
As described in the paper, there are two ways to determine how to partition the mean-field bands into subspaces/slices in which stochastic pseudobands are constructed. In this implementation, one specifies the total number of slices with the input parameters **N_S_val** and/or **N_S_cond**, and the energy fraction ***F*** is deduced (see paper for explanation). To do so, we enforce an exponential ansatz on the energy spanned by the bands in each subspace, i.e. [energy spanned subspace S_j] = α * exp(j * β) for some fitting parameters α, β. By doing this, we can deduce ***F*** as being roughly proportional to β.

## <span style="color:red">**Important Warnings**</span>
- Do not use valence pseudobands when computing the self-energy, as this leads to errors on the **1 eV** scale! Only the dielectric computation can use valence pseudobands!
- Explicit frequency-dependent dielectric calculations require extra input parameters! See optional flags **max_freq** and **uniform_width** corresponding to this case.


## **Input Files**
- **WFN.h5**: WFN file in *h5* format. Should contain many unoccupied bands from either ParaBands or nscf calculation from a DFT code that has been run through pw2bgw.x and then wfn2hdf.x
- **WFNq.h5**: q-shifted WFN file in *h5* format. Should contain all occupied and at least one unoccupied band.

## **Output Files**
- **WFN_SPB.h5**: WFN file containing stochastic pseudobands
- **WFN_SPB_q.h5**: WFNq file containing stochastic pseudobands
- **WFN_SPB.log**: Log file containing all input parameters, filenames, and errors, if any.
- **phases.h5**: h5 file containing the random coefficients used to construct the *valence* pseudobands. Used for logging and constructing SPBs for **WFNq** if symmetries are not used. 

## **Input paramaters**
***Required***
- **fname_in (str)**: Input WFN.h5, in BerkeleyGW HDF5 format
- **fname_out (str)**: Output WFN_SPB.h5 with pseudobands, in HDF5 format
- **fname_in_q (str)**: Input WFNq.h5, in BerkeleyGW HDF5 format
- **fname_out_q (str)**: Output WFN_SPB_q.h5 with pseudobands, in HDF5 format

***(Strongly) Recommended***
- **N_P_val (int >= -1)**: Number of protected valence bands counting from VBM. If N_P_val == -1 then all valence bands are copied (i.e. no valence SPBs), which is preferable if there are less than ~100 valence states. Default == -1.
- **N_P_cond (int >= -1)**: Number of protected conduction bands counting from CBM. If N_P_cond == -1 then all conduction bands are copied (i.e. no conduction SPBs). Default == 100
- **N_S_val (int >= 0)**:Number of subspaces spanning the total energy range of the valence bands. Default == 10
- **N_S_cond (int >= 0)**: Number of subspaces spanning the total energy range of the conduction bands. Default == 100
- **N_xi_val (int >= 2)**: Number of stochastic pseudobands constructed per valence slice. **You must set this value to be at least 2!!** Default == 2.
- **N_xi_cond (int >= 2)**: Number of stochastic pseudobands constructed per conduction slice. **You must set this value to be at least 2!!** Default == 2.

***Optional***
- **NNS (int=0,1)**: If using a separate WFNq.h5 for nonuniform-neck-subsampling (NNS, see http://manual.berkeleygw.org/3.0/NNS/), set **NNS = 1**. The NNS WFNq fname_in/fname_out flags will be ignored without this flag. Default == 0.
- **fname_in_NNS (str)**: Input NNS WFNq.h5, in HDF5 format. Default == None
- **fname_out_NNS (str)**: Output NNS WFNq.h5 with pseudobands, in HDF5 format. Default == None.
- **max_freq (float)**: Maximum energy (Ry) for usage of uniform_freq, for full frequency calculations. Should be greater than the largest energy in the full-frequency calculation in epsilon. If **uniform_width** is not set then this flag does nothing. Default == 1.0.
- **uniform_width (float)**: energy width of constant slice slices. Required for use of **max_frequency**. Should be half of the spacing of the freqency grid. Default == None.
- **fname_phases (str)**: Phases.h5 file, containing random coefficients for SPB construction for valence states. Should be consistent with all other parameters. Intended for use with WFNq calculation. Default == None.
- **verbosity (int 0-3)**: Verbosity of output, written to the log file. Default == 2.


### **Example Run Command**
```bash
python pseudobands.py --fname_in WFN.h5 --fname_in_q WFNq.h5 --fname_out WFN_SPB.h5 --fname_out_q WFN_SPB_q.h5 --N_P_val 10 --N_P_cond 10 --N_S_val 10 --N_S_cond 150 --N_xi_val 2 --N_xi_cond 2
```
