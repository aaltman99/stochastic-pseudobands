#!/usr/bin/env python

# Pseudobands!

### ARA: mod for nanoparticles, val and cond SPBs, phases.h5 (01/22)
### NPs: need 2 efracs (fine + coarse) for cond bands to properly capture frequency dependence in epsilon
### val SPB: need higher nspbps, otherwise error from exchange self energy is large
### however, need to ensure there are still fewer total val SPBs, otherwise no savings

### TODO: automate selection of these parameters; optimize parameters 
### TODO: fix copydirectly=False case in pseudobands


### current scheme only works for gapped materials (?)
### current scheme only works for not-too-large spin-orbit splitting




from __future__ import print_function
import h5py
import numpy as np
import scipy as sp
from scipy import signal
import logging 
from inspect import currentframe, getframeinfo
import time
import sys


# valence bands aggregate in energy steps, almost in a piecewise-constant fashion. Slices should respect the positions of these steps
def get_steps(el):
    
    grad = np.gradient(el)
    kernel_size = [len(el)//300 if len(el)>300 else 2][0]
    kernel = np.ones(kernel_size)
    grad_convolved = np.convolve(grad, kernel, mode='same')
    
    var = np.var(grad_convolved)
    skew = sp.stats.skew(grad_convolved)
    
    if skew <= var: # no large steps found
        return [-1, len(el)]
    
    thresh = 0.1 
    
    edges = signal.find_peaks(grad_convolved, height=thresh*np.max(grad_convolved), distance=3*kernel_size)
    
    return list(edges[0])


def construct_blocks(el, n_copy, efrac, efrac_fine, max_freq):
    
    if efrac_fine is None:
        efrac_fine = efrac
        
    blocks = []
    nb_out = n_copy
    first_idx = n_copy
    
    edges = get_steps(el)
    
    while True:
        first_en = el[first_idx]
        assert first_en > 0
        
        if first_en <= max_freq:
            delta_en = first_en * efrac_fine
        else:
            delta_en = first_en * efrac
            
        last_en = first_en + delta_en

        try:
            last_idx = list(np.where(el > last_en))[0][0]
            
            if last_idx >= edges[0]:
                last_idx = edges[0]
                edges.pop(0)
            
            blocks.append([first_idx, last_idx - 1])

            nb_out += 1
            first_idx = last_idx
            
        except:
            last_idx = len(el)
            blocks.append([first_idx, last_idx - 1])
            nb_out += 1
            break
            
    return np.asarray(blocks)


# fix signs/relative alignment
def fix_blocks(blocks, vc, ifmax, nb_orig):
    if vc == 'v':
        np.flip(blocks)
        blocks *= -1
        blocks += ifmax-1
        assert blocks[-1][-1] == 0
    
    else:
        blocks += ifmax
        assert blocks[-1][-1] == nb_orig-1
        

# bunch of sanity checks, block construction for WFN and WFNq simultaneously 
def check_and_block(fname_in = None, fname_in_q = None, nv=-1, nc=1, efrac_v=0.01, efrac_c=0.02, efrac_c_fine=None, max_freq=1.0, nspbps_v=1, nspbps_c=1, verbosity=0, **kwargs):
    
    nv = int(nv)
    nc = int(nc)
    
    if efrac_c_fine is None:
        efrac_c_fine = efrac_c
    
    assert nv >= -1
    assert nc >= 0
    assert efrac_v > 0
    assert efrac_c_fine > 0
    assert efrac_c > 0
    
    
    f_in = h5py.File(fname_in, 'r')
    f_in_q = h5py.File(fname_in_q, 'r')
    
    en_orig = f_in['mf_header/kpoints/el'][()]
    nb_orig = f_in['mf_header/kpoints/mnband'][()]
    ifmax = np.asarray(f_in['mf_header/kpoints/ifmax'][()], dtype=int)[0]
    
    en_orig_q = f_in_q['mf_header/kpoints/el'][()]
    nb_orig_q = f_in_q['mf_header/kpoints/mnband'][()]
    ifmax_q = np.asarray(f_in_q['mf_header/kpoints/ifmax'][()], dtype=int)[0]
    
    try:
        assert np.array_equal(ifmax/ifmax[0], np.ones(ifmax.shape[0]))
    except AssertionError as err:
        logger.error('ifmax is not the same across all kpoints. This script does not currently support non-uniform band occupations.')
        raise err
        
    try:
        assert np.array_equal(ifmax_q/ifmax_q[0], np.ones(ifmax_q.shape[0]))
    except AssertionError as err:
        logger.error('ifmax_q is not the same across all kpoints. This script does not currently support non-uniform band occupations.')
        raise err
    
    try:
        assert np.array_equal(ifmax, ifmax_q)
    except AssertionError as err:
        logger.error('ifmax is different from ifmax_q. This is nonphysical, make sure your WFN and WFNq files are correct.')
        raise err
    
    ifmax = ifmax[0] # want a scalar, TODO: add array support?
    
    logger.info(f'ifmax = {ifmax}')
    
    if nv == -1:
        nv = ifmax
    
    lastcopy_en_v = en_orig[:,:,ifmax-1-nv]
    lastcopy_en_c = en_orig[:,:,ifmax-1+nc]
          
    
    if any(np.in1d(lastcopy_en_v, en_orig[...,ifmax-nv:ifmax])):
        logger.error('Chosen nv cuts degenerate bands, use different nv. Try degeneracy_check.x')
        raise Exception('Chosen nv cuts degenerate bands, use different nv. Try degeneracy_check.x')
    if any(np.in1d(lastcopy_en_c, en_orig[...,ifmax:ifmax+nc-2])):
        logger.error('Chosen nc cuts degenerate bands, use different nc. Try degeneracy_check.x')
        raise Exception('Chosen nc cuts degenerate bands, use different nc. Try degeneracy_check.x')
    
    
    VBMCBM = np.concatenate((en_orig[...,ifmax-1:ifmax+1], en_orig_q[...,ifmax-1:ifmax+1]), axis=1)
    
    
    # ifmax is 1-indexed, en_orig is 0-indexed
    # also fortran order...
    fermi = np.mean([np.max(VBMCBM[...,0]), np.min(VBMCBM[...,1])]) 
    logger.info(f'fermi = {fermi} Ry')

    
    assert nv <= ifmax
    assert nc <= nb_orig - ifmax
    
    en_orig -= fermi
    en_orig_q -= fermi
   
    # average over kpoints
    el_all_v = -np.flip(np.mean(np.concatenate((en_orig[...,:ifmax], en_orig_q[...,:ifmax]), axis=1)[0], axis=0))

    el_c = np.mean(en_orig[0,:,ifmax:None], axis=0)    
    # hopefully there are no states at 0 energy... add assert statement

    
    blocks_v = construct_blocks(el_all_v, nv, efrac_v, None, max_freq)
    blocks_c = construct_blocks(el_c, nc, efrac_c, efrac_c_fine, max_freq)
    
    fix_blocks(blocks_v, 'v', ifmax, nb_orig)
    fix_blocks(blocks_c, 'c', ifmax, nb_orig)
    
    
    
    if verbosity > 0:
        logger.info(f'valence slices: {blocks_v}')
        logger.info(f'conduction slices: {blocks_c}')
        
    f_in.close()
    f_in_q.close()
    
    nslices_c = len(blocks_c)
    nslices_v = len(blocks_v)
    nspb_c = nslices_c * nspbps_c
    nspb_v = nslices_v * nspbps_v
    nb_out_c = nc + nspb_c
    nb_out_v = nv + nspb_v
    
    try:
        assert nb_out_v <= ifmax
    except AssertionError as err:
        logger.error('More total valence states (copied + SPBs) than original valence bands ==> no computational savings in GW steps. Choose smaller nv, nspbps_v, or larger efrac_v')
        raise err
    
    assert len(blocks_c) == nslices_c
    assert len(blocks_v) == nslices_v
    
    logger.info(f'''nslices_c = {nslices_c}\n nslices_v = {nslices_v}\n nspb_c_total = {nspb_c}\n nspb_v_total = {nspb_v}\n nb_out_c_total = {nb_out_c}\n nb_out_v_total = {nb_out_v}\n''') 
    
    
    return blocks_v, blocks_c, ifmax
    
    
# nv and nc bands are copied. SPBS are constructed outside this range
def pseudoband(qshift, blocks_v, blocks_c, ifmax, fname_in = None, fname_out = None, fname_in_q = None, fname_out_q = None, nv=-1, nc=1, nspbps_v=1, nspbps_c=1, single_band=False, copydirectly=True, verbosity=0, fname_phases=None, **kwargs):
    
    start = time.time()
    
    nv = int(nv)
    nc = int(nc)
    
    assert nv >= -1
    assert nc >= 0
    
    nslices_c = len(blocks_c)
    nslices_v = len(blocks_v)
    nspb_c = nslices_c * nspbps_c
    nspb_v = nslices_v * nspbps_v
    nb_out_c = nc + nspb_c
    nb_out_v = nv + nspb_v
    
    
    if single_band:
        nspbps_v = 1
        nspbps_c = 1
        
    if not qshift:
        f_in = h5py.File(fname_in, 'r')
        f_out = h5py.File(fname_out, 'w')
        logger.info(f'fname_in = {fname_in}\n fname_out = {fname_out}\n ')
    else:
        f_in = h5py.File(fname_in_q, 'r')
        f_out = h5py.File(fname_out_q, 'w')
        logger.info(f'fname_in = {fname_in_q}\n fname_out = {fname_out_q}\n ')
    
    nk = f_in['mf_header/kpoints/nrk'][()]
    ngk = f_in['mf_header/kpoints/ngk'][()]
    cum_ngk = np.insert(np.cumsum(ngk), 0, 0)
    
    if fname_phases is None:
        phases_file = h5py.File('phases'+fname_out[0:-3]+'.h5', 'w')
    else:
        phases_file = h5py.File(fname_phases, 'r')
    
   
    f_out.copy(f_in['mf_header'], 'mf_header')
    f_out.create_group('wfns')
    f_out.copy(f_in['wfns/gvecs'], 'wfns/gvecs')
    f_out['mf_header/kpoints/mnband'][()] = nb_out_v + nb_out_c
    f_out['mf_header/kpoints/ifmax'][()] = nb_out_v
    
    if fname_phases is None:
        phases_file.copy(f_in['mf_header'], 'mf_header')
        phases_file.create_group('phases')

        
    def resize(file, name):
        file.move(name, name + '_orig_pb')
        shape = list(file[name + '_orig_pb'].shape)
        shape[-1] = nb_out_v + nb_out_c
        file.create_dataset(name, shape, dtype='d')
        file[name][:, :, :nb_out_v+nb_out_c] = file[name + '_orig_pb'][:, :, :nb_out_v+nb_out_c]
        del file[name + '_orig_pb']

            
    resize(f_out, 'mf_header/kpoints/occ')
    resize(f_out, 'mf_header/kpoints/el')
    
    
    logger.info('Copying {} protected bands'.format(nv+nc))
    shape = list(f_in['wfns/coeffs'].shape)
    shape[0] = nb_out_v + nb_out_c
    f_out.create_dataset('wfns/coeffs', shape, 'd')
    
    if fname_phases is None:
        shape_phases = (ifmax-nv, nk, nspbps_v, 2,)
        phases_file.create_dataset('phases/coeffs', shape_phases, 'd')
   
    
    if copydirectly:
        # FHJ: this can take quite a while for some systems..
        logger.warning('copydirectly=True. Copying all protected bands at once, can be slow if copying >1000 bands')
        f_out['wfns/coeffs'][nb_out_v-nv:nb_out_v+nc, :, :] = f_in['wfns/coeffs'][ifmax-nv:ifmax+nc, :, :]
        f_out['mf_header/kpoints/el'][:,:,nb_out_v-nv:nb_out_v+nc] = f_in['mf_header/kpoints/el'][:,:,ifmax-nv:ifmax+nc]
        
        
    else: ### FIXME: indices not correct (if copydirectly indices are correct though); phases_file not included yet
        nbs_block = 1000
        ibs_start = range(ifmax-nv, ifmax + nc -1, nbs_block)
        if Bar is not None:
            bar = Bar('Copying protected bands', max=len(ibs_start), bar_prefix=' [', bar_suffix='] ',
                      fill='#', suffix='%(percent)d%% - Remaining: %(eta_td)s')
        for ib in ibs_start:
            if Bar is not None:
                bar.next()
            ib1, ib2 = ib, min(ib + nbs_block, ifmax+nv+nc-1)
            tmp = f_in['wfns/coeffs'][ib1:ib2, :, :, :]
            try:
                f_out['wfns/coeffs'][ib1:ib2, :, :] = tmp
            except Exception as err:
                frameinfo = getframeinfo(currentframe())
                logger.error(frameinfo.lineno)
                logger.error(err)
        if Bar is not None:
            bar.finish()
    
    
    logger.info('Creating {} valence pseudobands'.format(nspb_v))
    ib = nb_out_v-nv-1

    for b in blocks_v:
        if verbosity > 1:
            logger.info(f'slice_index, slice: {ib, b}')

        if single_band:
            band_avg = b[0] + (b[1] - b[0]) // 2
            f_out['wfns/coeffs'][ib, :, :] = (f_in['wfns/coeffs'][band_avg, :, :]
                                              * np.sqrt(float(b[1] - b[0] + 1)))
            f_out['mf_header/kpoints/el'][:, :, ib] = f_in['mf_header/kpoints/el'][:, :, band_avg].mean(axis=-1)
            ib -= 1
        elif nspbps_v == 1:
            if b[0] == 1:
                f_out['wfns/coeffs'][ib, :, :] = f_in['wfns/coeffs'][0, :, :].sum(axis=0)
                f_out['mf_header/kpoints/el'][:, :, ib] = f_in['mf_header/kpoints/el'][:, :, 0].mean(axis=-1)
            elif b[0] == 0:
                continue
            else:
                f_out['wfns/coeffs'][ib, :, :] = f_in['wfns/coeffs'][b[1]:b[0]+1, :, :].sum(axis=0)
                f_out['mf_header/kpoints/el'][:, :, ib] = f_in['mf_header/kpoints/el'][:, :, b[1]:b[0]+1].mean(axis=-1)
                
            ib -= 1
            
        else: # nspbps > 1
            if b[0] == 0:
                continue
            
            elif b[0] == b[1] and verbosity > 2:
                logger.warning(f'Current slice {b} has the same start and end bands. Using nspbps > 1 will introduce unnecessary error. Hack the code to automatically deal with this, or ignore it, but using valence SPBs in sigma may be an extremely poor approximation due to the exchange term.')
                # coeffs = f_in['wfns/coeffs'][b[0], :, :, :]
                # el = f_in['mf_header/kpoints/el'][:, :, b[0]]
                # ib -= 1
                
                ### WARNING: uncommenting ^^^ does not work, since the h5 dataset is a fixed size.
                ### (purpose is to just copy the band instead of putting extra SPBs)
                
            else:
                coeffs = f_in['wfns/coeffs'][b[1]:b[0]+1, :, :, :].view(np.complex128)
                el = f_in['mf_header/kpoints/el'][:, :, b[1]:b[0]+1].mean(axis=-1)
                    
                num_bands_in = b[0] - b[1] + 1
                
                for ispb in range(nspbps_v):
                    if fname_phases is None:
                        # Phases are normalized to the DOS / sqrt(number_pseudobands)
                        phases = np.random.random((num_bands_in, nk,))
                        phases = np.exp(2 * np.pi * 1.0j * phases) / np.sqrt(float(nspbps_v))

                        phases_file['phases/coeffs'][b[1]:b[0]+1, :, ispb, :] = phases.view(np.float64).reshape((phases.shape + (2,)))
                    
                    else:
                        phases = phases_file['phases/coeffs'][b[1]:b[0]+1, :, ispb, :].view(np.complex128)[:,:,0]
                        
                        
                    if verbosity > 1:
                        logger.info(f'phases.shape for SPB {ib}: {phases.shape}')
                        logger.info(f'coeffs.shape for SPB {ib}: {coeffs.shape}')
                    
                    # Make sure we do a complex mult., and then view result as float
                    
                    spb = np.concatenate([np.tensordot(coeffs[:,:,cum_ngk[k]:cum_ngk[k+1]], phases[:,k], axes=(0, 0)) for k in range(nk)], axis=-2)
                    
                    f_out['wfns/coeffs'][ib, :, :] = spb.view(np.float64)
                    f_out['mf_header/kpoints/el'][:, :, ib] = el

                    ib -= 1
                    
    
    if qshift:
        f_out.close()
        phases_file.close()

        f_in.close()

        end = time.time()
        logger.info(f'Done! Time taken: {round(end-start,2)} sec')
        return
        
    logger.info('Creating {} conduction pseudobands'.format(nspb_c))
    ib = nb_out_v + nc 
    
    for b in blocks_c:
        if verbosity > 1:
            logger.info(f'slice_index, slice: {ib, b}')

        if single_band:
            band_avg = b[0] + (b[1] - b[0]) // 2
            f_out['wfns/coeffs'][ib, :, :] = (f_in['wfns/coeffs'][band_avg, :, :]
                                              * np.sqrt(float(b[1] - b[0] + 1)))
            f_out['mf_header/kpoints/el'][:, :, ib] = f_in['mf_header/kpoints/el'][:, :, band_avg].mean(axis=-1)
            ib += 1
        elif nspbps_c == 1:
            f_out['wfns/coeffs'][ib, :, :] = f_in['wfns/coeffs'][b[0]:b[1] + 1, :, :].sum(axis=0)
            f_out['mf_header/kpoints/el'][:, :, ib] = f_in['mf_header/kpoints/el'][:, :, b[0]:b[1] + 1].mean(axis=-1)
                
            ib += 1
            
        else:
            coeffs = f_in['wfns/coeffs'][b[0]:b[1] + 1, :, :, :].view(np.complex128)
            el = f_in['mf_header/kpoints/el'][:, :, b[0]:b[1] + 1].mean(axis=-1)
            num_bands_in = b[1] - b[0] + 1
            for ispb in range(nspbps_c):
                # no phases file for conduction states
                # Phases are normalized to the DOS / sqrt(number_pseudobands)
                phases = np.random.random((num_bands_in, nk,))
                phases = np.exp(2 * np.pi * 1.0j * phases) / np.sqrt(float(nspbps_c))
                
                if verbosity > 1:
                        logger.info(f'phases.shape for SPB {ib}: {phases.shape}')
                        logger.info(f'coeffs.shape for SPB {ib}: {coeffs.shape}')
                
                # Make sure we do a complex mult., and then view result as float
                
                spb = np.concatenate([np.tensordot(coeffs[:,:,cum_ngk[k]:cum_ngk[k+1]], phases[:,k], axes=(0, 0)) for k in range(nk)], axis=-2)
                print(spb.shape)

                f_out['wfns/coeffs'][ib, :, :] = spb.view(np.float64)
                f_out['mf_header/kpoints/el'][:, :, ib] = el

                ib += 1
                
    
    f_out.close()
    phases_file.close()

    f_in.close()
    
    end = time.time()
    logger.info(f'Done! Time taken: {round(end-start,1)} sec \n \n')



if __name__ == "__main__":
    import argparse
    
    desc = '''Constructs stochastic pseudobands given an input WFN. Bands with
    energies up to a protection window E0 are copied without modification, and
    states higher in energy are aggregated into subspaces and represented with
    pseudobands, each one spanning the energy range [En, En * E_frac], where En
    is the energy of the first state in the subspace, and E_frac is a constant.'''
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('--fname_in', help='Input WFN.h5, in HDF5 format', required=True)
    parser.add_argument('--fname_in_q', help='Input WFNq.h5, in HDF5 format', required=True)
    parser.add_argument('--fname_out', help='Output WFN.h5 with pseudobands, in HDF5 format', required=True)
    parser.add_argument('--fname_out_q', help='Output WFNq.h5 with pseudobands, in HDF5 format', required=True)
    parser.add_argument('--nv', type=int, default=-1,
                        help='Number of protected valence bands counting from VBM.')
    parser.add_argument('--nc', type=int, default=100,
                        help='Number of protected conduction bands counting from CBM.')
    parser.add_argument('--efrac_v', type=float, default=0.01,
                        help=('Accumulation window for valence slices, as a fraction of the energy of the '
                              'band in each subspace.'))
    parser.add_argument('--efrac_c_fine', type=float, default=None,
                        help=('Accumulation window for conduction slices with energies <= max_freq, as a fraction of the energy of the '
                              'band in each subspace.'))
    parser.add_argument('--efrac_c', type=float, default=0.01,
                        help=('Accumulation window for conduction slices with energies > max_freq, as a fraction of the energy of the '
                              'band in each subspace.'))
    parser.add_argument('--max_freq', type=float, default=1.0,
                        help=('Maximum energy (Ry) before coarse slicing kicks in for conduction SPBs. This should be at least the maximum frequency for which you plan to evaluate epsilon.'))
    parser.add_argument('--nspbps_v', type=int, default=1,
                        help=('Number of stochastic pseudobands per valence slice'))
    parser.add_argument('--nspbps_c', type=int, default=1,
                        help=('Number of stochastic pseudobands per conduction slice'))
    parser.add_argument('--copydirectly', default=True,
                        help=('Direct copying for protected bands. If False, then copying is done in chunks to limit memory usage. Set to False is you have a large number of protected bands.'))
    parser.add_argument('--verbosity', type=int, default=0,
                        help='Set verbosity level')
    parser.add_argument('--single_band', default=False, action='store_true',
                        help='Use a single band instead of a stochastic combination')
    # parser.add_argument('--fname_phases', default=None,
    #                     help=('Phases.h5 file, containing random coefficients for SPB construction. Should be consistent with all other parameters. Intended for use with WFNq calculation.'))
    
    args = parser.parse_args()
    
    logging.basicConfig(filename=vars(args)['fname_out'][0:-3]+'.log', level=logging.DEBUG, 
                    format='%(asctime)s %(levelname)s %(name)s %(message)s')
    logger=logging.getLogger(__name__)
    
    logger.info(vars(args))
    
    out = check_and_block(**vars(args))
    
    pseudoband(False, *out, **vars(args))
    
    fname_phases = 'phases'+vars(args)['fname_out'][0:-3]+'.h5'
    
    pseudoband(True, *out, fname_phases = fname_phases, **vars(args))

    
