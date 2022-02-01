#!/usr/bin/env python

# Pseudobands!
# It's recommended that you install the "progress" python script with:
# "pip install --user progress"
#
# Felipe H da Jornada (2015)

from __future__ import print_function
import h5py
import numpy as np



try:
    from progress.bar import Bar
except:
    Bar = None


def pseudoband(fname_in, fname_out, e0, efrac, nspbps, single_band):
    if single_band:
        nspbps = 1

    f_in = h5py.File(fname_in, 'r')
    f_out = h5py.File(fname_out, 'w')

    print('Pseudoband parameters:')
    print('  Input file: {}'.format(fname_in))
    print('  Output file: {}'.format(fname_out))
    print('  Protection window: {} Ry'.format(e0))
    print('  Aggregation energy fraction: {} %'.format(efrac * 100))
    print('  Number of stochastic pseudobands per slice: {}'.format(nspbps))
    print()

    en_orig = f_in['mf_header/kpoints/el'][()]
    nb_orig = f_in['mf_header/kpoints/mnband'][()]
    print('Original number of bands: {}'.format(nb_orig))

    # Ok, so e0 and e_max doesn`t really work, we just take the average
    en_min = np.mean(en_orig[0, :, :], axis=0)
    en_max = np.mean(en_orig[0, :, :], axis=0)

    ind = np.arange(nb_orig)
    nb_fixed = ind[en_min > e0][0]
    print('Number of bands in the protection window: {}'.format(nb_fixed))

    blocks = []
    nb_out = nb_fixed
    first_idx = nb_fixed
    while True:
        first_en = en_min[first_idx]
        delta_en = first_en * efrac
        last_en = first_en + delta_en
        try:
            last_idx = ind[en_max > last_en][0]
            blocks.append((first_idx, last_idx - 1))
            nb_out += 1
            first_idx = last_idx
        except:
            last_idx = nb_orig
            blocks.append((first_idx, last_idx - 1))
            nb_out += 1
            break

    nslices = nb_out - nb_fixed
    nspb = nslices * nspbps
    nb_out = nb_fixed + nspb
    assert len(blocks) == nslices
    print('Number of slices: {}'.format(nslices))
    print('Number of pseudobands: {}'.format(nspb))
    print('Number of bands in the output file: {}'.format(nb_out))
    f_out.copy(f_in['mf_header'], 'mf_header')
    f_out.create_group('wfns')
    f_out.copy(f_in['wfns/gvecs'], 'wfns/gvecs')
    f_out['mf_header/kpoints/mnband'][()] = nb_out

    def resize(name):
        f_out.move(name, name + '_orig_pb')
        shape = list(f_out[name + '_orig_pb'].shape)
        shape[-1] = nb_out
        f_out.create_dataset(name, shape, dtype='d')
        f_out[name][:, :, :nb_out] = f_out[name + '_orig_pb'][:, :, :nb_out]
        del f_out[name + '_orig_pb']

    resize('mf_header/kpoints/occ')
    resize('mf_header/kpoints/el')

    print('Copying {} protected bands'.format(nb_fixed))
    shape = list(f_in['wfns/coeffs'].shape)
    shape[0] = nb_out
    f_out.create_dataset('wfns/coeffs', shape, 'd')
    # print(shape)
    if False:
        # FHJ: this can take quite a while for some systems..
        f_out['wfns/coeffs'][:nb_fixed, :, :] = f_in['wfns/coeffs'][:nb_fixed, :, :]
    else:
        nbs_block = 1000
        ibs_start = range(0, nb_fixed, nbs_block)
        if Bar is not None:
            bar = Bar('Copying protected bands', max=len(ibs_start), bar_prefix=' [', bar_suffix='] ',
                      fill='#', suffix='%(percent)d%% - Remaining: %(eta_td)s')
        for ib in ibs_start:
            if Bar is not None:
                bar.next()
            ib1, ib2 = ib, min(ib + nbs_block, nb_fixed)
            tmp = f_in['wfns/coeffs'][ib1:ib2, :, :, :]
            f_out['wfns/coeffs'][ib1:ib2, :, :] = tmp
    if Bar is not None:
        bar.finish()
    
    
    break
    
    
    print('Creating {} pseudobands'.format(nspb))
    ib = nb_fixed
    if Bar is not None:
        bar = Bar('Creating superbands', max=nb_orig - nb_fixed, bar_prefix=' [', bar_suffix='] ',
                  fill='#', suffix='%(percent)d%% - Remaining: %(eta_td)s')
    for b in blocks:
        if Bar is not None:
            bar.next(b[1] - b[0] + 1)
        if single_band:
            band_avg = b[0] + (b[1] - b[0]) // 2
            f_out['wfns/coeffs'][ib, :, :] = (f_in['wfns/coeffs'][band_avg, :, :]
                                              * np.sqrt(float(b[1] - b[0] + 1)))
            f_out['mf_header/kpoints/el'][:, :, ib] = f_in['mf_header/kpoints/el'][:, :, band_avg].mean(axis=-1)
            ib += 1
        elif nspbps == 1:
            f_out['wfns/coeffs'][ib, :, :] = f_in['wfns/coeffs'][b[0]:b[1] + 1, :, :].sum(axis=0)
            f_out['mf_header/kpoints/el'][:, :, ib] = f_in['mf_header/kpoints/el'][:, :, b[0]:b[1] + 1].mean(axis=-1)
            ib += 1
        else:
            coeffs = f_in['wfns/coeffs'][b[0]:b[1] + 1, :, :, :].view(np.complex128)
            el = f_in['mf_header/kpoints/el'][:, :, b[0]:b[1] + 1].mean(axis=-1)
            num_bands_in = b[1] - b[0] + 1
            for ispb in range(nspbps):
                # Phases are normalized to the DOS / sqrt(number_pseudobands)
                phases = np.random.random((num_bands_in,))
                phases = np.exp(2 * np.pi * 1.0j * phases) / np.sqrt(np.float(nspbps)) ### ARA: added the 2 pi
                # Make sure we do a complex mult., and then view result as float
                spb = np.tensordot(coeffs, phases, axes=(0, 0))
                f_out['wfns/coeffs'][ib, :, :] = spb.view(np.float64)
                f_out['mf_header/kpoints/el'][:, :, ib] = el
                ib += 1
    if Bar is not None:
        bar.finish()
    print('')
    print('All done!')

    f_out.close()
    f_in.close()


if __name__ == "__main__":
    import argparse
    
    desc = '''Constructs stochastic pseudobands given an input WFN. Bands with
    energies up to a protection window E0 are copied without modification, and
    states higher in energy are aggregated into subspaces and represented with
    pseudobands, each one spanning the energy range [En, En * E_frac], where En
    is the energy of the first state in the subspace, and E_frac is a constant.'''
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('fname_in', help='Input WFN.h5, in HDF5 format')
    parser.add_argument('fname_out', help='Output WFN.h5 with pseudobands, in HDF5 format')
    parser.add_argument('--e0', type=float, default=4.,
                        help='Protection window, in Ry, measured as an absolute energy.')
    parser.add_argument('--efrac', type=float, default=0.02,
                        help=('Accumulation window, as a fraction of the energy of the '
                              'band in each subspace.'))
    parser.add_argument('--nspbps', type=int, default=1,
                        help=('Number of stochastic pseudobands per slice'))
    parser.add_argument('--single_band', default=False, action='store_true',
                        help='Use a single band instead of a stochastic combination')
    args = parser.parse_args()

    pseudoband(**vars(args))

    
