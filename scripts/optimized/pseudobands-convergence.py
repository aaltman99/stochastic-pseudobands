#!/usr/bin/env python

# ARA 06/2022

# Purpose: convergence testing for GW with SPB

# this script reads in a WFN_SPB.h5 file and removes 1 SPB from
# each slice. The user can specify whether to do this for only
# valence states, conduction states, or both. Protected bands are 
# untouched.

# Note that if order to perform an epsilon convergence test, 
# this script must be run on both WFN_SPB.h5 and WFN_SPB_q.h5.


import numpy as np
import h5py


# extracts pseudobands parameters from WFN_SPB.h5 file
def parse_pseudobands_params(wfn):
    
    val_params = {}
    cond_params = {}
    
    try:
        params = wfn['parabands/pseudobands']
        
        cond_params['nprot']= params['nc'][()]
        cond_params['nslice'] = params['n_subspaces'][()]
        cond_params['nspbps'] = params['num_per_subspace'][()]
        cond_params['max_freq'] = 0.0
        cond_params['uniform_width'] = 0.0
        
    except:
        params_c = wfn['pseudobands/conduction']
        params_v = wfn['pseudobands/valence']
    
        cond_params['nprot']= params_c['nprot'][()]
        cond_params['nslice'] = params_c['nslice'][()]
        cond_params['nspbps'] = params_c['nspbps'][()]
        cond_params['max_freq'] = params_c['max_freq'][()]
        cond_params['uniform_width'] = params_c['uniform_width'][()]

        val_params['nprot']= params_v['nprot'][()]
        val_params['nslice'] = params_v['nslice'][()]
        val_params['nspbps'] = params_v['nspbps'][()]
        val_params['max_freq'] = params_v['max_freq'][()]
        val_params['uniform_width'] = params_v['uniform_width'][()]

    return val_params, cond_params


def fill_pseudoband_params(f_out, val_params, cond_params):
    
    groupv = 'pseudobands/valence'
    groupc = 'pseudobands/conduction'
    
    for p in val_params.keys():
        f_out[groupv].create_dataset(p, (), data=val_params[p])
        f_out[groupc].create_dataset(p, (), data=cond_params[p])
    

# function to remove one SPB from each slice
def remove_spb(f_in, f_out, params, nb_out_v_orig, nb_out_v, nprot, vc):
    
    if vc == 'valence':
        ib = nb_out_v - nprot - 1
        # need separate index for the input wfn
        ib_in = nb_out_v_orig - nprot - 1
    elif vc == 'conduction':
        ib = nb_out_v + nprot
        # need separate index for the input wfn
        ib_in = nb_out_v_orig + nprot
    
    for idx_b in range(params['nslice']):
        for idx_spb in range(params['nspbps']): # referencing original nspbps

            # skip the last SPB in this slice, move on to next slice
            if idx_spb == params['nspbps'] - 1:
                # only ib_in gets iterated here
                if vc == 'valence':
                    ib_in -= 1
                elif vc == 'conduction':
                    ib_in += 1
                continue

            f_out['wfns/coeffs'][ib, :, :] = f_in['wfns/coeffs'][ib_in, :, :]
            f_out['mf_header/kpoints/el'][:, :, ib] = f_in['mf_header/kpoints/el'][:, :, ib_in]

            if vc == 'valence':
                ib -= 1
                ib_in -= 1
            elif vc == 'conduction':
                ib += 1
                ib_in += 1


# reads WFN_SPB.h5 file and creates a new WFN.h5 file with 1 band per slice removed
# updates pseudobands parameters as well
def reduce_wfn(fname_in, do_val, do_cond):
    
    assert do_val or do_cond, 'Reducing neither valence nor conduction states, output will be the same as input.' 
    
    f_in = h5py.File(fname_in, 'r')
    
    ifmax = f_in['mf_header/kpoints/ifmax'][()][0][0]
    
    val_params, cond_params = parse_pseudobands_params(f_in)
    
    if do_val:
        assert val_params['nspbps'] >= 2
    if do_cond:
        assert cond_params['nspbps'] >= 2

    val_params_new = val_params.copy()
    cond_params_new = cond_params.copy()
    if do_val:
        val_params_new['nspbps'] -= 1
    if do_cond:
        cond_params_new['nspbps'] -= 1
    
    print('Original valence SPB parameters: ', val_params)
    print('Original conduction SPB parameters: ', cond_params)
    print()
    print('New valence SPB parameters: ', val_params_new)
    print('New conduction SPB parameters: ', cond_params_new)
    
    nv = val_params['nprot']
    nc = cond_params['nprot']
    nb_out_v = nv + val_params_new['nspbps'] * val_params_new['nslice']
    nb_out_c = nc + cond_params_new['nspbps'] * cond_params_new['nslice']
     
    fname_out = fname_in[:-3]
    if do_val: 
        fname_out += '-val'
    if do_cond:
        fname_out += '-cond'
    fname_out += '-reduce.h5'
    f_out = h5py.File(fname_out, 'w')
    
    f_out.copy(f_in['mf_header'], 'mf_header')
    
    if do_val and do_cond:
        f_out['mf_header/kpoints/mnband'][()] -= val_params['nslice'] + cond_params['nslice']
        f_out['mf_header/kpoints/ifmax'][()] -= val_params['nslice']
    elif not do_val and do_cond:
        f_out['mf_header/kpoints/mnband'][()] -= cond_params['nslice']
    elif do_val and not do_cond:
        f_out['mf_header/kpoints/mnband'][()] -= val_params['nslice']
        f_out['mf_header/kpoints/ifmax'][()] -= val_params['nslice']
   
    
    
    f_out.create_group('pseudobands')
    f_out.create_group('pseudobands/conduction')
    f_out.create_group('pseudobands/valence')

    fill_pseudoband_params(f_out, val_params_new, cond_params_new)
        
    f_out.create_group('wfns')
    f_out.copy(f_in['wfns/gvecs'], 'wfns/gvecs')
    
    def resize(file, name):
        file.move(name, name + '_orig_pb')
        shape = list(file[name + '_orig_pb'].shape)
        shape[-1] = nb_out_v + nb_out_c
        file.create_dataset(name, shape, dtype='d')
        del file[name + '_orig_pb']
        
            
    resize(f_out, 'mf_header/kpoints/occ')
    resize(f_out, 'mf_header/kpoints/el')
    
    shape = list(f_in['wfns/coeffs'].shape)
    shape[0] = nb_out_v + nb_out_c
    f_out.create_dataset('wfns/coeffs', shape, 'd')
    
    # copy the unmodified states
    if do_val and do_cond:
        f_out['wfns/coeffs'][nb_out_v-nv:nb_out_v+nc, :, :] = f_in['wfns/coeffs'][ifmax-nv:ifmax+nc, :, :]
        f_out['mf_header/kpoints/el'][:,:,nb_out_v-nv:nb_out_v+nc] = f_in['mf_header/kpoints/el'][:,:,ifmax-nv:ifmax+nc]
    elif not do_val and do_cond:
        f_out['wfns/coeffs'][0:nb_out_v+nc, :, :] = f_in['wfns/coeffs'][0:ifmax+nc, :, :]
        f_out['mf_header/kpoints/el'][:,:,0:nb_out_v+nc] = f_in['mf_header/kpoints/el'][:,:,0:ifmax+nc]
    elif do_val and not do_cond:
        f_out['wfns/coeffs'][nb_out_v-nv:None, :, :] = f_in['wfns/coeffs'][ifmax-nv:None, :, :]
        f_out['mf_header/kpoints/el'][:,:,nb_out_v-nv:None] = f_in['mf_header/kpoints/el'][:,:,ifmax-nv:None]
        
        
    # deal with the modified states
    nb_out_v_orig = val_params['nprot'] + val_params['nspbps'] * val_params['nslice']
    if do_val:
        remove_spb(f_in, f_out, val_params, nb_out_v_orig, nb_out_v, nv, 'valence')   
    if do_cond:
        remove_spb(f_in, f_out, cond_params, nb_out_v_orig, nb_out_v, nc, 'conduction')
      
    
    # that's it, close the files
    f_in.close()
    f_out.close()
    
    return



if __name__ == "__main__":
    import argparse
    
    desc = '''Removes 1 PseudoBand from each slice in a given SPB WFN.h5. 
    Intended for simple convergence testing of the GW step with PseudoBands.'''
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('--fname_in', help='Input WFN.h5, in HDF5 format', required=True)
    
    def str_to_bool(value):
        if value.lower() in {'false', 'f', '0', 'no', 'n'}:
            return False
        elif value.lower() in {'true', 't', '1', 'yes', 'y'}:
            return True
        raise ValueError(f'{value} is not a valid boolean value')

    parser.add_argument('--do_val', help='Remove SPBs from valence states?', type=str_to_bool, nargs='?', const=True, default=False)
    parser.add_argument('--do_cond', help='Remove SPBs from conduction states?', type=str_to_bool, nargs='?', const=True, default=True)
    
    args = parser.parse_args()

    reduce_wfn(**vars(args))
    
  
