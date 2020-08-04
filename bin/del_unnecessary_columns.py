#!/usr/bin/env python
import gc
import os

from astropy.table import Table, join
import sim_analysis_utilities as utilities

# Number of blinded catalogs to produce
key_info    =   'xiangchong.li@ipmu.jp'
blind_level =   1
n_blind_cat =   3

inDir   =       './catalog_obs_reGaus_calibrated/fields/'
outDir  =       './catalog_obs_reGaus_blinded/fields/level-%d/' %(blind_level)
if not os.path.exists(outDir):
    os.system('mkdir -p %s' %outDir)

del_columns = [
    'i_hsmshaperegauss_derived_bias_m',
    'i_hsmshaperegauss_derived_bias_c1',
    'i_hsmshaperegauss_derived_bias_c2']

def main():
    for fieldname in utilities.field_names:
        fname =   os.path.join(inDir,'%s_calibrated.fits' %(fieldname))
        data    =   Table.read(fname)
        for entry in del_columns:
            del data[entry]
        data.write(fname,overwrite=True)
        del data
        del fname
        gc.collect()
        for ib in range(n_blind_cat):
            fname   =   os.path.join(outDir,'%s_blinded_%d.fits' %(fieldname,ib))
            data    =   Table.read(fname)
            for entry in del_columns:
                del data[entry]
            data.write(fname,overwrite=True)
            del data
            del fname
            gc.collect()
    return

if __name__ == "__main__":
    main()
