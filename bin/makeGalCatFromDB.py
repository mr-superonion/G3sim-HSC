#!/usr/bin/env python
import os
import gc
import numpy as np
import healpy as hp
from sys import getsizeof
import astropy.table as astTable
import astropy.io.ascii as ascii
import sim_analysis_utilities as utilities


prefix      =   'database/s18/tracts/'
outfix      =   'database/s18/fields/'  
tractList   =   ascii.read('fieldTractInfo.csv')['tract']
fieldList0  =   ascii.read('fieldTractInfo.csv')['field']
fieldList   =   utilities.field_keys.keys()

if not os.path.exists(outfix):
    os.system('mkdir -p %s' %outfix)

def main():
    for field in fieldList:
        processField(field)
    return

def processField(field):
    outfname=   os.path.join(outfix,'%s.fits' %field)
    print('processing field %s:' %field)
    tractLU =   tractList[fieldList0==field]
    dataList=   list(map(process,tractLU))
    dataList=   [dd for dd in dataList if dd is not None]
    dataT   =   astTable.vstack(dataList)
    dataT.write(outfname)
    del outfname
    del dataList
    del dataT
    gc.collect()
    return


def process(tract):
    print('processing tract: %s' %tract)

    """
    # Read data
    """
    infname =   os.path.join(prefix,'%s.fits' %tract)
    if not os.path.exists(infname):
        return None
    fits    =   astTable.Table.read(infname)
    print(getsizeof(fits))
    colns   =   fits.colnames
    colns2  =   [cn for cn in colns if "_isnull" not in cn]
    assert len(colns2) == len(colns)//2, \
        "*_insull columns are not properly removed."
    fits    =   fits[colns2]
    del colns
    del colns2

    """
    # Get npass
    """
    fits['npass']  =   getNPass(fits)

    """
    # FDFC, WL, BrightStar Mask
    """
    sel     =   selectGalaxyCatalog(fits)

    return fits[sel] 

def selectGalaxyCatalog(d):
    # apply FDFC mask
    m = hp.read_map("/gpfs02/work/miyatake/S18A_fdfc/s18a_fdfc_hp_contarea.fits", nest = True, dtype = np.bool)
    indices_map = np.where(m)[0]
    nside = hp.get_nside(m)
    phi = np.radians(d["i_ra"])
    theta = np.radians(90. - d["i_dec"])
    indices_obj = hp.ang2pix(nside, theta, phi, nest = True)
    sel = np.in1d(indices_obj, indices_map)
    # weak lensing flag
    sel = sel & d["weak_lensing_flag"]
    #bright star mask
    sel = sel & (d["i_mask_s18a_bright_objectcenter"] == False)
    return sel

def getNPass(fits):
    # multi-band detection to remove junk
    min_multiband_snr_data  = 5.0
    g_snr   =   (2.5/np.log(10.0)) / fits['forced_g_cmodel_magsigma']
    r_snr   =   (2.5/np.log(10.0)) / fits['forced_r_cmodel_magsigma']
    z_snr   =   (2.5/np.log(10.0)) / fits['forced_z_cmodel_magsigma']
    y_snr   =   (2.5/np.log(10.0)) / fits['forced_y_cmodel_magsigma']
    gcountinputs=   fits['g_inputcount_value']
    rcountinputs=   fits['r_inputcount_value']
    zcountinputs=   fits['z_inputcount_value']
    ycountinputs=   fits['y_inputcount_value']
    g_mask  =   (g_snr >= min_multiband_snr_data) & (~np.isnan(g_snr) & (gcountinputs >= 2))
    r_mask  =   (r_snr >= min_multiband_snr_data) & (~np.isnan(r_snr) & (rcountinputs >= 2))
    z_mask  =   (z_snr >= min_multiband_snr_data) & (~np.isnan(z_snr) & (zcountinputs >= 2))
    y_mask  =   (y_snr >= min_multiband_snr_data) & (~np.isnan(y_snr) & (ycountinputs >= 2))
    n_pass  =   g_mask.astype(int) + r_mask.astype(int) + z_mask.astype(int) + y_mask.astype(int)
    return n_pass

if __name__=='__main__':
    main()
