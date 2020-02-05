#!/usr/bin/env python
import os
import gc
import glob
import numpy as np
import healpy as hp
from sys import getsizeof
import astropy.table as astTable
import astropy.io.ascii as ascii
import sim_analysis_utilities as utilities

version     =   's18'
prefix      =   'database/%s/tracts/'%version
outfix      =   'database/%s/fields/' %version
fieldList   =   utilities.field_keys.values()
if 's16' in version:
    infoname=   'fieldTractInfoS16.csv'
elif 's17' in version:
    infoname=   'fieldTractInfoS16.csv'
elif 's18' in version:
    infoname=   'fieldTractInfo.csv'
elif 's19' in version:
    infoname=   'fieldTractInfoS16.csv'

tractList   =   ascii.read(infoname)['tract']
fieldList0  =   ascii.read(infoname)['field']

def main():
    if not os.path.isdir(prefix):
        raise RuntimeError('Cannot find input directory')
        return
    if not os.path.exists(outfix):
        os.system('mkdir -p %s' %outfix)
    for field in fieldList:
        processField(field)
    return

def processField(field):
    outfname=   os.path.join(outfix,'%s.fits' %field)
    if os.path.exists(outfname):
        print('already have output for field: %s' %field)
    print('processing field %s:' %field)
    tractLU =   tractList[fieldList0==field]
    dataList=   list(map(process,tractLU))
    dataList=   [dd for dd in dataList if dd is not None]
    dataT   =   astTable.vstack(dataList)
    dataT.write(outfname,overwrite=True)
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
    """
    # Load diff Patch
    """
    pDirL   =   glob.glob('database/%s/tracts-diff*/' %version)
    for pp in pDirL:
        pfname  =   os.path.join(pp,'%s.fits' %tract)
        pf      =   astTable.Table.read(pfname)
        fits    =   astTable.join(fits,pf)

    """
    # Remove _isnull columns
    """
    colns   =   fits.colnames
    colns2  =   [cn for cn in colns if "_isnull" not in cn]
    assert len(colns2) == len(colns)//2, \
        "*_insull columns are not properly removed."
    fits    =   fits[colns2]
    del colns
    del colns2

    """
    # FDFC, WL, BrightStar Mask
    """
    #sel     =   selectGalaxyCatalog(fits)
    return fits

def selectGalaxyCatalog(d):
    # apply FDFC mask
    m = hp.read_map("database/s18a_fdfc_hp_contarea.fits", nest = True, dtype = np.bool)
    indices_map = np.where(m)[0]
    nside = hp.get_nside(m)
    phi = np.radians(d["i_ra"])
    theta = np.radians(90. - d["i_dec"])
    indices_obj = hp.ang2pix(nside, theta, phi, nest = True)
    sel = np.in1d(indices_obj, indices_map)
    # weak lensing flag
    sel = sel & d["weak_lensing_flag"]
    #bright star mask
    if 's18' in version or 's19' in version:
        sel = sel & (d["i_mask_s18a_bright_objectcenter"] == False)
    return sel

if __name__=='__main__':
    main()
