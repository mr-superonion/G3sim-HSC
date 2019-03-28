try:
    import astropy.io.fits as pyfits
except:
    import pyfits

from astropy.table import Table, join
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['ps.useafm'] = True
matplotlib.rcParams['pdf.use14corefonts'] = True
matplotlib.rcParams['text.usetex'] = True
matplotlib.use('Agg')
import math
import warnings


def prepare_sim_catalog(dir=None,subfield=None, rescale_est_variance=False):
    # Load the catalog from fits file
    filename = os.path.join(dir, 'outcomeFPFS/src-%04d.fits' % subfield)
    # If the file does not exist, an IOError will be raised
    catalog = Table.read(filename)
    
    # Add columns with names matching the real data
    cols = []
    data_fields = {'ext_shapeHSM_HsmShapeRegauss_e1': 'ishape_hsm_regauss_e1',
                   'ext_shapeHSM_HsmShapeRegauss_e2': 'ishape_hsm_regauss_e2',
                   'ext_shapeHSM_HsmShapeRegauss_sigma': 'ishape_hsm_regauss_sigma',
                   'ext_shapeHSM_HsmShapeRegauss_resolution': 'ishape_hsm_regauss_resolution'}
    for i in data_fields:
        catalog[data_fields[i]] = catalog[i]
    
    # Load the header to get proper name of flags
    header = pyfits.getheader(filename, 1)
    n_flag = catalog["flags"].shape[1]
    for i in range(n_flag):
        catalog[header["TFLAG%s" % (i+1)]] = catalog["flags"][:, i]
    
    # Load information from noiVarPsf
    noiPsfFname =   './noiVar_rPsf_s18a.csv'
    noiPsf      =   Table.read(noiPsfFname)
    inoiPsf     =   noiPsf[subfield]
    catalog['est_over_true_var'] = inoiPsf['var_meas'] / inoiPsf['var_target']
    catalog['reweight_s18a']     = inoiPsf['weightNoi']*inoiPsf['weightNoi']
    # Load information from psfPre
    psfFname    =   os.path.join(dir, 'psfPre/psf%04d.fits' % subfield)
    psfInfo     =   pyfits.getheader(psfFname)
    catalog['g1']=  psfInfo['g1']
    catalog['g2']=  psfInfo['g2']
    # Add a gridIndex column instead of using the index column, for compatibility
    catalog['ipos']     =   (catalog['base_SdssCentroid_y']//64)*50 +(catalog['base_SdssCentroid_x']//64)
    catalog['ipos']     =   catalog['ipos'].astype(np.int)
    catalog['centDist']  =   ((catalog['base_SdssCentroid_y']%64-32)**2. +
                            (catalog['base_SdssCentroid_x']%64-32)**2.)
    catalog['centDist']  =   np.sqrt(catalog['centDist'])
    # First, keep only detections that are the closest to the grid point
    # Get sorted index by grid index and grid distance
    catalog     =   catalog[catalog['deblend_nChild']==0]
    inds        =   np.lexsort([catalog['centDist'], catalog['ipos']])
    inds_unique =   np.unique(catalog['ipos'][inds], return_index=True)[1]
    catalog     =   catalog[inds[inds_unique]]
    
    # Then, apply mask for permissive cuts
    mask =  (~(catalog['base_SdssCentroid_flag'])) & \
            (~catalog['ext_shapeHSM_HsmShapeRegauss_flag']) & \
            (catalog['base_ClassificationExtendedness_value'] != 0) & \
            (~np.isnan(catalog['ext_shapeHSM_HsmShapeRegauss_sigma']))
    catalog = catalog[mask]


    # Load information from catPre
    catFname    =   os.path.join(dir, 'catPre/catalog-%04d.fits' % subfield)
    catInfo     =   Table.read(catFname)
    # Join catalogs on grid index
    catalog = join(catalog, catInfo, 'ipos')
    # Now, identify pairs based on their input ipos values
    n = 50*50 # Number of grid points
    mask_even = np.in1d(catalog['ipos'] % n, 2*np.arange(n/2).astype('int'), assume_unique=True)
    mask_odd  = np.in1d(catalog['ipos'] % n, 2*np.arange(n/2).astype('int')+1, assume_unique=True)

    # Find index of broken pairs
    broken_pairs = np.setxor1d(catalog['ipos'][mask_even] % n,
                               catalog['ipos'][mask_odd] % n - 1)
    mask_pairs = ~np.in1d(catalog['ipos'] % n,
                              np.concatenate([broken_pairs, broken_pairs+1]), assume_unique=True)
    catalog['paired'] = mask_pairs
    catalog.sort(['paired', 'ipos'])
    return catalog


