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
    # Add a gridIndex column instead of using the index column, for compatibility
    catalog['ipos']     =   ((catalog['base_SdssCentroid_y']//64)*100 +
                            (catalog['base_SdssCentroid_x']//64))
    catalog['cenDist']  =   ((catalog['base_SdssCentroid_y']%64-32)**2. +
                            (catalog['base_SdssCentroid_x']%64-32)**2.)
    catalog['cenDist']  =   np.sqrt(catalog['cenDist'])
    # First, keep only detections that are the closest to the grid point
    # Get sorted index by grid index and grid distance
    inds        =   np.lexsort([catalog['centDist'], catalog['ipos']])
    inds_unique =   np.unique(catalog['ipos'][inds], return_index=True)[1]
    catalog     =   catalog[inds[inds_unique]]
    
    return catalog


