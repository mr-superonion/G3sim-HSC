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
    
    
    # Add extinction file
    filename_extinction = os.path.join(dir, 'cosmos_extinction.fits')
    extinction_table = Table.read(filename_extinction)
    
    # Join tables on cosmos_index
    catalog     =   join(catalog, extinction_table,'IDENT')
    catalog.sort(['paired', 'ipos'])
    
    # Add column to make a WL flag for the simulation data:
    #     imag_cmodel - a_i < 24.5
    #     sigma_e > 0, < 0.4
    #     |e|<2
    #     resolution factor
    #     i-band SNR
    sig_e = get_sigma_e(catalog)
    e= np.sqrt(catalog['ishape_hsm_regauss_e1']**2 +
               catalog['ishape_hsm_regauss_e2']**2)

    catalog['weak_lensing_flag'] = ((get_imag(catalog) - catalog['a_i']) < 24.5) & \
        (sig_e > 0.0) & \
        (sig_e < 0.4) & \
        (e <= 2.) & \
        (catalog['ishape_hsm_regauss_resolution'] >= 0.3) & \
        (get_snr(catalog) >= 10.)

    # add a bunch of columns that get used for WL analysis - fake it for now!
    catalog['ishape_hsm_regauss_derived_sigma_e'] = catalog['ishape_hsm_regauss_sigma']
    catalog['ishape_hsm_regauss_derived_rms_e'] = 0.365

    wt =  1./(catalog['ishape_hsm_regauss_derived_rms_e']**2 + \
                 catalog['ishape_hsm_regauss_derived_sigma_e']**2)
    catalog['ishape_hsm_regauss_derived_shape_weight'] = wt
    # Unpleasant discovery: the lines below used to say "= 0" and it turns out that meant those
    # columns were ints, which led to all kinds of issues later on.  This fix is for the next time a
    # sims catalog gets turned into HDF5.
    catalog['ishape_hsm_regauss_derived_shear_bias_m'] = 0.
    catalog['ishape_hsm_regauss_derived_shear_bias_c1'] = 0.
    catalog['ishape_hsm_regauss_derived_shear_bias_c2'] = 0.
    
    
    return catalog


def get_snr(catalog):
    """
    This utility computes the S/N for each object in the catalog, either based on
    magnitude error for the real data or on cmodel_flux for the simulations.
    It does not impose any cuts and returns NaNs for invalid S/N values.
    In the case of GREAT3-like simulations, the noise rescaling factor is applied
    to match the data.
    """

    snr = catalog['modelfit_CModel_flux'] / (catalog['modelfit_CModel_fluxSigma']
    return snr


def get_imag(catalog):
    """
    This utility returns the i-band magnitude of the objects in the input data or
    simulation catalog. Does not apply any cuts and returns NaNs for invalid
    values.
    """
    
    return -2.5 * np.log10(catalog['modelfit_CModel_flux']) + 27.


def get_abs_ellip(catalog):
    """
    Returns the modulus of galaxy distortions.
    """
    return np.sqrt(catalog['ishape_hsm_regauss_e1']**2 +
                   catalog['ishape_hsm_regauss_e2']**2)


def get_res(catalog):
    """
    Returns the resolution
    """
    return catalog['ishape_hsm_regauss_resolution']
    
    
def get_sigma_e(catalog):
    """
    This utility returns the hsm_regauss_sigma values for the catalog, without
    imposing any additional flag cuts.
    In the case of GREAT3-like simulations, the noise rescaling factor is applied
    to match the data.
    """
    sigma_e = catalog['ishape_hsm_regauss_sigma']

    return sigma_e

