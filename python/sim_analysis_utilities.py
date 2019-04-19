try:
    import astropy.io.fits as pyfits
except:
    import pyfits

from astropy.table import Table, join
import os
import numpy as np
import matplotlib
matplotlib.rcParams['ps.useafm'] = True
matplotlib.rcParams['pdf.use14corefonts'] = True
matplotlib.rcParams['text.usetex'] = True
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import math
import warnings

# Field keys used in the noise variance file for simulations
field_keys = {'aegis': 'AEGIS',
              'gama9': 'GAMA09H',
              'gama15': 'GAMA15H',
              'hect': 'HECTOMAP',
              'vvds': 'VVDS',
              'wide12': 'WIDE12H',
              'xmm': 'XMM'}
# Hex keys for the color scheme for shear catalog paper
plot_colors = {
    'GAMA09H' : '#d73027',
    'GAMA15H' : '#fc8d59',
    'HECTOMAP' : '#fee090',
    'VVDS' : '#000000',
    'WIDE12H' : '#91bfdb',
    'XMM' : '#4575b4'
    }

def prepare_sim_catalog(catalog,isSim=True):
    colNames    =   catalog.colnames
    
    # Load the header to get proper name of flags
    header = pyfits.getheader(filename, 1)
    n_flag = catalog["flags"].shape[1]
    for i in range(n_flag):
        catalog[header["TFLAG%s" % (i+1)]] = catalog["flags"][:, i]
    
    # Then, apply mask for permissive cuts
    mask =  (~(catalog['base_SdssCentroid_flag'])) &\
            (~catalog['ext_shapeHSM_HsmShapeRegauss_flag']) &\
            (catalog['base_ClassificationExtendedness_value'] > 0) &\
            (~np.isnan(catalog['modelfit_CModel_flux']))&\
            (~np.isnan(catalog['modelfit_CModel_fluxSigma'])) &\
            (~np.isnan(catalog['ext_shapeHSM_HsmShapeRegauss_resolution'])) &\
            (~np.isnan(catalog['ext_shapeHSM_HsmPsfMoments_xx'])) &\
            (~np.isnan(catalog['ext_shapeHSM_HsmPsfMoments_yy'])) &\
            (~np.isnan(catalog['ext_shapeHSM_HsmPsfMoments_xy'])) &\
            (~np.isnan(catalog['base_Variance_value']))&\
            (~np.isnan(catalog['modelfit_CModel_flux']))&\
            (~np.isnan(catalog['modelfit_CModel_fluxSigma']))&\
            (~np.isnan(catalog['ext_shapeHSM_HsmShapeRegauss_resolution']))
    catalog =   catalog[colNames]
    # Add columns with names matching the real data
    cols = []
    data_fields = {'ext_shapeHSM_HsmShapeRegauss_e1': 'ishape_hsm_regauss_e1',
                   'ext_shapeHSM_HsmShapeRegauss_e2': 'ishape_hsm_regauss_e2',
                   'ext_shapeHSM_HsmShapeRegauss_sigma': 'ishape_hsm_regauss_sigma',
                   'ext_shapeHSM_HsmShapeRegauss_resolution': 'ishape_hsm_regauss_resolution'}
    for i in data_fields:
        catalog[data_fields[i]] = catalog[i]
    
    # For simulation
    if isSim:
        catalog['reweight_s16a']    =  0.
        # Load information from psfPre
        psfFname    =   os.path.join(dir, 'psfPre/psf%04d.fits' % isim)
        psfInfo     =   pyfits.getheader(psfFname)
        catalog['g1']=  psfInfo['g1']
        catalog['g2']=  psfInfo['g2']
        catalog['field_name']   = np.string_(field_keys[psfInfo['field']])
        # Add a gridIndex column instead of using the index column, for compatibility
        catalog['ipos']     =   (catalog['base_SdssCentroid_y']//64)*50 +(catalog['base_SdssCentroid_x']//64)
        catalog['ipos']     =   catalog['ipos'].astype(np.int)
        catalog['centDist'] =   ((catalog['base_SdssCentroid_y']%64-32)**2. +
                                (catalog['base_SdssCentroid_x']%64-32)**2.)
        catalog['centDist'] =   np.sqrt(catalog['centDist'])
        # First, keep only detections that are the closest to the grid point
        # Get sorted index by grid index and grid distance
        inds        =   np.lexsort([catalog['centDist'], catalog['ipos']])
        inds_unique =   np.unique(catalog['ipos'][inds], return_index=True)[1]
        catalog     =   catalog[inds[inds_unique]]
        mask    =   mask&(catalog['deblend_nChild']==0)&(catalog['centDist']<5.)
        # Load information from catPre
        catFname    =   os.path.join(dir, 'catPre/catalog-%04d.fits' % isim)
        catInfo     =   Table.read(catFname)
        # Join catalogs on grid index
        catalog = join(catalog, catInfo, 'ipos')
        catalog['BAND'].dtype   =   'S20'
        # Add extinction file
        #filename_extinction = os.path.join(dir, 'cosmos_extinction.fits')
        #extinction_table = Table.read(filename_extinction)
        
        # Join tables on cosmos_index
        #catalog     =   join(catalog, extinction_table,'IDENT')   
           
    catalog =   catalog[mask]
    if len(catalog)==0:
        return None
    
    if isSim:
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
    
    # Add column to make a WL flag for the simulation data:
    #     imag_cmodel - a_i < 24.5
    #     sigma_e > 0, < 0.4
    #     |e|<2
    #     resolution factor
    #     i-band SNR
    sig_e = get_sigma_e(catalog)
    e= np.sqrt(catalog['ishape_hsm_regauss_e1']**2 +
               catalog['ishape_hsm_regauss_e2']**2)

    catalog['weak_lensing_flag'] = ((get_imag(catalog)) < 24.5) & \
        (e <= 2.) & \
        (catalog['ishape_hsm_regauss_resolution'] >= 0.3) & \
        (get_snr(catalog) >= 10.)
        #(sig_e > 0.0) & \
        #(sig_e < 0.4) &

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


def galaxy_selector(catalog, min_snr=None, max_snr=None,
                    min_res=None, max_res=None,
                    max_sigma_e=None, max_e=None, max_mag=None,
                    min_multiband_snr_data=None, n_pass_multiband=2,
                    cancel_shape_noise=True, bright_star_mask=False,
                    cut_griddist=5.0, background_cut=-0.015,
                    use_wl_flag=False, max_blendedness=None,
                    applyboth="", use_res_round=False):
    """
    This script takes a catalog and returns a boolean NumPy
    array with a mask that says which galaxies to use for analysis.

    Warning: Cuts are applied in extinction corrected magnitude. For simulations,
    extinction values are extracted from the COSMOS field and correpond to the
    extinction at the position of the COSMOS galaxies used as the input of the
    simulation.

    Parameters
    ----------

    catalog: recarray
        Catalog on which the cuts apply.

    min_snr: float (optional)
        Minimum i-band unforced cmodel SNR cut

    max_snr: float (optional)
        Maximum i-band unforced cmodel SNR cut.

    min_res: float (optional)
        Minimum resolution cut

    max_res: float (optional)
        Maximum resolution cut

    max_sigma_e: float (optional)
        Maximum sigma_e cut. High values can result from various
        measurement failures.

    max_e: float (optional)
        Maximum ellipticity cut (2 is a good number to use).

    max_mag: float (optional)
        i-band magnitude cut.

    min_multiband_snr_data: float (optional)
        Minimum multiband snr cut; can only be applied in data, and will be
        ignored otherwise (require snr above this for r or z or y). Uses forced
        cmodel photometry for these bands.

    cancel_shape_noise: bool (default: True)
        If set to true, impose selection in a way that preserves 90 degree
        rotated pais. Only applies to the simulated data loaded with
        `load_sims_catalog`.

    bright_star_mask: bool (default: False)
        Impose bright star mask?

    cut_griddist: float (default: 5.0)
        Require simulation matches to lie within this many pixels of the
        expected location.  Only really relevant for cases where there are
        off-centered detections and the HSC pipeline is run with detection and
        deblending on (as in e.g. v4).

    use_wl_flag: bool (default: False)
        Use all of the default WL cuts?

    max_blendedness: float (optional)
        Maximum value of iblendedness_abs_flux

    background_cut: float (optional)
        Maximum allowed value of abs(background estimate in postage stamp)

    applyboth: string (optional)
        Which cuts should be applied to both galaxies instead of just 1?  To be
        used for tests of selection bias.  If "" then this is not done.

    use_res_round: bool (default: False)
        Apply resolution cut using round-ified version, which is always stricter?

    Notes
    -----

    There are some subtleties here: if we want to take advantage of shape
    noise cancellation, then we have to fairly impose the selection on both
    members of each 90 degree rotated pair (which are always adjacent to
    each other in the catalog -- i.e., catalog indices 0 and 1 are 90 degree
    rotated pairs, as are indices 2 and 3, etc.).

    `cancel_shape_noise` is used to decide whether we impose the selection
    in a way that maintains 90 degree rotated pairs or not. This option only
    """
    mask = np.ones(len(catalog)).astype(bool)
    # Tests if the catalog contains the `paired` column which allows us to
    # apply cuts and keep the noise cancellation.
    cols = catalog.dtype.names
    if cancel_shape_noise and 'paired' not in cols:
        warnings.warn("Noise cancellation requested on real data")
        cancel_shape_noise = False

    if cancel_shape_noise:
        ind_pair = np.where(catalog['paired'])[0]
        ind_pair1 = ind_pair[0::2]
        ind_pair2 = ind_pair[1::2]
        # Masking elements not in a pair to begin with
        mask &= catalog['paired']

    # Basic sanity checks: none of the quantities that we want to use should
    # be NaN. We also impose the sigma_e and e cuts here, because even just
    # a few bad values can totally mess stuff up.
    mask &= (~np.isnan(catalog['ishape_hsm_regauss_e1'])) & \
            (~np.isnan(catalog['ishape_hsm_regauss_e2'])) & \
            (~np.isnan(catalog['ishape_hsm_regauss_sigma'])) & \
            (~np.isnan(catalog['ishape_hsm_regauss_resolution'])) & \
            (~np.isnan(get_snr(catalog))) & \
            (~np.isnan(get_imag(catalog)))
    # For real data, also check valid psf measurements
    if 'ishape_sdss_psf_ixx' in catalog.dtype.names:
        mask &= (~np.isnan(catalog['ishape_sdss_psf_ixx'])) & \
                (~np.isnan(catalog['ishape_sdss_psf_iyy'])) & \
                (~np.isnan(catalog['ishape_sdss_psf_iyy']))

    if max_sigma_e is not None:
        sigma_e = get_sigma_e(catalog[mask])
        mask[mask] &= sigma_e <= max_sigma_e

    if max_e is not None:
        e = np.sqrt(catalog['ishape_hsm_regauss_e1'][mask]**2 +
                    catalog['ishape_hsm_regauss_e2'][mask]**2)
        mask[mask] &= (e <= max_e)

    if min_multiband_snr_data is not None:
        try:
            # try to access other-band SNR info.  If this is a sim, it won't
            # find that info and will throw an exception, so we know to skip it
            g_snr = (2.5/np.log(10.0)) / catalog['gmag_forced_cmodel_err'][mask]
            r_snr = (2.5/np.log(10.0)) / catalog['rmag_forced_cmodel_err'][mask]
            z_snr = (2.5/np.log(10.0)) / catalog['zmag_forced_cmodel_err'][mask]
            y_snr = (2.5/np.log(10.0)) / catalog['ymag_forced_cmodel_err'][mask]
            gcountinputs = catalog['gcountinputs'][mask]
            rcountinputs = catalog['rcountinputs'][mask]
            zcountinputs = catalog['zcountinputs'][mask]
            ycountinputs = catalog['ycountinputs'][mask]
            g_mask = (g_snr >= min_multiband_snr_data) & (~np.isnan(g_snr) & (gcountinputs >= 2))
            r_mask = (r_snr >= min_multiband_snr_data) & (~np.isnan(r_snr) & (rcountinputs >= 2))
            z_mask = (z_snr >= min_multiband_snr_data) & (~np.isnan(z_snr) & (zcountinputs >= 2))
            y_mask = (y_snr >= min_multiband_snr_data) & (~np.isnan(y_snr) & (ycountinputs >= 2))
            n_pass = g_mask.astype(int) + r_mask.astype(int) + \
                z_mask.astype(int) + y_mask.astype(int)
            mask[mask] &= (n_pass >= n_pass_multiband)
        except:
            warnings.warn("Attempted to apply multiband snr cut on simulations.")

    if bright_star_mask:
        try:
            mask[mask] &= ~catalog['iflags_pixel_bright_object_center'][mask]
        except:
            warnings.warn("Attempted to apply bright star mask on simulations.")

    if cut_griddist > 0:
        # Only do this cut on sims.
        if 'gridDist' in catalog.dtype.names:
            mask[mask] &= catalog['gridDist'][mask] < cut_griddist

    # The sanity checks are systematically applied to both members of the pair
    # Since NaN's make the data unusable, we have to self-consistently impose
    # this across both 90-degree-rotated pairs (if we're using shape noise
    # cancellation):
    if cancel_shape_noise:
        mask[ind_pair1] = mask[ind_pair1] & mask[ind_pair2]
        mask[ind_pair2] = mask[ind_pair1]

    if use_wl_flag:
        # Technically, this flag only exists in the data, but we created it for
        # the sims in load_sims_catalog.
        mask[mask] &= catalog['weak_lensing_flag'][mask]

    if use_res_round:
        mask[mask] &= get_res_round(catalog[mask]) >= 0.3

    # Now we move on to the cuts on resolution, S/N and magnitude.
    # We should only try to make the cuts on the galaxies that passed the flag
    # cuts, otherwise we get lots of warnings about invalid values etc.
    if min_res is not None:
        mask[mask] &= catalog['ishape_hsm_regauss_resolution'][mask] >= min_res

    if max_res is not None:
        mask[mask] &= catalog['ishape_hsm_regauss_resolution'][mask] <= max_res

    if min_snr is not None:
        snr = get_snr(catalog[mask])
        mask[mask] &= snr >= min_snr

    if max_snr is not None:
        snr = get_snr(catalog[mask])
        mask[mask] &= snr <= max_snr

    if max_mag is not None:
        mag = get_imag(catalog[mask]) - catalog['a_i'][mask]
        mask[mask] &= mag <= max_mag

    if max_blendedness is not None:
        mask[mask] &= catalog['iblendedness_abs_flux'][mask] <= max_blendedness

    # For sims only.
    if background_cut is not None and 'background' in catalog.dtype.names:
        mask[mask] &= catalog['background'][mask] > background_cut

    # For noise cancellation, apply the cut only to the first element of the
    # pair to avoid funny selection biases.
    # But enable introduction of selection bias due to single cuts through applyboth.
    if cancel_shape_noise:
        mask[ind_pair2] = mask[ind_pair1]
        if applyboth != "":
            if 'resolution' in applyboth:
                if min_res is None and use_wl_flag:
                    min_res = 0.3
                if min_res is not None:
                    print('Applying resolution cut to both objects in pair!')
                    mask[mask] &= (get_res(catalog[mask]) >= min_res)
            if 'snr' in applyboth:
                if min_snr is None and use_wl_flag:
                    min_snr = 10.0
                if min_snr is not None:
                    print('Applying SNR cut to both objects in pair!')
                    snr = get_snr(catalog[mask])
                    mask[mask] &= snr >= min_snr
            if 'mag' in applyboth:
                if max_mag is None and use_wl_flag:
                    max_mag = 24.5
                if max_mag is not None:
                    print('Applying magnitude cut to both objects in pair!')
                    mag = get_imag(catalog[mask]) - catalog['a_i'][mask]
                    mask[mask] &= mag < 24.5

    return mask


def sliding_window_def(snr, res, n_gal=100*100, n_bin_1d=20, snr_overlap=1.4):
    """
    A routine to define 2d bins in a sample based on its (snr, res)
    using a sliding window.  `n_gal` sets the minimum number of
    galaxies in a 2d bin, while `n_bin_1d` gives the number of bins to
    use in 1 dimension.  It returns an array of lower and upper bin
    edges for the overlapping SNR bins, each of length n_bin_1d.  It
    also returns resolution bin edges within each SNR bin, i.e., the
    lower and upper resolution bin edges have shape (n_bin_1d,
    n_bin_1d).
    """
    # First define bin edges in SNR.  To do this, we find percentiles
    # for `n_bin_1d` bins with overlap defined by `snr_overlap`.
    lower_percentiles = np.linspace(0.,100.,n_bin_1d+1)[:-1]
    delta_perc = lower_percentiles[1]-lower_percentiles[0]
    lower_percentiles -= (0.5*(snr_overlap-1.))*delta_perc
    lower_percentiles[0] = 0
    upper_percentiles = np.linspace(0.,100.,n_bin_1d+1)[1:]
    upper_percentiles += (0.5*(snr_overlap-1.))*delta_perc
    upper_percentiles[-1] = 100.
    snr_lower_edges = np.percentile(snr, lower_percentiles)
    snr_upper_edges = np.percentile(snr, upper_percentiles)
    res_lower_edges = np.zeros((n_bin_1d, n_bin_1d))
    res_upper_edges = np.zeros((n_bin_1d, n_bin_1d))

    # Now loop over the bins in SNR.
    for snr_bin in range(n_bin_1d):
        snr_bin_mask = (snr >= snr_lower_edges[snr_bin]) & \
            (snr < snr_upper_edges[snr_bin])
        tmp_res = res[snr_bin_mask]
        if len(tmp_res) < 1.1*n_gal:
            raise RuntimeError("Not enough objects to make "+
                               "sliding window")
        delta_perc = 100*float(n_gal) / len(tmp_res) # width of bin in percentiles
        lower_percentiles = np.linspace(0.,100.-delta_perc,n_bin_1d)
        upper_percentiles = lower_percentiles + delta_perc
        res_lower_edges[snr_bin, :] = np.percentile(tmp_res,
                                                    lower_percentiles)
        res_upper_edges[snr_bin, :] = np.percentile(tmp_res,
                                                    upper_percentiles)

    return snr_lower_edges, snr_upper_edges, res_lower_edges, res_upper_edges



def get_snr(catalog):
    """
    This utility computes the S/N for each object in the catalog, either based on
    magnitude error for the real data or on cmodel_flux for the simulations.
    It does not impose any cuts and returns NaNs for invalid S/N values.
    In the case of GREAT3-like simulations, the noise rescaling factor is applied
    to match the data.
    """
    snr     =   catalog['modelfit_CModel_flux']/catalog['modelfit_CModel_fluxSigma']
    return snr


def get_imag(catalog):
    """
    This utility returns the i-band magnitude of the objects in the input data or
    simulation catalog. Does not apply any cuts and returns NaNs for invalid
    values.
    """
    mag     =   catalog['Mag']
    #return -2.5 * np.log10(catalog['modelfit_CModel_flux']) + 27.
    return mag


def get_abs_ellip(catalog):
    """
    Returns the modulus of galaxy distortions.
    """
    absE    =   catalog['ext_shapeHSM_HsmShapeRegauss_e1']**2.+catalog['ext_shapeHSM_HsmShapeRegauss_e2']**2.
    return absE 


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

def get_true_shear(catalog):
    """
    This routine accesses the truth tables to get the true shear.
    """
    return catalog['g1'], catalog['g2']
    
def get_psf_size(catalog, type='trace'):
    """
    This utility gets the PSF size from a data or sims catalog using the
    specified size definition from the second moments matrix.

    Parameters
    ----------

    catalog: recarray
        Simulation or data catalog

    type: string
        Type of psf size measurement in ['trace', 'det', 'fwhm-arcsec'], default: 'trace'

    """
    # Convert into pixel units for the data
    psf_mxx = catalog['base_SdssShape_psf_xx']
    psf_myy = catalog['base_SdssShape_psf_yy']
    psf_mxy = catalog['base_SdssShape_psf_xy']

    if type == 'trace':
        size = psf_mxx + psf_myy
    elif type == 'det':
        size = (psf_mxx * psf_myy - psf_mxy**2)**(0.25)
    elif type == 'fwhm-arcsec':
        size = 0.168 * 2.355 * (psf_mxx * psf_myy - psf_mxy**2)**(0.25)
    else:
        raise ValueError("Unknown PSF size type: %s"%type)
    return size


def get_sigma_e_model(catalog):
    """
    This utility returns a model for the shape measurement uncertainty as a function of SNR and
    resolution.  It uses the catalog directly to get the SNR and resolution values.

    The file storing the data used to build the approximate correction is expected to be found in
    plot/sigmae_ratio.dat

    """
    # Get the necessary quantities.
    snr = get_snr(catalog)
    log_snr = np.log10(snr)
    res = get_res(catalog)
    # Fix interpolation for NaN/inf
    m = np.isnan(log_snr) | np.isinf(log_snr)
    log_snr[m] = 1. # set to minimum value that passes nominal cut
    snr[m] = 10.

    # Build the baseline model for sigma_e.
    #sigma_e = 0.264 * ((snr/20.)**-0.891) * ((res/0.5)**-1.015)
    sigma_e = 0.269 * ((snr/20.)**-0.956) * ((res/0.5)**-1.085)
    
    # Get the corrections from interpolation amongst saved values.
    data_file = "plot/sigmae/sigmae_ratio.dat"
    dat = np.loadtxt(data_file).transpose()
    saved_snr = dat[0,:]
    log_saved_snr = np.log10(saved_snr)
    saved_res = dat[1,:]
    saved_corr = dat[2,:]

    # Interpolate the corrections (which multiply the power-law results).
    result = grid_interpolate(log_saved_snr, saved_res, saved_corr, log_snr, res)
    return result*sigma_e 

def get_erms_model(catalog):
    """
    This utility returns a model for the RMS ellipticity as a function of SNR and
    resolution.  It uses the catalog directly to get the SNR and resolution values.

    The file storing the data used to build the model is expected to be found in
    plot/eRMS/intrinsicshape_2d.dat

    """
    # Get the necessary quantities.
    snr = get_snr(catalog)
    log_snr = np.log10(snr)
    res = get_res(catalog)
    # Fix interpolation for NaN/inf
    m = np.isnan(log_snr) | np.isinf(log_snr)
    log_snr[m] = 1. # set to minimum value that passes nominal cut
    snr[m] = 10.

    # Get saved model values.
    data_file = "plot/eRMS/intrinsicshape_2d.dat"
    dat = np.loadtxt(data_file).transpose()
    saved_snr = dat[0,:]
    log_saved_snr = np.log10(saved_snr)
    saved_res = dat[1,:]
    saved_model = dat[2,:]

    # Interpolate the e_rms values and return them.
    result = grid_interpolate(log_saved_snr, saved_res, saved_model, log_snr, res)
    return result


def get_weight_model(catalog):
    """
    This utility returns a model for the shape measurement weight as a
    function of SNR and resolution.  It relies on two other routines
    to get models for the intrinsic shape RMS and measurement error.
    """
    sigmae_meas = get_sigma_e_model(catalog)
    erms = get_erms_model(catalog)
    return 1./(sigmae_meas**2 + erms**2)


def grid_interpolate(x, y, z, eval_x, eval_y):
    """
    This is a utility for interpolating a 2D function z(x, y) linearly to values
    (x, y) = (eval_x, eval_y), but also enabling extrapolation beyond the (x, y)
    bounds using the nearest neighbor method.
    """
    from scipy.interpolate import griddata
    result = griddata((x, y), z, (eval_x, eval_y), method='linear')
    nn_result = griddata((x, y), z, (eval_x, eval_y), method='nearest')
    mask = np.isnan(result)
    result[mask] = nn_result[mask]
    return result
