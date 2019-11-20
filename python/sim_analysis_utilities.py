try:
    import astropy.io.fits as pyfits
except:
    import pyfits

from astropy.table import Table, join,vstack
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
field_keys = {'W02': 'XMM',
              'W03': 'GAMA09H',
              'W04-12': 'WIDE12H',
              'W04-15': 'GAMA15H',
              'W05': 'VVDS',
              'W06': 'HECTOMAP'}


field_names = {'XMM': 'W02',
              'GAMA09H': 'W03',
              'WIDE12H': 'W04_WIDE12H',
              'GAMA15H': 'W04_GAMA15H',
              'VVDS': 'W05',
              'HECTOMAP': 'W06'}


# Hex keys for the color scheme for shear catalog paper
plot_colors = {
    'GAMA09H' : '#d73027',
    'GAMA15H' : '#fc8d59',
    'HECTOMAP' : '#fee090',
    'VVDS' : '#000000',
    'WIDE12H' : '#91bfdb',
    'XMM' : '#4575b4'
    }

def m_func(x, b, c, d, e):
    """
    Empirically-motivated model we are trying to fit for m(SNR, res).

    x[0,:] is SNR, x[1,:] is resolution.
    """
    model = (x[0,:]/20.)**d
    model *= b*((x[1,:]/0.5)**c)
    return model+e


def a_func(x, b, c, d):
    """
    Empirically-motivated model we are trying to fit for a(SNR, res).

    x[0,:] is SNR, x[1,:] is resolution.
    """
    model = b*(x[1,:]-c)
    model *= (x[0,:]/20.)**d
    return model

def prepare_slide_catalog(anaDir):
    outFname    =   os.path.join(anaDir,'slide_cat.fits')
    if os.path.exists(outFname):
        return Table.read(outFname)
    colnames    =   ['modelfit_CModel_flux', 'modelfit_CModel_fluxSigma','Mag',
        'i_hsmshaperegauss_e1','i_hsmshaperegauss_e2','weak_lensing_flag',
        'i_hsmshaperegauss_resolution','i_hsmshaperegauss_sigma',
        'base_SdssShape_psf_xx','base_SdssShape_psf_yy','base_SdssShape_psf_xy',
        'reweight_s16a','i_hsmshaperegauss_derived_shape_weight','paired']
    catList =   []
    for key in field_keys:
        fieldname=  field_keys[key]
        if fieldname==  'AEGIS':
            continue
        fname   =   os.path.join(anaDir,fieldname,'%s.fits'%(fieldname))
        catalog =   Table.read(fname)[colnames]
        catList.append(catalog)
    slide_cat   =   vstack(catList)
    slide_cat.write(outFname)
    return slide_cat


def prepare_catalog(filename,fieldname,isim=-1):
    """
    # Read the catalog data
    """
    catalog =   Table.read(filename)
    colNames=   catalog.colnames
 
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
            (~np.isnan(catalog['ext_shapeHSM_HsmShapeRegauss_resolution']))&\
            (catalog['deblend_nChild']  ==  0)
    catalog =   catalog[mask]
    if len(catalog)==0:
        return None
    catalog =   catalog[colNames]

    # Add columns with names matching the real data
    cols = []
    data_fields = {'ext_shapeHSM_HsmShapeRegauss_e1': 'i_hsmshaperegauss_e1',
                   'ext_shapeHSM_HsmShapeRegauss_e2': 'i_hsmshaperegauss_e2',
                   'ext_shapeHSM_HsmShapeRegauss_sigma': 'i_hsmshaperegauss_sigma',
                   'ext_shapeHSM_HsmShapeRegauss_resolution': 'i_hsmshaperegauss_resolution'}
    for i in data_fields:
        catalog[data_fields[i]] = catalog[i]

    # Add field name infomation
    catalog['field_name']   = np.string_(fieldname)
    if isim>=0:
        """
        # For simulation
        """
        #Reweight factor to be determined
        catalog['reweight_s16a']    =  0.

        # Add a gridIndex column instead of using the index column, for compatibility
        catalog['ipos']     =   (catalog['base_SdssCentroid_y']//64)*100 +(catalog['base_SdssCentroid_x']//64)
        catalog['ipos']     =   catalog['ipos'].astype(np.int)
        catalog['centDist'] =   ((catalog['base_SdssCentroid_y']%64-32)**2. +
                                (catalog['base_SdssCentroid_x']%64-32)**2.)
        catalog['centDist'] =   np.sqrt(catalog['centDist'])

        # First, keep only detections that are the closest to the grid point
        # Get sorted index by grid index and grid distance
        inds        =   np.lexsort([catalog['centDist'], catalog['ipos']])
        inds_unique =   np.unique(catalog['ipos'][inds], return_index=True)[1]
        catalog     =   catalog[inds[inds_unique]]
        catalog     =   catalog[(catalog['centDist']<5.)]

        """
        # Load information from psfPre
        """
        psfFname    =   'psfPre/psf-%05d.fits' % isim
        psfInfo     =   pyfits.getheader(psfFname)
        catalog['g1']=  psfInfo['g1']
        catalog['g2']=  psfInfo['g2']

        """
        # Load information from catPre
        """
        catFname    =   'catPre/catalog-%05d.fits' % isim
        catInfo     =   Table.read(catFname)
        # Join catalogs on grid index
        # Extinction information is contained in catInfo 
        catalog = join(catalog, catInfo, 'ipos')

        """
        ***
        # Identify pairs based on their input ipos values
        """
        n = 100*100 # Number of grid points
        mask_even = np.in1d(catalog['ipos'] % n, 2*np.arange(n/2).astype('int'), assume_unique=True)
        mask_odd  = np.in1d(catalog['ipos'] % n, 2*np.arange(n/2).astype('int')+1, assume_unique=True)

        # Find index of broken pairs
        broken_pairs = np.setxor1d(catalog['ipos'][mask_even] % n,
                                   catalog['ipos'][mask_odd] % n - 1)
        mask_pairs = ~np.in1d(catalog['ipos'] % n,
                                  np.concatenate([broken_pairs, broken_pairs+1]), assume_unique=True)
        catalog['paired'] = mask_pairs
        catalog.sort(['paired', 'ipos'])
        catalog['isim'] =   isim
    
    # Add column to make a WL flag for the simulation data:
    #     imag_cmodel - a_i < 24.5
    #     sigma_e > 0, < 0.4
    #     |e|<2
    #     resolution factor
    #     i-band SNR
    sig_e   =   get_sigma_e(catalog)
    e       =   np.sqrt(catalog['i_hsmshaperegauss_e1']**2 +
               catalog['i_hsmshaperegauss_e2']**2)

    catalog['weak_lensing_flag'] = ((get_imag(catalog)-catalog['a_i']) < 24.5) & \
        (e <= 2.) & \
        (catalog['i_hsmshaperegauss_resolution'] >= 0.3) & \
        (get_snr(catalog) >= 10.) &\
        (sig_e > 0.0) & \
        (sig_e < 0.4)

    # add a bunch of columns that get used for WL analysis - fake it for now!
    catalog['i_hsmshaperegauss_derived_sigma_e'] = catalog['i_hsmshaperegauss_sigma']
    catalog['i_hsmshaperegauss_derived_rms_e'] = 0.365

    wt =  1./(catalog['i_hsmshaperegauss_derived_rms_e']**2 + \
                 catalog['i_hsmshaperegauss_derived_sigma_e']**2)
    catalog['i_hsmshaperegauss_derived_shape_weight'] = wt
    # Unpleasant discovery: the lines below used to say "= 0" and it turns out that meant those
    # columns were ints, which led to all kinds of issues later on.  This fix is for the next time a
    # sims catalog gets turned into HDF5.
    catalog['i_hsmshaperegauss_derived_shear_bias_m'] = 0.
    catalog['i_hsmshaperegauss_derived_shear_bias_c1'] = 0.
    catalog['i_hsmshaperegauss_derived_shear_bias_c2'] = 0.
    return catalog

def galaxy_selector(catalog, min_snr=None, max_snr=None,
                    min_res=None, max_res=None,
                    max_sigma_e=None, max_e=None, max_mag=None,
                    lower_psf_size=0.,upper_psf_size=100.,
                    min_multiband_snr_data=None, n_pass_multiband=2,
                    force_90_pair=True, bright_star_mask=False,
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

    force_90_pair: bool (default: True)
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

    `force_90_pair` is used to decide whether we impose the selection
    in a way that maintains 90 degree rotated pairs or not. This option only
    """
    mask = np.ones(len(catalog)).astype(bool)
    # Tests if the catalog contains the `paired` column which allows us to
    # apply cuts and keep the noise cancellation.
    cols = catalog.dtype.names
    if force_90_pair and 'paired' not in cols:
        warnings.warn("Noise cancellation requested on real data")
        force_90_pair = False

    if force_90_pair:
        ind_pair = np.where(catalog['paired'])[0]
        ind_pair1 = ind_pair[0::2]
        ind_pair2 = ind_pair[1::2]
        # Masking elements not in a pair to begin with
        mask &= catalog['paired']

    # Basic sanity checks: none of the quantities that we want to use should
    # be NaN. We also impose the sigma_e and e cuts here, because even just
    # a few bad values can totally mess stuff up.
    mask &= (~np.isnan(catalog['i_hsmshaperegauss_e1'])) & \
            (~np.isnan(catalog['i_hsmshaperegauss_e2'])) & \
            (~np.isnan(catalog['i_hsmshaperegauss_sigma'])) & \
            (~np.isnan(catalog['i_hsmshaperegauss_resolution'])) & \
            (~np.isnan(get_snr(catalog))) & \
            (~np.isnan(get_imag(catalog)))
    # For real data, also check valid psf measurements
    if 'i_sdssshape_psf_shape11' in catalog.dtype.names:
        mask &= (~np.isnan(catalog['i_sdssshape_psf_shape11'])) & \
                (~np.isnan(catalog['i_sdssshape_psf_shape22'])) & \
                (~np.isnan(catalog['i_sdssshape_psf_shape22']))

    if max_sigma_e is not None:
        sigma_e = get_sigma_e(catalog[mask])
        mask[mask] &= sigma_e <= max_sigma_e

    if max_e is not None:
        e = np.sqrt(catalog['i_hsmshaperegauss_e1'][mask]**2 +
                    catalog['i_hsmshaperegauss_e2'][mask]**2)
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
    if force_90_pair:
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
        mask[mask] &= catalog['i_hsmshaperegauss_resolution'][mask] >= min_res

    if max_res is not None:
        mask[mask] &= catalog['i_hsmshaperegauss_resolution'][mask] <= max_res

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
        #mask[mask] &= catalog['iblendedness_abs_flux'][mask] <= max_blendedness
        mask[mask] &= catalog['base_Blendedness_abs_flux'][mask] <= max_blendedness

    # For sims only.
    if background_cut is not None and 'background' in catalog.dtype.names:
        mask[mask] &= catalog['background'][mask] > background_cut

    #For PSF size 
    psf_size = get_psf_size(catalog[mask])
    if lower_psf_size > 0.:
        cut_val = np.percentile(psf_size, lower_psf_size)
        mask[mask] &= (psf_size >= cut_val)
    if upper_psf_size < 100.:
        cut_val = np.percentile(psf_size, upper_psf_size)
        mask[mask] &= (psf_size < cut_val)

    # For noise cancellation, apply the cut only to the first element of the
    # pair to avoid funny selection biases.
    # But enable introduction of selection bias due to single cuts through applyboth.
    if force_90_pair:
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


def weighted_percentile(data,percents,weights=None):
    ''' percents in units of 1%
        weights specifies the frequency (count) of data.
    '''
    if weights is None:
        return np.percentile(data, percents)
    ind=np.argsort(data)
    d=data[ind]
    w=weights[ind]
    p=1.*w.cumsum()/w.sum()*100
    y=np.interp(percents, p, d)
    return y

def sliding_window_def(snr,res,weights=None,n_gal=100*100, n_bin_1d=20, snr_overlap=1.4):
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
    n_per_bin = len(snr) / (n_bin_1d**2)
    if n_per_bin < n_gal:
        # If we can't divide up into that many bins, then use a sliding window.
        # First define bin edges in SNR.  To do this, we find percentiles
        # for `n_bin_1d` bins with overlap defined by `snr_overlap`.
        lower_percentiles = np.linspace(0.,100.,n_bin_1d+1)[:-1]
        delta_perc = lower_percentiles[1]-lower_percentiles[0]
        lower_percentiles -= (0.5*(snr_overlap-1.))*delta_perc
        lower_percentiles[0] = 0
        upper_percentiles = np.linspace(0.,100.,n_bin_1d+1)[1:]
        upper_percentiles += (0.5*(snr_overlap-1.))*delta_perc
        upper_percentiles[-1] = 100.
        snr_lower = weighted_percentile(snr,lower_percentiles,weights)
        snr_upper = weighted_percentile(snr,upper_percentiles,weights)
        res_lower = np.zeros((n_bin_1d, n_bin_1d))
        res_upper = np.zeros((n_bin_1d, n_bin_1d))
        # Now loop over the bins in SNR.
        for snr_bin in range(n_bin_1d):
            snr_bin_mask = (snr >= snr_lower[snr_bin]) & \
                (snr < snr_upper[snr_bin])
            tmp_res = res[snr_bin_mask]
            if weights  is not None:
                tmp_wt  =   weights[snr_bin_mask]
            else:
                tmp_wt  =   None
            if len(tmp_res) < 1.1*n_gal:
                raise RuntimeError("Not enough objects to make "+
                                   "sliding window")
            delta_perc = 100*float(n_gal) / len(tmp_res) # width of bin in percentiles
            lower_percentiles = np.linspace(0.,100.-delta_perc,n_bin_1d)
            upper_percentiles = lower_percentiles + delta_perc
            res_lower[snr_bin, :] = weighted_percentile(tmp_res,
                                                lower_percentiles,tmp_wt)
            res_upper[snr_bin, :] = weighted_percentile(tmp_res,
                                                upper_percentiles,tmp_wt)

    else:
        # Otherwise, just divide them up.
        lower_percentiles = np.linspace(0.,100.,n_bin_1d+1)[:-1]
        upper_percentiles = np.linspace(0.,100.,n_bin_1d+1)[1:]
        snr_lower = np.percentile(snr, lower_percentiles)
        snr_upper = np.percentile(snr, upper_percentiles)
        res_lower = np.zeros((n_bin_1d,n_bin_1d))
        res_upper = np.zeros((n_bin_1d,n_bin_1d))
        for ind in range(n_bin_1d):
            m = (snr >= snr_lower[ind]) & (snr < snr_upper[ind])
            tmp_res = res[m]
            res_lower[ind,:] = np.percentile(tmp_res, lower_percentiles)
            res_upper[ind,:] = np.percentile(tmp_res, upper_percentiles)
            
    return snr_lower, snr_upper, res_lower, res_upper

def estimate_subfield_shear(catalog, use_model=False, verbose=False,
                            simple_resp=True, report_resp=False, reweight_s16a=False,
                            roundtrip_data=None,force_weight = True):
    """
    The purpose of this routine is to estimate the shear for a given subfield.

    This function expects the true shear to be the same for all entries.  If it
    finds >1 subfield index it will raise a RuntimeError.

    Parameters
    ----------

    catalog: recarray
        Catalog to use for shear estimation.

    use_model: bool (default: False)
        Use get_X_model routines for weights, RMS ellipticity, and sigma_e?  If
        False, then will use catalog entries.

    verbose: bool (default: False)
        Emit diagnostic statements?


    simple_resp: bool (default: True)
        Use the simpler formula for the responsivity?

    report_resp: bool (default: False)
        Return the value of responsivity along with the shears?

    reweight_s16a: bool (default: False)
        Use the S16A reweighting?  This will be the same value for each entry
        since the routine works at the subfield level.

    roundtrip_data: tuple or none (default: None)
        Data to be used for a round-trip analysis given some set of calibration
        corrections that differs from the ones in the catalogs.  This should be
        a tuple of (m, c1, c2) values with the same length as `catalog`.  If it
        is passed in, then the calibration values in the catalog are completely
        ignored.  If None, then the catalog calibration values are used.
    """

    if verbose:
        print('Using %d objects in input catalog' % (len(catalog)))

    if roundtrip_data is not None:
        model_m, model_c1, model_c2 = roundtrip_data

    # The calculation below follows the agreed-upon procedure from
    # http://hscsurvey.pbworks.com/w/page/113068663/notes%20on%20lensing%20signal%20calculation%20with%20S16A%20catalog
    if not use_model:
        weights = catalog['i_hsmshaperegauss_derived_shape_weight']
        e_rms = catalog['i_hsmshaperegauss_derived_rms_e']
        if 'i_hsmshaperegauss_derived_sigma_e' in catalog.dtype.names:
            sigma_e = catalog['i_hsmshaperegauss_derived_sigma_e']
        else:
            sigma_e = get_sigma_e(catalog)
    else:
        weights = get_weight_model(catalog)
        e_rms = get_erms_model(catalog)
        sigma_e = get_sigma_e_model(catalog)
    if reweight_s16a:
        weights *= catalog['reweight_s16a']
    # Enforce equality of weights etc. for galaxies in a pair if shape noise cancellation is requested.
    if force_weight:
        weights[1::2] = weights[0::2]
        e_rms[1::2] = e_rms[0::2]
        sigma_e[1::2] = sigma_e[0::2]

    # Deal gracefully with the case of zero weights (generally due to reweight_s16a option).
    if np.sum(weights)==0:
        if not report_resp:
            return 0, 0
        else:
            return 0, 0, 0

    e1 = catalog['i_hsmshaperegauss_e1']
    e2 = catalog['i_hsmshaperegauss_e2']
    shear_num_1 = np.sum(weights * e1)
    shear_num_2 = np.sum(weights * e2)
    sum_weight = np.sum(weights)
    if roundtrip_data is None:
        bias_num =  0. 
        add_num_1 = 0. 
        add_num_2 = 0. 
    else:
        bias_num = np.sum(weights*model_m)
        add_num_1 = np.sum(weights*model_c1)
        add_num_2 = np.sum(weights*model_c2)
    bias = bias_num / sum_weight

    if simple_resp:
        resp_num = np.sum(weights * e_rms**2)
        responsivity = 1.0 - resp_num/sum_weight
    else:
        f = e_rms**2 * weights # e_rms^2 / (e_rms^2 + sigma_e^2)
        k0 = (1.0-f) * e_rms**2
        k1 = f**2
        resp_num = np.sum(weights * (1.0 - k0 - 0.5*k1*(e1**2 + e2**2)))
        responsivity = resp_num / sum_weight


    if responsivity > 2. or responsivity < 0.:
        raise RuntimeError(
            "Error: responsivity is %f with sum_weight %f"%(responsivity,sum_weight))
    if verbose:
        print("Shear responsivity: %f" % responsivity)
    if verbose:
        print("Bias: %f, %f, %f" %(bias,bias_num,sum_weight))
    g1 = shear_num_1 / (2.0*responsivity*sum_weight*(1.0+bias)) \
        - add_num_1/((1.0+bias)*sum_weight)
    g2 = shear_num_2 / (2.0*responsivity*sum_weight*(1.0+bias)) \
        - add_num_2/((1.0+bias)*sum_weight)
    if not report_resp:
        return g1, g2
    else:
        return g1, g2, responsivity



def getMeanAbsEllip(catalog, weight=None):
    """Routine for estimating the mean absolute ellipticity from a catalog.
    Should be used on catalogs that have already been suitably masked so
    they don't have NaN values or other garbage. Can optionally use weights.

    Parameters
    ----------

    catalog: recarray
        Catalog to process

    weight: ndarray (optional)
        Weight array of the same size as the catalog

    Returns
    -------

    mean_e: float
        Mean absolute ellipticity
    """
    # If there are no weights given, then simply use a set of weights that are
    # equal to 1.
    if weight is None:
        weight = np.ones(len(catalog))

    mean_e = np.mean(weight * np.sqrt(catalog['i_hsmshaperegauss_e1']**2 +
                                      catalog['i_hsmshaperegauss_e2']**2))

    return mean_e


def getMeanRes(catalog, weight=None):
    """
    Routine for estimating the mean resolution from a catalog.
    Should be used on catalogs that have already been suitably masked so
    they don't have NaN values or other garbage. Can optionally use weights.

    Parameters
    ----------

    catalog: recarray
        Catalog to process

    weight: ndarray (optional)
        Weight array of the same size as the catalog

    Returns
    -------

    mean_res: float
        Mean resolution
    """
    if weight is None:
        weight = np.ones(len(catalog))

    return np.mean(weight * catalog['i_hsmshaperegauss_resolution'])


def getMeanSnr(catalog, weight=None):
    """
    Routine for estimating the mean S/N from a catalog.
    Should be used on catalogs that have already been suitably masked so
    they don't have NaN values or other garbage. Can optionally use weights.

    Parameters
    ----------

    catalog: recarray
        Catalog to process

    weight: ndarray (optional)
        Weight array of the same size as the catalog

    Returns
    -------

    mean_snr: float
        Mean snr
    """
    if weight is None:
        weight = np.ones(len(catalog))

    tmp_snr = weight * get_snr(catalog)

    return np.mean(tmp_snr)


def fitline(xarr, yarr):
    """Fit a line y = a + b * x to input x and y arrays by least squares.

    Returns the tuple (a, b, Var(a), Cov(a, b), Var(b)), after performing an internal estimate of
    measurement errors from the best-fitting model residuals.

    See Numerical Recipes (Press et al 1992; Section 15.2) for a clear description of the details
    of this simple regression.
    """
    # Get the S values (use default sigma2, best fit a and b still valid for
    # stationary data)
    S, Sx, Sy, Sxx, Sxy = _calculateSvalues(xarr, yarr)
    # Get the best fit a and b
    Del = S * Sxx - Sx * Sx
    a = (Sxx * Sy - Sx * Sxy) / Del
    b = (S * Sxy - Sx * Sy) / Del
    # Use these to estimate the sigma^2 by residuals from the best-fitting
    # model
    ymodel = a + b * xarr
    sigma2 = np.mean((yarr - ymodel)**2)
    # And use this to get model parameter error estimates
    var_a = sigma2 * Sxx / Del
    cov_ab = - sigma2 * Sx / Del
    var_b = sigma2 * S / Del
    return (a, b, var_a, cov_ab, var_b)


def _calculateSvalues(xarr, yarr, sigma2=1.):
    """Calculates the intermediate S values required for basic linear regression.

    See, e.g., Numerical Recipes (Press et al 1992) Section 15.2.
    """
    if len(xarr) != len(yarr):
        raise ValueError("Input xarr and yarr differ in length!")
    if len(xarr) <= 1:
        raise ValueError("Input arrays must have 2 or more values elements.")

    S = len(xarr) / sigma2
    Sx = np.sum(xarr / sigma2)
    Sy = np.sum(yarr / sigma2)
    Sxx = np.sum(xarr * xarr / sigma2)
    Sxy = np.sum(xarr * yarr / sigma2)
    return (S, Sx, Sy, Sxx, Sxy)


def show_shear_bias(g, g_true, m=None, merr=None, c=None):
    """
    This is a routine to plot estimated vs. true shear values (for one component).  Optionally, a
    trend-line showing some calibration bias value can be shown.
    """
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.scatter(g_true, g - g_true)
    if m is not None:
        x_vals = np.sort(g_true)
        y_vals = m * x_vals
        if c is not None:
            y_vals += c
        ax.plot(x_vals, y_vals, 'r-')
        if merr is not None:
            ax.set_title('m=%.3f +/- %.3f'%(m,merr))
    ax.plot(x_vals, np.zeros_like(x_vals), 'k--')
    plt.xlabel('true g')
    plt.ylabel('g - true g')
    return fig


def show_shear_additive_bias(residual, e_psf, a, aerr):
    """
    This is a routine for showing additive bias.  It takes a NumPy array of the residuals,
    g-(1+m)g_true, and plots it against the PSF ellipticity.  The slope should be the "a" value.
    """
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.scatter(e_psf, residual)
    x_vals = np.sort(e_psf)
    y_vals = a*x_vals
    ax.set_title('a=%.3f +/- %.3f'%(a,aerr))
    ax.plot(x_vals, y_vals, 'r-')
    ax.plot(x_vals, np.zeros_like(x_vals), 'k--')
    plt.xlabel('PSF e (this component)')
    plt.ylabel('Residual additive error')
    return fig


def get_true_shear(catalog):
    """
    This routine accesses the truth tables to get the true shear.
    """
    return catalog['g1'], catalog['g2']


def systematics_func(x, m, a):
    return m * x[0, :] + a * x[1, :]


def get_m_c_a(g1, g2, g1_true, g2_true, psf_e1, psf_e2, sigma_vals, separate_components=True):
    """
    This routine takes a list/tuple/NumPy array of values for estimated shear (`g1`, `g2`), true
    shear (`g1_true`, `g2_true`), and PSF ellipticity (`psf_e1`, `psf_e2`).  Then it does a
    GREAT3-like analysis for calibration bias and additive bias.

    Depending on the value of `separate_components`, it does this
    either separately for each component, or for both combined.
    """
    # First, just make sure the inputs are NumPy arrays to enable the math that has to happen later.
    g1 = np.array(g1)
    g2 = np.array(g2)
    g1_true = np.array(g1_true)
    g2_true = np.array(g2_true)
    psf_e1 = np.array(psf_e1)
    psf_e2 = np.array(psf_e2)
    sigma_vals = np.array(sigma_vals)

    if not separate_components:
        g = np.concatenate([g1, g2])
        g_true = np.concatenate([g1_true, g2_true])
        psf_e = np.concatenate([psf_e1, psf_e2])
        sigma_vals = np.concatenate([sigma_vals, sigma_vals])

    # identify which ones we don't want to plot
    to_use = (sigma_vals < 1.e10)

    # Our default algorithm is going to be to fit for
    # g-g_true = m g_true + a e_PSF
    # for each component, using scipy.optimize.curve_fit.
    # This can eventually accomodate weights as needed.
    import scipy.optimize
    if separate_components:
        popt1, pcov1 = scipy.optimize.curve_fit(systematics_func,
                                                np.vstack((g1_true,psf_e1)), 
                                                g1-g1_true,sigma=sigma_vals)
        popt2, pcov2 = scipy.optimize.curve_fit(systematics_func,
                                                np.vstack((g2_true,psf_e2)),
                                                g2-g2_true,sigma=sigma_vals)
        # Return everything: (m, a, sigma_m, sigma_a) first for component 1 then for
        # component 2.
        return popt1[0], popt1[1], math.sqrt(pcov1[0][0]), math.sqrt(pcov1[1][1]), \
            popt2[0], popt2[1], math.sqrt(pcov2[0][0]), math.sqrt(pcov2[1][1])
    else:
        popt, pcov = scipy.optimize.curve_fit(systematics_func,
                                              np.vstack((g_true,psf_e)),
                                              g-g_true,sigma=sigma_vals)
        # Return everything: (m, a, sigma_m, sigma_a).
        return popt[0], popt[1], math.sqrt(pcov[0][0]), math.sqrt(pcov[1][1]) 

def get_snr(catalog):
    """
    This utility computes the S/N for each object in the catalog, either based on
    magnitude error for the real data or on cmodel_flux for the simulations.
    It does not impose any cuts and returns NaNs for invalid S/N values.
    In the case of GREAT3-like simulations, the noise rescaling factor is applied
    to match the data.
    """
    if 'i_cmodel_flux' in catalog.dtype.names:
        snr     =   catalog['i_cmodel_flux']/catalog['i_cmodel_fluxsigma']
    else:
        snr     =   catalog['modelfit_CModel_flux']/catalog['modelfit_CModel_fluxSigma']
    return snr

def get_imag(catalog):
    """
    This utility returns the i-band magnitude of the objects in the input data or
    simulation catalog. Does not apply any cuts and returns NaNs for invalid
    values.
    """
    if 'i_cmodel_mag' in catalog.dtype.names:
        magnitude   =   catalog['i_cmodel_mag']
    else:
        magnitude   =   -2.5 * np.log10(catalog['modelfit_CModel_flux']) + 27.
    return magnitude

def get_abs_ellip(catalog):
    """
    Returns the modulus of galaxy distortions.
    """
    if 'i_hsmshaperegauss_e1' in catalog.dtype.names:
        absE    =   catalog['i_hsmshaperegauss_e1']**2.+catalog['i_hsmshaperegauss_e2']**2.
    else:
        absE    =   catalog['ext_shapeHSM_HsmShapeRegauss_e1']**2.+catalog['ext_shapeHSM_HsmShapeRegauss_e2']**2.
    absE    =   np.sqrt(absE)
    return absE 

def get_res(catalog):
    """
    Returns the resolution
    """
    if 'i_hsmshaperegauss_resolution' in catalog.dtype.names:
        resolution  =   catalog['i_hsmshaperegauss_resolution']
    else:
        resolution  =   catalog['ext_shapeHSM_HsmShapeRegauss_resolution']
    return resolution 

def get_sigma_e(catalog):
    """
    This utility returns the hsm_regauss_sigma values for the catalog, without
    imposing any additional flag cuts.
    In the case of GREAT3-like simulations, the noise rescaling factor is applied
    to match the data.
    """
    sigma_e = catalog['i_hsmshaperegauss_sigma']
    return sigma_e

def get_true_shear(catalog):
    """
    This routine accesses the truth tables to get the true shear.
    """
    return catalog['g1'], catalog['g2']

def get_psf_size(catalog, type='fwhm-arcsec'):
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
    if 'base_SdssShape_psf_xx' in catalog.dtype.names:
        psf_mxx = catalog['base_SdssShape_psf_xx']*0.168**2.
        psf_myy = catalog['base_SdssShape_psf_yy']*0.168**2.
        psf_mxy = catalog['base_SdssShape_psf_xy']*0.168**2.
    else:
        psf_mxx = catalog['i_sdssshape_psf_shape11']
        psf_myy = catalog['i_sdssshape_psf_shape22']
        psf_mxy = catalog['i_sdssshape_psf_shape12']

    if type == 'trace':
        size = psf_mxx + psf_myy
    elif type == 'fwhm-arcsec':
        size = 2.355 * (psf_mxx * psf_myy - psf_mxy**2)**(0.25)
    else:
        raise ValueError("Unknown PSF size type: %s"%type)
    return size

def get_psf_ellip(catalog):
    """
    This utility gets the PSF ellipticity (distortion) from a data or
    sims catalog.
    It does not impose flag cuts, but rather assumes that has already been done
    """
    if 'base_SdssShape_psf_xx' in catalog.dtype.names:
        psf_mxx = catalog['base_SdssShape_psf_xx']*0.168
        psf_myy = catalog['base_SdssShape_psf_yy']*0.168
        psf_mxy = catalog['base_SdssShape_psf_xy']*0.168
    else:
        psf_mxx = catalog['i_sdssshape_psf_shape11']
        psf_myy = catalog['i_sdssshape_psf_shape22']
        psf_mxy = catalog['i_sdssshape_psf_shape12']
    return (psf_mxx - psf_myy) / (psf_mxx + psf_myy), 2. * \
        psf_mxy / (psf_mxx + psf_myy)


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
    sigma_e = 0.271 * ((snr/20.)**-0.935) * ((res/0.5)**-0.978)
    
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

def get_m_model(catalog, weight_bias=True):
    """
    Routine to get a model for calibration bias m given some input snr
    and resolution values or arrays.
    """
    # Get the necessary quantities.
    snr = get_snr(catalog)
    log_snr = np.log10(snr)
    res = get_res(catalog)
    # Fix interpolation for NaN/inf
    m = np.isnan(log_snr) | np.isinf(log_snr)
    log_snr[m] = 1. # set to minimum value that passes nominal cut
    snr[m] = 10.

    #m = -0.1408*((snr/20.)**-1.23)*((res/0.5)**1.76) - 0.0214
    maFname =   './plot/reGausCalib/shear_m_a_model_par.npy'
    m_opt   =   np.load(maFname).item()['m_opt']
    fake_x  =   np.vstack((snr,res))
    model_m =   m_func(fake_x,*m_opt)
    
    data_file = "plot/reGausCalib/shear_dm_da_2d.csv"
    dat = Table.read(data_file)
    saved_snr = dat['snr']
    log_saved_snr = np.log10(saved_snr)
    saved_res = dat['res']
    saved_m = dat['dm']
    # Interpolate the model residuals and return them.  These are additional
    # biases beyond the power-law model and so should be added to the power-law.
    result = grid_interpolate(log_saved_snr, saved_res, saved_m, log_snr, res)

    if weight_bias:
        result += get_mwt_model(catalog)

    model_m   = model_m + result
    return model_m

def get_mwt_model(catalog):
    """
    Routine to get a model for calibration bias m due to weight bias, given some input snr
    and resolution values or arrays.
    """
    # Get the necessary quantities.
    snr = get_snr(catalog)
    log_snr = np.log10(snr)
    res = get_res(catalog)
    # Fix interpolation for NaN/inf
    m = np.isnan(log_snr) | np.isinf(log_snr)
    log_snr[m] = 1. # set to minimum value that passes nominal cut
    snr[m] = 10.

    #m = -1.31 + (27.26 + (snr/20.)**-1.22) / (res + 20.8)
    maFname =   './plot/reGausCalib/weightBias_m_a_model_par.npy'
    m_opt   =   np.load(maFname).item()['m_opt']
    fake_x  =   np.vstack((snr,res))
    model_m =   m_func(fake_x,*m_opt)
    
    data_file = "plot/reGausCalib/weightBias_dm_da_2d.csv"
    dat = Table.read(data_file)
    saved_snr = dat['snr']
    log_saved_snr = np.log10(saved_snr)
    saved_res = dat['res']
    saved_m = dat['dm']
    # Interpolate the model residuals and return them.  These are additional
    # biases beyond the power-law model and so should be added to the power-law.
    result = grid_interpolate(log_saved_snr, saved_res, saved_m, log_snr, res)
    return result + model_m

def get_c_model(catalog,weight_bias=True):
    """
    Routine to get a model for additive bias coefficient a given some input snr
    and resolution values or arrays.
    """
    # Get the necessary quantities.
    snr = get_snr(catalog)
    log_snr = np.log10(snr)
    res = get_res(catalog)
    # Fix interpolation for NaN/inf
    m = np.isnan(log_snr) | np.isinf(log_snr)
    
    
    log_snr[m] = 1. # set to minimum value that passes nominal cut
    snr[m] = 10.

    #a = 0.175 * ((snr/20.)**-1.07) * (res - 0.508)

    maFname =   './plot/reGausCalib/shear_m_a_model_par.npy'
    a_opt   =   np.load(maFname).item()['a_opt']
    fake_x  =   np.vstack((snr,res))
    model_a =   a_func(fake_x,*a_opt)
    
    data_file = "plot/reGausCalib/shear_dm_da_2d.csv"
    dat = Table.read(data_file)
    saved_snr = dat['snr']
    log_saved_snr = np.log10(saved_snr)
    saved_res = dat['res']
    saved_a = dat['da']
    # Interpolate the model residuals and return them.  These are additional
    # biases beyond the power-law model and so should be added to the power-law.
    result = grid_interpolate(log_saved_snr, saved_res, saved_a, log_snr, res)
    if weight_bias:
        result += get_awt_model(catalog)
    model_a =  model_a +  result
    psf_e1, psf_e2 = get_psf_ellip(catalog)
    model_c1 = model_a*psf_e1
    model_c2 = model_a*psf_e2
    return model_c1,model_c2 

def get_awt_model(catalog):
    """
    Routine to get a model for additive bias coefficient a due to weight bias given some input snr
    and resolution values or arrays.
    """
    # Get the necessary quantities.
    snr = get_snr(catalog)
    log_snr = np.log10(snr)
    res = get_res(catalog)
    # Fix interpolation for NaN/inf
    m = np.isnan(log_snr) | np.isinf(log_snr)
    log_snr[m] = 1. # set to minimum value that passes nominal cut
    snr[m] = 10.

    #a = -0.089 * (res-0.71) * ((snr/20.)**-2.2)
    maFname =   './plot/reGausCalib/weightBias_m_a_model_par.npy'
    a_opt   =   np.load(maFname).item()['a_opt']
    fake_x  =   np.vstack((snr,res))
    model_a =   a_func(fake_x,*a_opt)
    
    data_file = "plot/reGausCalib/weightBias_dm_da_2d.csv"
    dat = Table.read(data_file)
    saved_snr = dat['snr']
    log_saved_snr = np.log10(saved_snr)
    saved_res = dat['res']
    saved_a = dat['da']
    # Interpolate the model residuals and return them.  These are additional
    # biases beyond the power-law model and so should be added to the power-law.
    result = grid_interpolate(log_saved_snr, saved_res, saved_a, log_snr, res)
    return result + model_a

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

def update_wl_flags_sim(catalog):
    sig_e   =   get_sigma_e(catalog)
    e       =   np.sqrt(catalog['i_hsmshaperegauss_e1']**2 +
                catalog['i_hsmshaperegauss_e2']**2)
    catalog['weak_lensing_flag'] = ((get_imag(catalog)-catalog['a_i']) < 24.5) & \
        (e <= 2.) & \
        (catalog['i_hsmshaperegauss_resolution'] >= 0.3) & \
        (get_snr(catalog) >= 10.) &\
        (sig_e > 0.0) & \
        (sig_e < 0.4)
    return catalog 
