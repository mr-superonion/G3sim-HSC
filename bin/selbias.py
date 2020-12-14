import numpy as np

def get_psf_ellip(catalog):
    """
    This utility gets the PSF ellipticity (distortion) from a data or sims
    catalog.  It does not impose flag cuts, but rather assumes that has already
    been done

    Parameters:
    ----------
    catalog: array_like
        input galaxy catalog

    Returns:
    e1_psf: array_like
        the first component of ellipticity
    e2_psf: array_like
        the second component of ellipticity
    -------
    """
    if 'i_sdssshape_psf_shape11' in catalog.dtype.names:
        psf_mxx = catalog['i_sdssshape_psf_shape11']
        psf_myy = catalog['i_sdssshape_psf_shape22']
        psf_mxy = catalog['i_sdssshape_psf_shape12']
    else:
        raise ValueError("Do not have PSF shape measurement columns")

    return (psf_mxx - psf_myy) / (psf_mxx + psf_myy), 2. * \
        psf_mxy / (psf_mxx + psf_myy)

def get_bias(weight,magA10,res):
    """
    This utility gets the selection bias (multiplicative and additive)

    Parameters:
    ----------
    weight: array_like
        Weight for dataset.  E.g., lensing shape weight, Sigma_c^-2 weight, ...
    res: array_like
        Resolution factor for dataset
    magA10: array_like
        aperture magnitude (1 arcsec) for dataset

    Returns:
    -------
    m_sel (float) :
        Multiplicative edge-selection bias
    a_sel (float) :
        additive edge-selection bias (c1)
    m_sel_err (float) :
        1-sigma uncertainty in m_sel
    a_sel_err (float) :
        1-sigma uncertainty in a_sel
    """

    if not(np.all(np.isfinite(weight))):
        raise ValueError("Non-finite weight")
    if not(np.all(weight) >= 0.0):
        raise ValueError("Negative weight")
    wSum    =   np.sum(weight)

    bin_magA=   0.025
    # magA10  =   get_imag_A10(data)
    pedgeM  =   np.sum(weight[(magA10>= 25.5-bin_magA)])/wSum/bin_magA

    bin_res =   0.01
    # res     =   get_ires(data)
    pedgeR  =   np.sum(weight[(res<= 0.3+bin_res)])/wSum/bin_res
    m_sel   =   -0.059*pedgeM+0.019*pedgeR
    a_sel   =   0.0064*pedgeM+0.0063*pedgeR

    # assume the errors for 2 cuts are independent.
    m_err   =   np.sqrt((0.0089*pedgeM)**2.+(0.0013*pedgeR)**2.)
    a_err   =   np.sqrt((0.0034*pedgeM)**2.+(0.0009*pedgeR)**2.)
    return m_sel,a_sel,m_err,a_err
