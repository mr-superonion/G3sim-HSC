import numpy as np

def get_res(catalog):
    """
    Returns the resolution
    """
    if 'res' in catalog.dtype.names:
        return catalog['res']
    elif 'i_hsmshaperegauss_resolution' in catalog.dtype.names:   #s18 & s19
        resolution  =   catalog['i_hsmshaperegauss_resolution']
    elif 'ishape_hsm_regauss_resolution' in catalog.dtype.names:  #s15
        resolution  =   catalog['ishape_hsm_regauss_resolution']
    elif 'ext_shapeHSM_HsmShapeRegauss_resolution' in catalog.dtype.names:#pipe 7
        resolution  =   catalog['ext_shapeHSM_HsmShapeRegauss_resolution']
    else:
        resolution  =   np.empty(len(catalog))
        resolution.fill(np.nan)
    return resolution

def get_imag_A10(catalog):
    """
    This utility returns the i-band magnitude of the objects in the input data or
    simulation catalog. Does not apply any cuts and returns NaNs for invalid
    values.
    """
    if 'magA10' in catalog.dtype.names:
        return catalog['magA10']
    elif 'i_apertureflux_10_mag' in catalog.dtype.names: #s18 s19
        magnitude   =   catalog['i_apertureflux_10_mag']
    elif 'base_CircularApertureFlux_3_0_instFlux' in catalog.dtype.names: #pipe 7
        magnitude   =   -2.5 * np.log10(catalog['base_CircularApertureFlux_3_0_instFlux']) + 27.
    else:
        magnitude   =   np.empty(len(catalog))
        magnitude.fill(np.nan)
    return magnitude
