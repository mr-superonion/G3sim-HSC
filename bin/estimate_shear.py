import numpy as np
import sim_analysis_utilities as utilities

def estimate_pergal_shear(catalog):
    """
    Parameters
    ----------
    catalog: recarray
        Catalog to use for shear estimation.
    """
    # Weights
    if 'i_hsmshaperegauss_derived_weight' in catalog.dtype.names:
        weights =   catalog['i_hsmshaperegauss_derived_weight']#s19
    else:
        weights =   catalog['ishape_hsm_regauss_derived_shape_weight']#s16
    sum_weight  = np.sum(weights)
    # Responsivity
    if 'i_hsmshaperegauss_derived_rms_e' in catalog.dtype.names:
        e_rms   =   catalog['i_hsmshaperegauss_derived_rms_e']#s19
    else:
        e_rms   =   catalog['ishape_hsm_regauss_derived_rms_e']#s16
    resp_num = np.sum(weights * e_rms**2)
    responsivity = 1.0 - resp_num/sum_weight
    if responsivity > 2. or responsivity < 0.:
        raise RuntimeError(
            "Error: responsivity is %f with sum_weight %f"%(responsivity,sum_weight))
    # Multiplicative bias
    if 'i_hsmshaperegauss_derived_shear_bias_m' in catalog.dtype.names:
        model_m=catalog['i_hsmshaperegauss_derived_shear_bias_m']#s19
    else:
        model_m=catalog['ishape_hsm_regauss_derived_shear_bias_m']#s16
    bias_num= np.sum(weights*model_m)
    bias = bias_num / sum_weight

    # Shear
    e1,e2   =   utilities.get_gal_ellip(catalog)
    # Additive bias
    if 'i_hsmshaperegauss_derived_shear_bias_c1' in catalog.dtype.names:
        c1=catalog['i_hsmshaperegauss_derived_shear_bias_c1'] #s19
        c2=catalog['i_hsmshaperegauss_derived_shear_bias_c2']
    else:
        c1=catalog['ishape_hsm_regauss_derived_shear_bias_c1']#s16
        c2=catalog['ishape_hsm_regauss_derived_shear_bias_c2']

    g1 = e1 / (2.0*responsivity*(1.0+bias)) - c1/((1.0+bias))
    g2 = e2 / (2.0*responsivity*(1.0+bias)) - c2/((1.0+bias))
    return g1, g2, weights
