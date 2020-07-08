SELECT
  meas.object_id
, meas.parent_id
, meas.i_ra
, meas.i_dec
, meas.i_variance_value
, meas.i_footprintarea_value
, meas.i_deblend_parenttoobig

-- unforced CModel magnitudes
, meas.i_cmodel_mag
, meas.i_cmodel_magerr
, meas.i_cmodel_flag
, meas.i_cmodel_flag_badcentroid
, meas.i_cmodel_flux
, meas.i_cmodel_fluxerr
, meas.i_cmodel_objective

-- unforced Aperture fluxes
,meas3.i_apertureflux_10_flux
,meas3.i_apertureflux_10_fluxerr
,meas3.i_apertureflux_15_flux
,meas3.i_apertureflux_15_fluxerr
,meas3.i_apertureflux_20_flux
,meas3.i_apertureflux_20_fluxerr
,meas3.i_apertureflux_30_flux
,meas3.i_apertureflux_30_fluxerr
,meas3.i_apertureflux_40_flux
,meas3.i_apertureflux_40_fluxerr

-- forced measurement
-- extinction
, forced.a_g
, forced.a_r
, forced.a_i
, forced.a_z
, forced.a_y

-- forced CModel magnitudes and fluxes
, forced.g_cmodel_mag       as forced_g_cmodel_mag
, forced.g_cmodel_magerr    as forced_g_cmodel_magerr
, forced.g_cmodel_flux      as forced_g_cmodel_flux
, forced.g_cmodel_fluxerr   as forced_g_cmodel_fluxerr
, forced.g_cmodel_flag      as forced_g_cmodel_flag

, forced.r_cmodel_mag       as forced_r_cmodel_mag
, forced.r_cmodel_magerr    as forced_r_cmodel_magerr
, forced.r_cmodel_flux      as forced_r_cmodel_flux
, forced.r_cmodel_fluxerr   as forced_r_cmodel_fluxerr
, forced.r_cmodel_flag      as forced_r_cmodel_flag

, forced.i_cmodel_mag       as forced_i_cmodel_mag
, forced.i_cmodel_magerr    as forced_i_cmodel_magerr
, forced.i_cmodel_flux      as forced_i_cmodel_flux
, forced.i_cmodel_fluxerr   as forced_i_cmodel_fluxerr
, forced.i_cmodel_flag      as forced_i_cmodel_flag

, forced.z_cmodel_mag       as forced_z_cmodel_mag
, forced.z_cmodel_magerr    as forced_z_cmodel_magerr
, forced.z_cmodel_flux      as forced_z_cmodel_flux
, forced.z_cmodel_fluxerr   as forced_z_cmodel_fluxerr
, forced.z_cmodel_flag      as forced_z_cmodel_flag

, forced.y_cmodel_mag       as forced_y_cmodel_mag
, forced.y_cmodel_magerr    as forced_y_cmodel_magerr
, forced.y_cmodel_flux      as forced_y_cmodel_flux
, forced.y_cmodel_fluxerr   as forced_y_cmodel_fluxerr
, forced.y_cmodel_flag      as forced_y_cmodel_flag

,meas3.i_apertureflux_10_mag
,meas3.i_apertureflux_10_magerr
,meas3.i_apertureflux_15_mag
,meas3.i_apertureflux_15_magerr
,meas3.i_apertureflux_20_mag
,meas3.i_apertureflux_20_magerr
,meas3.i_apertureflux_30_mag
,meas3.i_apertureflux_30_magerr

,meas2.i_psfflux_mag
,meas2.i_psfflux_magerr

-- shapes
-- hsm regauss shapes
, meas2.i_hsmshaperegauss_e1
, meas2.i_hsmshaperegauss_e2
, meas2.i_hsmshaperegauss_sigma
, meas2.i_hsmshaperegauss_resolution
-- derived hsm regauss quantitites
-- weight assuming rms_e = 0.365
, 1./(meas2.i_hsmshaperegauss_sigma*meas2.i_hsmshaperegauss_sigma + 0.365*0.365) as i_hsmshaperegauss_derived_weight
, meas2.i_hsmshaperegauss_sigma     as i_hsmshaperegauss_derived_sigma_e
, CAST(0.365 as FLOAT)              as i_hsmshaperegauss_derived_rms_e
, CAST(0. as FLOAT)                 as i_hsmshaperegauss_derived_bias_m
, CAST(0. as FLOAT)                 as i_hsmshaperegauss_derived_bias_c1
, CAST(0. as FLOAT)                 as i_hsmshaperegauss_derived_bias_c2

-- sdss shapes (without PSF correction)
, meas2.i_sdssshape_shape11
, meas2.i_sdssshape_shape22
, meas2.i_sdssshape_shape12
, meas2.i_sdssshape_psf_shape11
, meas2.i_sdssshape_psf_shape22
, meas2.i_sdssshape_psf_shape12


-- columns which can be used for selection
, meas.tract
, meas.patch
, meas.merge_peak_g
, meas.merge_peak_r
, meas.merge_peak_r2
, meas.merge_peak_i
, meas.merge_peak_i2
, meas.merge_peak_z
, meas.merge_peak_y

-- Input exposures number
, meas.g_inputcount_value
, meas.r_inputcount_value
, meas.i_inputcount_value
, meas.z_inputcount_value
, meas.y_inputcount_value

-- Mask plane
, meas.i_pixelflags_inexact_psfcenter
, meas.i_pixelflags_inexact_psf
-- Bright object mask (Arcturus)
, masks.g_mask_s18a_bright_objectcenter
, masks.r_mask_s18a_bright_objectcenter
, masks.i_mask_s18a_bright_objectcenter
, masks.z_mask_s18a_bright_objectcenter
, masks.y_mask_s18a_bright_objectcenter

-- Blendedness
, meas2.i_blendedness_abs

-- Photo-z
, mizuki.photoz_best as mizuki_photoz_best
, mizuki.photoz_conf_best as mizuki_photoz_conf_best
, mizuki.photoz_std_best as mizuki_photoz_std_best
, mizuki.photoz_err95_min as mizuki_photoz_err95_min
, mizuki.photoz_err95_max as mizuki_photoz_err95_max

FROM
s19a_wide.meas as meas
LEFT JOIN s19a_wide.meas2 as meas2 using (object_id)
LEFT JOIN s19a_wide.meas3 as meas3 using (object_id)
LEFT JOIN s19a_wide.forced as forced using (object_id)
LEFT JOIN s19a_wide.masks as masks using (object_id)
LEFT JOIN s19a_wide.photoz_mizuki as mizuki using (object_id)

WHERE
meas.tract  IN  (9570,9571,9812,9813,9814,10054,10055,10056) AND
NOT meas.i_deblend_skipped                  AND
NOT meas2.i_sdsscentroid_flag               AND
NOT meas.i_pixelflags_edge                  AND
NOT meas.i_pixelflags_interpolatedcenter    AND
NOT meas.i_pixelflags_saturatedcenter       AND
NOT meas.i_pixelflags_crcenter              AND
NOT meas.i_pixelflags_bad                   AND
NOT meas.i_pixelflags_suspectcenter         AND
NOT meas.i_pixelflags_clipped               AND
meas.i_detect_isprimary                     AND
NOT meas2.i_hsmshaperegauss_flag            AND
meas.i_extendedness_value != 0              AND
-- In postgres, all numbers are comparable including NaN
meas2.i_hsmshaperegauss_sigma != 'NaN'
ORDER BY meas.object_id
