SELECT
  meas.object_id
, meas.parent_id
, meas.i_ra
, meas.i_dec
, meas.i_variance_value
, meas.i_footprintarea_value
, meas.i_deblend_parenttoobig

-- photometry
-- unforced kron magnitudes
, meas2.i_kronflux_flux
, meas2.i_kronflux_fluxsigma
, meas2.i_kronflux_flag
, meas2.i_kronflux_mag
, meas2.i_kronflux_magsigma

-- unforced CModel magnitudes
, meas.i_cmodel_mag
, meas.i_cmodel_magsigma
, meas.i_cmodel_flag
, meas.i_cmodel_flag_badcentroid
, meas.i_cmodel_flux
, meas.i_cmodel_fluxsigma

-- unforced Aperture fluxes
,meas3.i_apertureflux_10_flux
,meas3.i_apertureflux_10_fluxsigma
,meas3.i_apertureflux_15_flux
,meas3.i_apertureflux_15_fluxsigma
,meas3.i_apertureflux_20_flux
,meas3.i_apertureflux_20_fluxsigma
,meas3.i_apertureflux_30_flux
,meas3.i_apertureflux_30_fluxsigma
,meas3.i_apertureflux_40_flux
,meas3.i_apertureflux_40_fluxsigma

-- forced measurement
, forced.a_g
, forced.a_r
, forced.a_i
, forced.a_z
, forced.a_y

-- Local background Flux
, meas.i_localbackground_flux
, meas.i_localbackground_fluxsigma
, meas.i_localbackground_mag
, meas.i_localbackground_flag

-- forced Kron magnitudes
, forced2.g_kronflux_mag    as forced_g_kronflux_mag
, forced2.g_kronflux_magsigma   as forced_g_kronflux_magsigma
, forced2.g_kronflux_flag       as forced_g_kronflux_flag

, forced2.r_kronflux_mag    as forced_r_kronflux_mag
, forced2.r_kronflux_magsigma   as forced_r_kronflux_magsigma
, forced2.r_kronflux_flag       as forced_r_kronflux_flag

, forced2.i_kronflux_mag    as forced_i_kronflux_mag
, forced2.i_kronflux_magsigma   as forced_i_kronflux_magsigma
, forced2.i_kronflux_flag       as forced_i_kronflux_flag

, forced2.z_kronflux_mag    as forced_z_kronflux_mag
, forced2.z_kronflux_magsigma   as forced_z_kronflux_magsigma
, forced2.z_kronflux_flag       as forced_z_kronflux_flag

, forced2.y_kronflux_mag    as forced_y_kronflux_mag
, forced2.y_kronflux_magsigma   as forced_y_kronflux_magsigma
, forced2.y_kronflux_flag       as forced_y_kronflux_flag

-- forced CModel magnitudes and fluxes
, forced.g_cmodel_mag           as forced_g_cmodel_mag
, forced.g_cmodel_magsigma      as forced_g_cmodel_magsigma
, forced.g_cmodel_flux      as forced_g_cmodel_flux
, forced.g_cmodel_fluxsigma     as forced_g_cmodel_fluxsigma
, forced.g_cmodel_flag      as forced_g_cmodel_flag

, forced.r_cmodel_mag           as forced_r_cmodel_mag
, forced.r_cmodel_magsigma      as forced_r_cmodel_magsigma
, forced.r_cmodel_flux      as forced_r_cmodel_flux
, forced.r_cmodel_fluxsigma     as forced_r_cmodel_fluxsigma
, forced.r_cmodel_flag      as forced_r_cmodel_flag

, forced.i_cmodel_mag           as forced_i_cmodel_mag
, forced.i_cmodel_magsigma      as forced_i_cmodel_magsigma
, forced.i_cmodel_flux      as forced_i_cmodel_flux
, forced.i_cmodel_fluxsigma     as forced_i_cmodel_fluxsigma
, forced.i_cmodel_flag      as forced_i_cmodel_flag

, forced.z_cmodel_mag           as forced_z_cmodel_mag
, forced.z_cmodel_magsigma      as forced_z_cmodel_magsigma
, forced.z_cmodel_flux      as forced_z_cmodel_flux
, forced.z_cmodel_fluxsigma     as forced_z_cmodel_fluxsigma
, forced.z_cmodel_flag      as forced_z_cmodel_flag

, forced.y_cmodel_mag           as forced_y_cmodel_mag
, forced.y_cmodel_magsigma      as forced_y_cmodel_magsigma
, forced.y_cmodel_flux      as forced_y_cmodel_flux
, forced.y_cmodel_fluxsigma     as forced_y_cmodel_fluxsigma
, forced.y_cmodel_flag      as forced_y_cmodel_flag

-- shapes
-- hsm regauss shapes
, meas2.i_hsmshaperegauss_e1
, meas2.i_hsmshaperegauss_e2
, meas2.i_hsmshaperegauss_sigma
, meas2.i_hsmshaperegauss_resolution
-- derived hsm regauss quantitites
-- weight assuming rms_e = 0.365
, 1./(meas2.i_hsmshaperegauss_sigma*meas2.i_hsmshaperegauss_sigma + 0.365*0.365) as i_hsmshaperegauss_derived_weight
, meas2.i_hsmshaperegauss_sigma         as i_hsmshaperegauss_derived_sigma_e
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

, meas2.i_hsmsourcemoments_shape11
, meas2.i_hsmsourcemoments_shape22
, meas2.i_hsmsourcemoments_shape12
, meas2.i_hsmpsfmoments_shape11
, meas2.i_hsmpsfmoments_shape22
, meas2.i_hsmpsfmoments_shape12

-- flags
-- weak lensing cut
, (meas.i_cmodel_flux/meas.i_cmodel_fluxsigma > 10) AND
(meas2.i_hsmshaperegauss_resolution > 0.3) AND
(meas2.i_hsmshaperegauss_e1*meas2.i_hsmshaperegauss_e1 +  meas2.i_hsmshaperegauss_e2*meas2.i_hsmshaperegauss_e2 < 4) AND
(meas2.i_hsmshaperegauss_sigma > 0.) AND
(meas2.i_hsmshaperegauss_sigma < 0.4) AND
(meas.i_cmodel_mag - forced.a_i < 24.5) as weak_lensing_flag

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
, meas.i_pixelflags_bright_objectcenter
, meas.i_pixelflags_bright_object
, meas.i_pixelflags_inexact_psfcenter
, meas.i_pixelflags_inexact_psf
-- Bright object mask (Arcturus)
, masks.g_mask_s18a_bright_objectcenter
, masks.r_mask_s18a_bright_objectcenter
, masks.i_mask_s18a_bright_objectcenter
, masks.z_mask_s18a_bright_objectcenter
, masks.y_mask_s18a_bright_objectcenter

-- Blendedness
, meas.i_blendedness_old
, meas.i_blendedness_raw_flux
, meas.i_blendedness_abs_flux
, meas.i_blendedness_flag

, meas4.g_convolvedflux_2_15_mag
, meas4.r_convolvedflux_2_15_mag
, meas4.i_convolvedflux_2_15_mag
, meas4.z_convolvedflux_2_15_mag
, meas4.y_convolvedflux_2_15_mag

FROM
s18a_wide.meas as meas
LEFT JOIN s18a_wide.meas2 as meas2 using (object_id)
LEFT JOIN s18a_wide.meas3 as meas3 using (object_id)
LEFT JOIN s18a_wide.meas4 as meas4 using (object_id)
LEFT JOIN s18a_wide.forced as forced using (object_id)
LEFT JOIN s18a_wide.forced2 as forced2 using (object_id)
LEFT JOIN s18a_wide.masks as masks using (object_id)

WHERE
meas.tract  IN  ({$tract})                  AND
NOT meas.i_deblend_skipped                  AND
 -- NOT meas.i_cmodel_flag_badcentroid      AND # This flag was not included in S16A, but if we find anything wrong, we can get it back.
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
