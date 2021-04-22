SELECT
  meas.object_id
, meas.parent_id
, meas.i_ra
, meas.i_dec
, meas.i_variance_value
, meas.i_footprintarea_value
, meas.i_deblend_parenttoobig
, meas.i_extendedness_value
, meas.tract

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

,meas3.i_apertureflux_10_mag
,meas3.i_apertureflux_10_magerr
,meas3.i_apertureflux_15_mag
,meas3.i_apertureflux_15_magerr

,meas2.i_hsmpsfmoments_shape11
,meas2.i_hsmpsfmoments_shape22
,meas2.i_hsmpsfmoments_shape12

,meas2.i_hsmsourcemoments_shape11
,meas2.i_hsmsourcemoments_shape22
,meas2.i_hsmsourcemoments_shape12

-- forced measurement
, forced.a_g
, forced.a_r
, forced.a_i
, forced.a_z
, forced.a_y

,meas2.i_psfflux_mag
,meas2.i_psfflux_magerr

-- shapes
-- hsm regauss shapes
, meas2.i_hsmshaperegauss_e1
, meas2.i_hsmshaperegauss_e2
, meas2.i_hsmshaperegauss_sigma
, meas2.i_hsmshaperegauss_resolution

-- sdss shapes (without PSF correction)
, meas2.i_sdssshape_shape11
, meas2.i_sdssshape_shape22
, meas2.i_sdssshape_shape12
, meas2.i_sdssshape_psf_shape11
, meas2.i_sdssshape_psf_shape22
, meas2.i_sdssshape_psf_shape12

-- Input exposures number
, meas.g_inputcount_value
, meas.r_inputcount_value
, meas.i_inputcount_value
, meas.z_inputcount_value
, meas.y_inputcount_value

-- Mask plane
, meas.i_pixelflags_inexact_psfcenter
, meas.i_pixelflags_inexact_psf
, meas.i_pixelflags_clipped
, meas.i_pixelflags_edge
, meas.i_pixelflags_interpolatedcenter
, meas.i_pixelflags_saturatedcenter
, meas.i_pixelflags_crcenter
, meas.i_pixelflags_bad
, meas.i_pixelflags_suspectcenter

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

, mizukiv2.photoz_best as mizukiv2_photoz_best
, mizukiv2.photoz_conf_best as mizukiv2_photoz_conf_best
, mizukiv2.photoz_std_best as mizukiv2_photoz_std_best
, mizukiv2.photoz_err95_min as mizukiv2_photoz_err95_min
, mizukiv2.photoz_err95_max as mizukiv2_photoz_err95_max

, dempv2.photoz_best as dempv2_photoz_best
, dempv2.photoz_conf_best as dempv2_photoz_conf_best
, dempv2.photoz_std_best as dempv2_photoz_std_best
, dempv2.photoz_err95_min as dempv2_photoz_err95_min
, dempv2.photoz_err95_max as dempv2_photoz_err95_max

, dnnzv2.photoz_best as dnnzv2_photoz_best
, dnnzv2.photoz_conf_best as dnnzv2_photoz_conf_best
, dnnzv2.photoz_std_best as dnnzv2_photoz_std_best
, dnnzv2.photoz_err95_min as dnnzv2_photoz_err95_min
, dnnzv2.photoz_err95_max as dnnzv2_photoz_err95_max

FROM
s19a_dud_wide_depth_best.meas as meas
LEFT JOIN s19a_dud_wide_depth_best.meas2 as meas2 using (object_id)
LEFT JOIN s19a_dud_wide_depth_best.meas3 as meas3 using (object_id)
LEFT JOIN s19a_dud_wide_depth_best.forced as forced using (object_id)
LEFT JOIN s19a_wide.masks as masks using (object_id)
LEFT JOIN s19a_wide.photoz_mizuki as mizuki using (object_id)
LEFT JOIN s19a_wide.photoz_v2_mizuki as mizukiv2 using (object_id)
LEFT JOIN s19a_wide.photoz_v2_dnnz as dnnzv2 using (object_id)
LEFT JOIN s19a_wide.photoz_v2_demp as dempv2 using (object_id)

WHERE
meas.tract  IN  (9570,9571,9812,9813,9814,10054,10055,10056) AND
boxSearch(coord, 149, 151, 1.5, 3.0)        AND
NOT meas.i_deblend_skipped                  AND
NOT meas2.i_sdsscentroid_flag               AND

meas.i_detect_isprimary                     AND
NOT meas2.i_hsmshaperegauss_flag            AND
-- In postgres, all numbers are comparable including NaN
meas2.i_hsmshaperegauss_sigma != 'NaN'
ORDER BY meas.object_id
