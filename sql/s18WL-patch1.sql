SELECT
  meas.object_id
, meas.parent_id
, meas.i_ra
, meas.i_dec

-- columns which can be used for selection
, meas.tract
, meas.patch
, meas.i_variance_value
, meas.i_footprintarea_value
, meas.i_cmodel_exp_niter
, meas.i_cmodel_dev_niter
, meas.i_cmodel_ellipse_11
, meas.i_cmodel_ellipse_22
, meas.i_cmodel_ellipse_12
, meas_ireprocess.i_deblend_parenttoobig

, forced3.g_apertureflux_10_mag as forced_g_apertureflux_10_mag
, forced3.g_apertureflux_10_magsigma as forced_g_apertureflux_10_magsigma
, forced3.r_apertureflux_10_mag as forced_r_apertureflux_10_mag
, forced3.r_apertureflux_10_magsigma as forced_r_apertureflux_10_magsigma
, forced3.i_apertureflux_10_mag as forced_i_apertureflux_10_mag
, forced3.i_apertureflux_10_magsigma as forced_i_apertureflux_10_magsigma
, forced3.z_apertureflux_10_mag as forced_z_apertureflux_10_mag
, forced3.z_apertureflux_10_magsigma as forced_z_apertureflux_10_magsigma
, forced3.y_apertureflux_10_mag as forced_y_apertureflux_10_mag
, forced3.y_apertureflux_10_magsigma as forced_y_apertureflux_10_magsigma

, forced3.g_apertureflux_15_mag as forced_g_apertureflux_15_mag
, forced3.g_apertureflux_15_magsigma as forced_g_apertureflux_15_magsigma
, forced3.r_apertureflux_15_mag as forced_r_apertureflux_15_mag
, forced3.r_apertureflux_15_magsigma as forced_r_apertureflux_15_magsigma
, forced3.i_apertureflux_15_mag as forced_i_apertureflux_15_mag
, forced3.i_apertureflux_15_magsigma as forced_i_apertureflux_15_magsigma
, forced3.z_apertureflux_15_mag as forced_z_apertureflux_15_mag
, forced3.z_apertureflux_15_magsigma as forced_z_apertureflux_15_magsigma
, forced3.y_apertureflux_15_mag as forced_y_apertureflux_15_mag
, forced3.y_apertureflux_15_magsigma as forced_y_apertureflux_15_magsigma

, forced3.g_apertureflux_20_mag as forced_g_apertureflux_20_mag
, forced3.g_apertureflux_20_magsigma as forced_g_apertureflux_20_magsigma
, forced3.r_apertureflux_20_mag as forced_r_apertureflux_20_mag
, forced3.r_apertureflux_20_magsigma as forced_r_apertureflux_20_magsigma
, forced3.i_apertureflux_20_mag as forced_i_apertureflux_20_mag
, forced3.i_apertureflux_20_magsigma as forced_i_apertureflux_20_magsigma
, forced3.z_apertureflux_20_mag as forced_z_apertureflux_20_mag
, forced3.z_apertureflux_20_magsigma as forced_z_apertureflux_20_magsigma
, forced3.y_apertureflux_20_mag as forced_y_apertureflux_20_mag
, forced3.y_apertureflux_20_magsigma as forced_y_apertureflux_20_magsigma

, forced3.g_apertureflux_30_mag as forced_g_apertureflux_30_mag
, forced3.g_apertureflux_30_magsigma as forced_g_apertureflux_30_magsigma
, forced3.r_apertureflux_30_mag as forced_r_apertureflux_30_mag
, forced3.r_apertureflux_30_magsigma as forced_r_apertureflux_30_magsigma
, forced3.i_apertureflux_30_mag as forced_i_apertureflux_30_mag
, forced3.i_apertureflux_30_magsigma as forced_i_apertureflux_30_magsigma
, forced3.z_apertureflux_30_mag as forced_z_apertureflux_30_mag
, forced3.z_apertureflux_30_magsigma as forced_z_apertureflux_30_magsigma
, forced3.y_apertureflux_30_mag as forced_y_apertureflux_30_mag
, forced3.y_apertureflux_30_magsigma as forced_y_apertureflux_30_magsigma

,meas3.i_apertureflux_10_mag
,meas3.i_apertureflux_10_magsigma
,meas3.i_apertureflux_15_mag
,meas3.i_apertureflux_15_magsigma
,meas3.i_apertureflux_20_mag
,meas3.i_apertureflux_20_magsigma
,meas3.i_apertureflux_30_mag
,meas3.i_apertureflux_30_magsigma

,meas2.i_psfflux_mag
,meas2.i_psfflux_magsigma

FROM
s18a_wide.meas as meas
LEFT JOIN s18a_wide.meas_ireprocess as meas_ireprocess using (object_id)
LEFT JOIN s18a_wide.forced3 as forced3 using (object_id)
LEFT JOIN s18a_wide.meas2 as meas2 using (object_id)
LEFT JOIN s18a_wide.meas3 as meas3 using (object_id)

WHERE
meas.tract  IN  ({$tract})                              AND
NOT meas_ireprocess.i_deblend_skipped                   AND
--NOT meas.i_cmodel_flag_badcentroid                   AND
--# This flag was not included in S16A, but if we find anything wrong, we can get it back.
NOT meas_ireprocess.i_sdsscentroid_flag                 AND
NOT meas_ireprocess.i_pixelflags_edge                   AND
NOT meas_ireprocess.i_pixelflags_interpolatedcenter     AND
NOT meas_ireprocess.i_pixelflags_saturatedcenter        AND
NOT meas_ireprocess.i_pixelflags_crcenter               AND
NOT meas_ireprocess.i_pixelflags_bad                    AND
NOT meas_ireprocess.i_pixelflags_suspectcenter          AND
NOT meas_ireprocess.i_pixelflags_clipped                AND
meas_ireprocess.i_detect_isprimary                      AND
NOT meas_ireprocess.i_hsmshaperegauss_flag              AND
meas.i_extendedness_value != 0                          AND
 -- In postgres, all numbers are comparable including NaN
meas_ireprocess.i_hsmshaperegauss_sigma != 'NaN'
ORDER BY meas.object_id
