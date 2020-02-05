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
, meas.i_deblend_parenttoobig

, forced3.g_apertureflux_10_mag as forced_g_apertureflux_10_mag
, forced3.g_apertureflux_10_magerr as forced_g_apertureflux_10_magerr
, forced3.r_apertureflux_10_mag as forced_r_apertureflux_10_mag
, forced3.r_apertureflux_10_magerr as forced_r_apertureflux_10_magerr
, forced3.i_apertureflux_10_mag as forced_i_apertureflux_10_mag
, forced3.i_apertureflux_10_magerr as forced_i_apertureflux_10_magerr
, forced3.z_apertureflux_10_mag as forced_z_apertureflux_10_mag
, forced3.z_apertureflux_10_magerr as forced_z_apertureflux_10_magerr
, forced3.y_apertureflux_10_mag as forced_y_apertureflux_10_mag
, forced3.y_apertureflux_10_magerr as forced_y_apertureflux_10_magerr

, forced3.g_apertureflux_15_mag as forced_g_apertureflux_15_mag
, forced3.g_apertureflux_15_magerr as forced_g_apertureflux_15_magerr
, forced3.r_apertureflux_15_mag as forced_r_apertureflux_15_mag
, forced3.r_apertureflux_15_magerr as forced_r_apertureflux_15_magerr
, forced3.i_apertureflux_15_mag as forced_i_apertureflux_15_mag
, forced3.i_apertureflux_15_magerr as forced_i_apertureflux_15_magerr
, forced3.z_apertureflux_15_mag as forced_z_apertureflux_15_mag
, forced3.z_apertureflux_15_magerr as forced_z_apertureflux_15_magerr
, forced3.y_apertureflux_15_mag as forced_y_apertureflux_15_mag
, forced3.y_apertureflux_15_magerr as forced_y_apertureflux_15_magerr

, forced3.g_apertureflux_20_mag as forced_g_apertureflux_20_mag
, forced3.g_apertureflux_20_magerr as forced_g_apertureflux_20_magerr
, forced3.r_apertureflux_20_mag as forced_r_apertureflux_20_mag
, forced3.r_apertureflux_20_magerr as forced_r_apertureflux_20_magerr
, forced3.i_apertureflux_20_mag as forced_i_apertureflux_20_mag
, forced3.i_apertureflux_20_magerr as forced_i_apertureflux_20_magerr
, forced3.z_apertureflux_20_mag as forced_z_apertureflux_20_mag
, forced3.z_apertureflux_20_magerr as forced_z_apertureflux_20_magerr
, forced3.y_apertureflux_20_mag as forced_y_apertureflux_20_mag
, forced3.y_apertureflux_20_magerr as forced_y_apertureflux_20_magerr

, forced3.g_apertureflux_30_mag as forced_g_apertureflux_30_mag
, forced3.g_apertureflux_30_magerr as forced_g_apertureflux_30_magerr
, forced3.r_apertureflux_30_mag as forced_r_apertureflux_30_mag
, forced3.r_apertureflux_30_magerr as forced_r_apertureflux_30_magerr
, forced3.i_apertureflux_30_mag as forced_i_apertureflux_30_mag
, forced3.i_apertureflux_30_magerr as forced_i_apertureflux_30_magerr
, forced3.z_apertureflux_30_mag as forced_z_apertureflux_30_mag
, forced3.z_apertureflux_30_magerr as forced_z_apertureflux_30_magerr
, forced3.y_apertureflux_30_mag as forced_y_apertureflux_30_mag
, forced3.y_apertureflux_30_magerr as forced_y_apertureflux_30_magerr

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

FROM
s19a_wide.meas as meas
LEFT JOIN s19a_wide.forced3 as forced3 using (object_id)
LEFT JOIN s19a_wide.meas3 as meas3 using (object_id)
LEFT JOIN s19a_wide.meas2 as meas2 using (object_id)

WHERE
meas.tract  IN  ({$tract})                  AND
NOT meas.i_deblend_skipped                  AND
-- NOT meas.i_cmodel_flag_badcentroid       AND # This flag was not included in S16A, but if we find anything wrong, we can get it back.
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
