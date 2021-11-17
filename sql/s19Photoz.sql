SELECT
  meas.object_id
, meas.parent_id
, meas.i_ra
, meas.i_dec
, meas.tract
, meas.patch

, mizuki.photoz_mean as mizuki_photoz_mean
, mizuki.photoz_mode as mizuki_photoz_mode
, mizuki.photoz_median as mizuki_photoz_median
, mizuki.photoz_best as mizuki_photoz_best
, mizuki.photoz_mc as mizuki_photoz_mc

, mizuki.photoz_conf_mean as mizuki_photoz_conf_mean
, mizuki.photoz_conf_mode as mizuki_photoz_conf_mode
, mizuki.photoz_conf_median as mizuki_photoz_conf_median
, mizuki.photoz_conf_best as mizuki_photoz_conf_best

, mizuki.photoz_risk_mean as mizuki_photoz_risk_mean
, mizuki.photoz_risk_mode as mizuki_photoz_risk_mode
, mizuki.photoz_risk_median as mizuki_photoz_risk_median
, mizuki.photoz_risk_best as mizuki_photoz_risk_best

, mizuki.photoz_std_mean as mizuki_photoz_std_mean
, mizuki.photoz_std_mode as mizuki_photoz_std_mode
, mizuki.photoz_std_median as mizuki_photoz_std_median
, mizuki.photoz_std_best as mizuki_photoz_std_best

, mizuki.photoz_err68_min as mizuki_photoz_err68_min
, mizuki.photoz_err68_max as mizuki_photoz_err68_max
, mizuki.photoz_err95_min as mizuki_photoz_err95_min
, mizuki.photoz_err95_max as mizuki_photoz_err95_max

, demp.photoz_mean as demp_photoz_mean
, demp.photoz_mode as demp_photoz_mode
, demp.photoz_median as demp_photoz_median
, demp.photoz_best as demp_photoz_best
, demp.photoz_mc as demp_photoz_mc

, demp.photoz_conf_mean as demp_photoz_conf_mean
, demp.photoz_conf_mode as demp_photoz_conf_mode
, demp.photoz_conf_median as demp_photoz_conf_median
, demp.photoz_conf_best as demp_photoz_conf_best

, demp.photoz_risk_mean as demp_photoz_risk_mean
, demp.photoz_risk_mode as demp_photoz_risk_mode
, demp.photoz_risk_median as demp_photoz_risk_median
, demp.photoz_risk_best as demp_photoz_risk_best

, demp.photoz_std_mean as demp_photoz_std_mean
, demp.photoz_std_mode as demp_photoz_std_mode
, demp.photoz_std_median as demp_photoz_std_median
, demp.photoz_std_best as demp_photoz_std_best

, demp.photoz_err68_min as demp_photoz_err68_min
, demp.photoz_err68_max as demp_photoz_err68_max
, demp.photoz_err95_min as demp_photoz_err95_min
, demp.photoz_err95_max as demp_photoz_err95_max

FROM
s19a_wide.meas as meas
LEFT JOIN s19a_wide.meas2 as meas2 using (object_id)
LEFT JOIN s19a_wide.photoz_v2_mizuki as mizuki using (object_id)
LEFT JOIN s19a_wide.photoz_v2_demp as demp using (object_id)

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
