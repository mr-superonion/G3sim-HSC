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

FROM
s18a_wide.meas as meas
LEFT JOIN s18a_wide.meas_ireprocess as meas_ireprocess using (object_id)

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
