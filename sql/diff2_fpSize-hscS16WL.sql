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
s16a_wide2.meas as meas
LEFT JOIN s16a_wide2.meas2 as meas2 using (object_id)
LEFT JOIN s16a_wide2.forced as forced using (object_id)

WHERE
meas.tract  IN  ({$tract})                  AND
NOT meas.ideblend_skipped                  AND
NOT meas.iflags_badcentroid                AND
NOT meas2.icentroid_sdss_flags             AND
NOT meas.iflags_pixel_edge                 AND
NOT meas.iflags_pixel_interpolated_center  AND
NOT meas.iflags_pixel_saturated_center     AND
NOT meas.iflags_pixel_cr_center            AND
NOT meas.iflags_pixel_bad                  AND
NOT meas.iflags_pixel_suspect_center       AND
NOT meas.iflags_pixel_clipped_any          AND
meas.idetect_is_primary                    AND
NOT meas2.ishape_hsm_regauss_flags         AND
meas.iclassification_extendedness != 0     AND
-- In postgres, all numbers are comparable including NaN
 meas2.ishape_hsm_regauss_sigma != 'NaN'    AND
forced.wl_fulldepth_fullcolor
ORDER BY meas.object_id
