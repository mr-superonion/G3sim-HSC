SELECT
  meas.object_id
, meas.parent_id
, meas.i_ra
, meas.i_dec
, meas.tract
, meas.patch
, meas.skymap_id

, meas.g_filterfraction_unweighted
, meas.r_filterfraction_unweighted
, meas.i_filterfraction_unweighted
, meas.z_filterfraction_unweighted
, meas.y_filterfraction_unweighted

, meas.g_filterfraction_weighted
, meas.r_filterfraction_weighted
, meas.i_filterfraction_weighted
, meas.z_filterfraction_weighted
, meas.y_filterfraction_weighted

, colorterm.i_offset
, colorterm.skymap_id
FROM
s19a_wide.meas as meas
LEFT JOIN s19a_wide.meas2 as meas2 using (object_id)
LEFT JOIN s19a_wide.colorterm_offset as colorterm using (skymap_id)

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
