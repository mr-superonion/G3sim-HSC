SELECT

  ran.object_id
, ran.ra
, ran.dec

-- columns which can be used for selection
, ran.tract
, ran.patch_s
, ran_masks.g_mask_s18a_bright_objectcenter
, ran_masks.r_mask_s18a_bright_objectcenter
, ran_masks.i_mask_s18a_bright_objectcenter
, ran_masks.z_mask_s18a_bright_objectcenter
, ran_masks.y_mask_s18a_bright_objectcenter
, ran.g_pixelflags_bright_objectcenter
, ran.r_pixelflags_bright_objectcenter
, ran.i_pixelflags_bright_objectcenter
, ran.z_pixelflags_bright_objectcenter
, ran.y_pixelflags_bright_objectcenter

FROM
s18a_wide.random as ran
LEFT JOIN s18a_wide.random_masks as ran_masks USING (object_id)

WHERE
meas.tract  IN  ({$tract})                              AND
NOT meas_ireprocess.i_sdsscentroid_flag                 AND
NOT meas_ireprocess.i_pixelflags_edge                   AND
NOT meas_ireprocess.i_pixelflags_interpolatedcenter     AND
NOT meas_ireprocess.i_pixelflags_saturatedcenter        AND
NOT meas_ireprocess.i_pixelflags_crcenter               AND
NOT meas_ireprocess.i_pixelflags_bad                    AND
NOT meas_ireprocess.i_pixelflags_suspectcenter          AND
meas_ireprocess.i_detect_isprimary                      AND
NOT meas_ireprocess.i_hsmshaperegauss_flag              AND
meas.i_extendedness_value != 0                          AND
 -- In postgres, all numbers are comparable including NaN
meas_ireprocess.i_hsmshaperegauss_sigma != 'NaN'
ORDER BY meas.object_id

NOT ran.i_pixelflags_edge               AND
NOT ran.i_pixelflags_interpolatedcenter AND
NOT ran.i_pixelflags_saturatedcenter    AND
NOT ran.i_pixelflags_crcenter           AND
NOT ran.i_pixelflags_bad                AND
NOT ran.i_pixelflags_suspectcenter      AND
ran.isprimary
ORDER BY ran.object_id
