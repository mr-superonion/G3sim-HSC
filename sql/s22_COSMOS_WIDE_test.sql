SELECT
        object_id
        ,f1.ra
        ,f1.dec
        ,f1.r_cmodel_flux
        ,f1.r_cmodel_fluxerr
        ,f1.i_cmodel_flux
        ,f1.i_cmodel_fluxerr
        ,f1.z_cmodel_flux
        ,f1.z_cmodel_fluxerr
        ,f2.r_psfflux_flux
        ,f2.r_psfflux_fluxerr
        ,f2.i_psfflux_flux
        ,f2.i_psfflux_fluxerr
        ,f2.z_psfflux_flux
        ,f2.z_psfflux_fluxerr
		,f3.r_apertureflux_10_flux
		,f3.r_apertureflux_10_fluxerr
		,f3.r_apertureflux_15_flux
		,f3.r_apertureflux_15_fluxerr
		,f3.i_apertureflux_10_flux
		,f3.i_apertureflux_10_fluxerr
		,f3.i_apertureflux_15_flux
		,f3.i_apertureflux_15_fluxerr
		,f3.z_apertureflux_10_flux
		,f3.z_apertureflux_10_fluxerr
		,f3.z_apertureflux_15_flux
		,f3.z_apertureflux_15_fluxerr
        ,f1.i_extendedness_value
        ,f1.a_r
        ,f1.a_i
        ,f1.a_z
        -- hsm regauss shapes
        , m2.r_hsmshaperegauss_e1
        , m2.r_hsmshaperegauss_e2
        , m2.r_hsmshaperegauss_sigma
        , m2.r_hsmshaperegauss_resolution
        , m2.i_hsmshaperegauss_e1
        , m2.i_hsmshaperegauss_e2
        , m2.i_hsmshaperegauss_sigma
        , m2.i_hsmshaperegauss_resolution
        , m2.z_hsmshaperegauss_e1
        , m2.z_hsmshaperegauss_e2
        , m2.z_hsmshaperegauss_sigma
        , m2.z_hsmshaperegauss_resolution
    FROM
        s22a_test.forced  AS f1
      LEFT JOIN
        s22a_test.forced2 AS f2 USING (object_id)
	  LEFT JOIN
		s22a_test.forced3 AS f3 USING (object_id)
	  LEFT JOIN
        s22a_test.meas    AS m1 USING (object_id)
	  LEFT JOIN
        s22a_test.meas2   AS m2 USING (object_id)
    WHERE
        f1.isprimary
		AND m1.tract  IN  (9570,9571,9812,9813,9814,10054,10055,10056)
        AND NOT m2.g_sdsscentroid_flag
        AND NOT m2.r_sdsscentroid_flag
        AND NOT m2.i_sdsscentroid_flag
        AND NOT m2.z_sdsscentroid_flag
        AND NOT m2.y_sdsscentroid_flag
        AND NOT f1.g_pixelflags_edge
        AND NOT f1.r_pixelflags_edge
        AND NOT f1.i_pixelflags_edge
        AND NOT f1.z_pixelflags_edge
        AND NOT f1.y_pixelflags_edge
        AND NOT f1.g_pixelflags_interpolatedcenter /**/
        AND NOT f1.r_pixelflags_interpolatedcenter /**/
        AND NOT f1.i_pixelflags_interpolatedcenter /**/
        AND NOT f1.z_pixelflags_interpolatedcenter /**/
        AND NOT f1.y_pixelflags_interpolatedcenter /**/
        AND NOT f1.g_pixelflags_saturatedcenter
        AND NOT f1.r_pixelflags_saturatedcenter
        AND NOT f1.i_pixelflags_saturatedcenter
        AND NOT f1.z_pixelflags_saturatedcenter
        AND NOT f1.y_pixelflags_saturatedcenter
        AND NOT f1.g_pixelflags_crcenter /**/
        AND NOT f1.r_pixelflags_crcenter /**/
        AND NOT f1.i_pixelflags_crcenter /**/
        AND NOT f1.z_pixelflags_crcenter /**/
        AND NOT f1.y_pixelflags_crcenter /**/
        AND NOT f1.g_pixelflags_bad
        AND NOT f1.r_pixelflags_bad
        AND NOT f1.i_pixelflags_bad
        AND NOT f1.z_pixelflags_bad
        AND NOT f1.y_pixelflags_bad
        AND NOT f1.g_cmodel_flag
        AND NOT f1.r_cmodel_flag
        AND NOT f1.i_cmodel_flag
        AND NOT f1.z_cmodel_flag
        AND NOT f1.y_cmodel_flag
		AND f1.i_extendedness_value > 0.01
;
