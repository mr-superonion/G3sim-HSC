SELECT
        object_id
        ,f1.ra
        ,f1.dec
        ,f1.g_cmodel_flux
        ,f1.g_cmodel_fluxerr
        ,f1.r_cmodel_flux
        ,f1.r_cmodel_fluxerr
        ,f1.i_cmodel_flux
        ,f1.i_cmodel_fluxerr
        ,f1.z_cmodel_flux
        ,f1.z_cmodel_fluxerr
        ,f1.y_cmodel_flux
        ,f1.y_cmodel_fluxerr
        ,f2.g_psfflux_flux
        ,f2.g_psfflux_fluxerr
        ,f2.r_psfflux_flux
        ,f2.r_psfflux_fluxerr
        ,f2.i_psfflux_flux
        ,f2.i_psfflux_fluxerr
        ,f2.z_psfflux_flux
        ,f2.z_psfflux_fluxerr
        ,f2.y_psfflux_flux
        ,f2.y_psfflux_fluxerr
        ,f1.i_extendedness_value
		,f2.i_sdssshape_shape11
		,f2.i_sdssshape_shape22
		,f2.i_sdssshape_shape12
		,f2.i_sdssshape_psf_shape11
		,f2.i_sdssshape_psf_shape22
		,f2.i_sdssshape_psf_shape12
		,f2.r_sdssshape_shape11
		,f2.r_sdssshape_shape22
		,f2.r_sdssshape_shape12
		,f2.r_sdssshape_psf_shape11
		,f2.r_sdssshape_psf_shape22
		,f2.r_sdssshape_psf_shape12
		,f2.z_sdssshape_shape11
		,f2.z_sdssshape_shape22
		,f2.z_sdssshape_shape12
		,f2.z_sdssshape_psf_shape11
		,f2.z_sdssshape_psf_shape22
		,f2.z_sdssshape_psf_shape12
		,m1.r_calib_psf_candidate
		,m1.r_calib_psf_reserved
		,m1.r_calib_psf_used
		,m1.i_calib_psf_candidate
		,m1.i_calib_psf_reserved
		,m1.i_calib_psf_used
		,m1.z_calib_psf_candidate
		,m1.z_calib_psf_reserved
		,m1.z_calib_psf_used
		,m1.tract
		,m1.patch
        ,f1.a_g
        ,f1.a_r
        ,f1.a_i
        ,f1.a_z
        ,f1.a_y
    FROM
        s22a_test.forced  AS f1
      LEFT JOIN
        s22a_test.forced2 AS f2 USING (object_id)
      LEFT JOIN
        s22a_test.meas2   AS m2 USING (object_id)
	  LEFT JOIN
        s22a_test.meas    AS m1 USING (object_id)
	WHERE
        f1.isprimary
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
		AND f1.i_extendedness_value < 0.01
;
