"""
This code checks whether a galaxy in a particular postage stamp passes the loose S/N cuts that we
want to use for HSC (S/N<100).

It also stores some information about the PSF for the truth tables.
"""
import galsim
import numpy as np

def HSCReject(config, base, value_type):
    """See if a drawn postage stamp should be rejected based on the desired HSC cuts on S/N.
    """
    # Note: this is the pre-noise rendering of the object
    im = base['current_stamp']
    psf = galsim.config.GetCurrentValue('psf',base)

    reject = False

    # Check the resolution.  This depends on the PSF-convolved galaxy image and the PSF image.  The
    # latter will be the same for all galaxies in the same image, so we can draw the PSF image and
    # cache it, or retrieve it from the cache if it's already there.
    if 'cached_psf_image' in base and psf is base['cached_psf']:
        psf_image = base['cached_psf_image']
    else:
        psf_image = psf.drawImage(scale=im.scale, method='no_pixel')
        base['cached_psf_image'] = psf_image
        base['cached_psf'] = psf
        base['cached_mom'] = psf_image.FindAdaptiveMom()
    # While we are at it, save some more information.
    base['psf_sigma'] = base['cached_mom'].moments_sigma
    base['psf_e1'] = base['cached_mom'].observed_shape.e1
    base['psf_e2'] = base['cached_mom'].observed_shape.e2
    # We will define resolution using the resolution factor from re-Gaussianization, which uses a
    # trace-based size estimate of the post-seeing object.  For simplicity, we run
    # re-Gaussianization directly.
    res = galsim.hsm.EstimateShear(im, psf_image, strict='False')
    # Cut based on the resolution failure, or if the measurement simply failed.
    #if res.error_message != "" or res.resolution_factor < 0.1:
    #    reject = True

    if not reject:
        # Next (if the object is still considered viable) we check the SNR using the idealized SNR
        # estimator used for GREAT3.  For a derivation of this calculation, see the docstring for
        # the addNoiseSNR method.  This result is too optimistic compared to one that might be used
        # in reality.  Also, for realistic galaxies, it will be confused by the noise in the image.
        # In GREAT3, we dealt with this by precomputing SNR values using parametric galaxies models
        # and using them even when calculating SNR for realistic galaxies.  It is hard to imagine
        # how to implement this here.  We could switch to using some kind of Gaussian-weighted SNR
        # estimator for the realistic galaxy branches perhaps.  But for now we just do the same
        # calculation, but subtract off the contribution of the original image noise from the
        # calculation of sum(I(x,y)^2).

        # Note: aside from the bit about the noise for real-galaxy branches, this calculation
        # could be replaced by setting
        #     min_snr: 5
        #     max_snr: 100
        # in the stamp field of the config file.

        # First get the sum of pixel values squared:
        sumsq = np.sum(im.array**2)
        # Get the noise variance for this stamp:
        var = galsim.config.CalculateNoiseVar(base)

        # If the galaxy has a noise attribute, then we are in a real-galaxy branch, so we need to
        # subtract off the contribution of the image noise to sumsq.
        gal = galsim.config.GetCurrentValue('gal',base)
        if hasattr(gal, 'noise'):
            final = galsim.Convolve(gal,psf)
            current_noise_var = final.noise.getVariance()
            npix = np.prod(im.array.shape)
            sumsq -= current_noise_var * npix

        snr = np.sqrt(sumsq / var)
        base['ideal_snr'] = snr
        base['resolution'] = res.resolution_factor

        if snr > 100.:
        #if snr < 5. or snr > 100.:
            reject = True

    # False here is a "safe" parameter, which tells GalSim not to cache this result to use again
    # for the next object.
    print(reject)
    return reject, False


galsim.config.RegisterValueType('HSCReject', HSCReject, [ bool ])

