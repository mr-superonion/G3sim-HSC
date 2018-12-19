import astropy.io.fits as pyfits



def main():
    nfield  =   800
    psfAll  =   pyfits.open('rand_hsc_wide_psfs_rand.fits.gz')
    nsample =   len(psfAll)
    varAll  =   pyfits.getdata('variance_v2_rand.fits')
    for ifield in range(nfield):
        #Noise variance, correlation function and psf
        np.random.seed(ifield)
        isp     =   np.random.randint(0,nsample)
        psf     =   psfAll[isp].data
        var     =   varAll[isp]
        hduOut  =   fits.PrimaryHDU(psf)
        for name in varAll.names:
            nameU   =   name.split('_')[0]
            hduOut.header[numeU]=   varAll[nameU]
        hduOut.writeto('psfPre/psf%04d.fits'%(ifield))
    return
    
if __name__ == '__main__':
    main()
