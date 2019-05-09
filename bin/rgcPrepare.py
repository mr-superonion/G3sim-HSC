#!/usr/bin/env python
import numpy as np
import galsim
import astropy.io.fits as pyfits
import astropy.io.ascii as ascii

nfield  =   3200
psfAll  =   pyfits.open('rand_hsc_wide_psfs_rand.fits.gz')
rscalAll=   ascii.read('rescale_v4-rgc-nr.csv')['scale_factor']
nsample =   len(psfAll)
varAll  =   pyfits.getdata('variance_v2_rand.fits')
for ifield in range(nfield):
    #Noise variance, correlation function and psf
    np.random.seed(ifield)
    g1      =   np.random.rand()*0.06-0.03
    g2      =   np.random.rand()*0.06-0.03
    isp     =   np.random.randint(0,nsample)
    psf     =   psfAll[isp].data
    nx      =   psf.shape[1]
    ny      =   psf.shape[0]
    nxy     =   min(nx,ny)
    shiftx  =   (nx-nxy)//2
    shifty  =   (ny-nxy)//2
    psfU    =   np.zeros((nxy,nxy))
    psfU    =   psf[shifty:shifty+nxy,shiftx:shiftx+nxy]
    psfImg  =   galsim.ImageF(psfU)
    psfImg.scale=0.168
    fwhm    =   psfImg.calculateFWHM()#(arcsec)
    var     =   varAll[isp]
    hduOut  =   pyfits.PrimaryHDU(psfU)
    for name in varAll.names:
        nameU   =   name.split('_')[0]
        hduOut.header[nameU]=   var[name]
    hduOut.header['g1'] =   g1
    hduOut.header['g2'] =   g2
    hduOut.header['fwhm']=  fwhm
    hduOut.header['rescale']=rscalAll[isp]
    hduOut.writeto('psfPre/psf%04d.fits'%(ifield),overwrite=True)
    
