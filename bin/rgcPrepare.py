#!/usr/bin/env python
import os
import galsim
import numpy as np
import astropy.io.fits as pyfits
from astropy.table import Table

nfield  =   800
psfAll  =   pyfits.open('rand_hsc_wide_psfs_rand.fits.gz')
rscalAll=   Table.read('rescale_v4-rgc-nr.csv')['scale_factor']
nsample =   len(psfAll)
assert len(psfAll)== nsample, 'number of rescaling parameters does not match number of PSF'
varAll  =   pyfits.getdata('variance_v2_rand.fits')
nums    =   range(nfield)
fields  =   [] 
fwhms   =   []
variances=  []
rescales=   []
for ifield in range(nfield):
    outFname=   'psfPre/psf%04d.fits'%(ifield)
    if not os.path.exists(outFname):
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
        header  =   hduOut.header
        for name in varAll.names:
            nameU   =   name.split('_')[0]
            header[nameU]=   var[name]
        header['g1'] =   g1
        header['g2'] =   g2
        header['fwhm']=  fwhm
        header['rescale']=rscalAll[isp]
        hduOut.writeto(outFname,overwrite=True)
    else:
        header  =   pyfits.getheader(outFname)
    fields.append(header['field'])
    fwhms.append(header['fwhm'])
    variances.append(header['var'])
    rescales.append(header['rescale'])
names   =   ('num','field','fwhm','var','rescale')
data    =   (nums,fields,fwhms,variances,rescales)
Table(data=data,names=names).write('psfPre/info.csv',overwrite=True)
    
