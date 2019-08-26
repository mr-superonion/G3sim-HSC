#!/usr/bin/env python
import os
import galsim
import numpy as np
import astropy.io.fits as pyfits
from astropy.table import Table

nameList=   ['field_name','patch','tract','ra','dec','var','var_varplane','count_inputs']
headList=   ['field','patch','tract','ra','dec','var','varPlane','count']

nums    =   [] 
fields  =   [] 
fwhms   =   []
var_meas=   []
var_plane=  []

tDir    =   './pre/tract/'
fNlist  =   os.popen('ls ./pre/tract/').readlines()
nPerT   =   10
for ifile,fname in enumerate(fNlist):
    fname   =   fname[:-1]
    fname   =   os.path.join(tDir,fname)
    preAll  =   pyfits.getdata(fname)
    corr    =   np.average(preAll['corr'],axis=0).reshape(21,21)
    corHDU  =   pyfits.ImageHDU(corr) 
    nsample =   len(preAll)
    for fid in range(nPerT):
        ifield  =   ifile*nPerT+fid
        outFname=   './pre/out/out%04d.fits'%(ifield)
        np.random.seed(ifield)
        g1      =   np.random.rand()*0.06-0.03
        g2      =   np.random.rand()*0.06-0.03
        isp     =   np.random.randint(0,nsample)
        pre     =   preAll[isp]
        fwhm    =   2.355*pre['psf_size']*0.168
        psfU    =   None
        hduOut  =   pyfits.PrimaryHDU(psfU)
        header  =   hduOut.header
        for name,hame in zip(nameList,headList):
            header[hame]=   pre[name]
        header['g1'] =   g1
        header['g2'] =   g2
        header['fwhm']=  fwhm
        hduList =   pyfits.HDUList([hduOut,corHDU])
        hduList.writeto(outFname,overwrite=True)
        nums.append(ifield)
        fields.append(header['field'])
        fwhms.append(header['fwhm'])
        var_meas.append(header['var'])
        var_plane.append(header['varPlane'])
names   =   ('num','field','fwhm','var','varPlane')
data    =   (nums,fields,fwhms,var_meas,var_plane)
Table(data=data,names=names).write('./pre/info.csv',overwrite=True)


"""
for ifield in range(nfield):
    if not os.path.exists(outFname):
        #Noise variance, correlation function and psf
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
        hduOut.writeto(outFname,overwrite=True)
    else:
        header  =   pyfits.getheader(outFname)
    
"""
