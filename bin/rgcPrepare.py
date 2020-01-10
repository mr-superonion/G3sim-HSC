#!/usr/bin/env python
import os
import galsim
import glob
import numpy as np
import astropy.io.fits as pyfits
from astropy.table import Table,vstack
import astropy.io.ascii as ascii

ngroup  =   8
Dir     =   '/work/xiangchong.li/work/S18ASimHSC/prepare2/'
psfDir  =   os.path.join(Dir,'psfPre')
pfname  =   os.path.join(Dir,'info.csv')

if not os.path.isdir(psfDir):
    os.system('mkdir -p %s' %psfDir)

def PreparePSF():
    print('Preparing PSF')
    nameList=   ['patch','tract','ra','dec','var','var_varplane','count_inputs']
    headList=   ['patch','tract','ra','dec','var','varPlane','count']

    nums    =   []
    fwhms   =   []
    var_meas=   []
    var_plane=  []
    tracts  =   []
    patches =   []
    groups  =   []

    fNlist  =   glob.glob('./from_Hironao/psf-catalogs/psf-catalog-*')
    nPerT   =   1
    for ifile,fname in enumerate(fNlist):
        preAll  =   pyfits.getdata(fname)
        if len(preAll)==0:
            os.popen('rm %s' %fname)
            continue
        tract   =   preAll[0]['tract']
        patch   =   preAll[0]['patch']
        imFname =   './from_Hironao/psf-images/psf-images-%s-%s.fits' %(tract,patch)
        if not os.path.exists(imFname):
            continue
        psfAll  =   pyfits.open(imFname)
        nsample =   len(preAll)
        for fid in range(nPerT):
            ifield  =   ifile*nPerT+fid
            outFname=   os.path.join(psfDir,'psf-%05d.fits' %ifield)
            np.random.seed(ifield)
            g1      =   np.random.rand()*0.06-0.03
            g2      =   np.random.rand()*0.06-0.03
            isp     =   np.random.randint(0,nsample)
            #For PSF
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
            hduOut  =   pyfits.PrimaryHDU(psfU)
            #Other information
            pre     =   preAll[isp]
            header  =   hduOut.header
            header['num']=  ifield
            header['group']=ifield%ngroup
            for name,hame in zip(nameList,headList):
                header[hame]=   pre[name]
            header['g1'] =   g1
            header['g2'] =   g2
            fwhm    =   2.355*psfImg.FindAdaptiveMom().moments_sigma
            header['fwhm']=  fwhm
            hduList =   pyfits.HDUList([hduOut])
            hduList.writeto(outFname,overwrite=True)
            nums.append(ifield)
            fwhms.append(header['fwhm'])
            var_meas.append(header['var'])
            var_plane.append(header['varPlane'])
            tracts.append(header['tract'])
            patches.append(header['patch'])
            groups.append(ifield%ngroup)
    names   =   ('num','fwhm','var','varPlane','tract','patch','group')
    data    =   (nums,fwhms,var_meas,var_plane,tracts,patches,groups)
    table   =   Table(data=data,names=names)
    table.write(pfname,overwrite=True)
    return

def PrepareNoiCor():
    print("preparing noise correaltions")
    #Do for tract
    info    =   Table.read('./info.csv')
    tractAll=   list(set(info['tract']))
    corAll  =   np.zeros((20,20))
    nAll    =   0
    for tract in tractAll:
        corTract    =   np.zeros((20,20))
        nInTract    =   0
        fnames  =   glob.glob('from_Hironao/psf-catalogs/psf-catalog-%s-*.fits' %(tract))
        for fn in fnames:
            preAll  =   pyfits.getdata(fn)
            corr    =   np.average(preAll['corr'],axis=0).reshape(21,21)
            corTract+=  corr[:-1,:-1]
            nInTract+=  1.
        corAll  +=  corTract
        nAll    +=  nInTract
        corTract/=  nInTract
        pyfits.writeto('./corPre/corTract-%s.fits' %(tract),corTract)
    pyfits.writeto('./corPre/corStack.fits',corAll/nAll)

    #Do for group
    corGroup=   np.zeros((ngroup,20,20))
    nPG     =   np.zeros(ngroup)
    for row in info:
        tract   =   row['tract']
        patch   =   row['patch']
        group   =   row['group']
        fn      =   'from_Hironao/psf-catalogs/psf-catalog-%s-%s.fits' %(tract,patch)
        preAll  =   pyfits.getdata(fn)
        corr    =   np.average(preAll['corr'],axis=0).reshape(21,21)
        corGroup[group]+=  corr[:-1,:-1]
        nPG[group]      +=  1.

    for ig in range(ngroup):
        pyfits.writeto('./corPre/corGroup-%02d.fits' %ig,corGroup[ig]/nPG[ig])
    return

def main():
    if not os.path.exists(pfname):
        PreparePSF()
    #PrepareNoiCor()
    return

if __name__=="__main__":
     main()
