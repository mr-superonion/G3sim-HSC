#!/usr/bin/env python
import os
import itertools
import numpy as np
import astropy.io.fits as pyfits
import astropy.io.ascii as ascii
from multiprocessing import Pool
from astropy.table import Table
import lsst.afw.image as afwImg

def histoGramWap(arg):
    return histoGram(*arg)

def histoGram(ifile):
    dataOut =   np.zeros((4,59))
    rescat =   '/work/xiangchong.li/work/S16ASimHSC/sample4/src-%04d.fits'%(ifile)
    if not os.path.exists(rescat):
        return dataOut 
    data    =   pyfits.getdata(rescat)
    data    =   data[data['cmodel_flux']>0.]
    snr     =   data['cmodel_flux']/data['cmodel_flux_err']
    Res     =   data['shape_hsm_regauss_resolution'] 
    #prepare the frc
    dataCalib = afwImg.Calib()
    dataCalib.setFluxMag0(63095734448.0194)
    mag     =   np.ones(len(data))*30
    for i in range(len(mag)):
        mag[i]=dataCalib.getMagnitude(data['cmodel_flux'][i]) 
    mask    =   np.ones(len(data)).astype(bool)
    mask    =   mask&(data['deblend_nChild']==0)
    #mask    =   mask&(data['fpfs_moments'][:,0]>600)
    mask    =   mask&(mag<24.5)
    #mask    =   mask&(data['classification_extendedness']>0)
    mask    =   mask&(~np.isnan(data['shape_hsm_regauss_resolution']))
    mask    =   mask&(~np.isnan(data['cmodel_flux']))
    mask    =   mask&(~np.isnan(data['cmodel_flux_err']))
    mask    =   mask&(~np.isnan(mag))
    mask    =   mask&(~np.isnan(data['shape_hsm_regauss_e1']))
    mask    =   mask&(~np.isnan(data['shape_hsm_regauss_e2']))
    mask    =   mask&(Res>0.3)
    mask    =   mask&(snr>10.)
    dx      =   abs(data['centroid_sdss'][:,0]%64-32)
    dy      =   abs(data['centroid_sdss'][:,1]%64-32)
    mask    =   mask&(dx**2.+dy**2.<400)
    data    =   data[mask]
    mag     =   mag[mask]
    Res     =   Res[mask]
    snr     =   snr[mask]
    """
    moments =   data['fpfs_moments']
    cRatio  =   4
    const   =   (np.load('./CBase.npy'))*cRatio
    # get weight
    weight  =   moments[:,0]+const
    # get magnitude
    e1      =   moments[:,1]/weight
    e2      =   moments[:,2]/weight
    # get response factor 
    R1      =   1./np.sqrt(2.)*(moments[:,0]-moments[:,3])/weight+np.sqrt(2)*(e1**2.)
    R2      =   1./np.sqrt(2.)*(moments[:,0]-moments[:,3])/weight+np.sqrt(2)*(e2**2.)
    mask2   =   abs((R1+R2)/2.)<5.
    data    =   data[mask2]
    """
    absE    =   data['shape_hsm_regauss_e1']**2.+data['shape_hsm_regauss_e2']**2.
    absE    =   np.sqrt(absE)
    ells    =   [Res,snr,mag,absE]
    mGrids  =   np.load('./histGridsHSM.npy')
    for i in range(4):
        dataOut[i,:]=np.histogram(ells[i],bins=mGrids[i])[0]
    return dataOut

if __name__ == '__main__':
    ndat    =   800
    files   =   range(ndat)
    p       =   Pool(8)
    ngrid   =   60
    gFname  =   './histGridsHSM.npy' 
    hgFname =   './histGramHSMS16a.npy' 
    mGrids  =   []
    mGrids.append(np.linspace(0.1,1.,ngrid))
    mGrids.append(np.linspace(0.,80.,ngrid))
    mGrids.append(np.linspace(22.,25.2,ngrid))
    mGrids.append(np.linspace(0.,2.,ngrid))
    mGrids  =   np.vstack(mGrids)
    np.save(gFname,mGrids)

    outDat  =   p.map(histoGram,files)
    hg      =   np.zeros((4,ngrid-1))
    outDat  =   [x for x in outDat if x is not None]
    weiFname=   'reweight_s16a_v4.txt'
    if os.path.exists(weiFname):
        weights =   ascii.read()['weight']
    else:
        weights  =   np.ones(len(outDat))

    for i,dat in enumerate(outDat):
        hg+=dat*weights[i]
    for i in range(4):
        hg[i,:]=hg[i,:]/sum(hg[i,:])
    np.save(hgFname,hg)
