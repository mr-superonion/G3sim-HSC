#!/usr/bin/env python
import os
import itertools
import numpy as np
import astropy.io.fits as pyfits
import astropy.io.ascii as ascii
from multiprocessing import Pool
from astropy.table import Table

def histoGramWap(arg):
    return histoGram(*arg)

def histoGram(ifile,sName):
    dataOut =   np.zeros((4,59))
    rescat =   './%s/src-%04d.fits'%(sName,ifile)
    if not os.path.exists(rescat):
        return dataOut 
    data    =   pyfits.getdata(rescat)
    mask    =   np.ones(len(data)).astype(bool)
    mask    =   mask&(np.all(~np.isnan(data['fps_momentsG']),axis=1))
    mask    =   mask&(data['deblend_nChild']==0)
    mask    =   mask&(data['fps_momentsG'][:,0]>600)
    mask    =   mask&(data['CModel_mag']<25.)
    mask    =   mask&(data['base_ClassificationExtendedness_value']>0)
    mask    =   mask&(~np.isnan(data['ext_shapeHSM_HsmSourceMoments_xx']))
    mask    =   mask&(~np.isnan(data['ext_shapeHSM_HsmSourceMoments_yy']))
    mask    =   mask&(~np.isnan(data['ext_shapeHSM_HsmPsfMoments_xx']))
    mask    =   mask&(~np.isnan(data['ext_shapeHSM_HsmPsfMoments_yy']))
    mask    =   mask&(~np.isnan(data['modelfit_CModel_flux']))
    mask    =   mask&(~np.isnan(data['modelfit_CModel_fluxSigma']))
    mask    =   mask&(~np.isnan(data['CModel_mag']))
    mask    =   mask&(~np.isnan(data['ext_shapeHSM_HsmShapeRegauss_e1']))
    mask    =   mask&(~np.isnan(data['ext_shapeHSM_HsmShapeRegauss_e2']))
    mask    =   mask&(data['base_Blendedness_abs_flux']<10**(-0.375))
    mask    =   mask&(data['matchMinD']<3.)
    #mask    =   mask&(data['matchResolution']>0.3)
    if '4' in sName:
        dx    =   abs(data['base_SdssCentroid_x']%64-32)
        dy    =   abs(data['base_SdssCentroid_y']%64-32)
        mask  = mask&(dx**2.+dy**2.<400)
    data    =   data[mask]
    moments =   data['fps_momentsG']
    cRatio  =   4
    const   =   (np.load('./s4Calib/CBase.npy'))*cRatio
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
    absE    =   data['ext_shapeHSM_HsmShapeRegauss_e1']**2.+data['ext_shapeHSM_HsmShapeRegauss_e2']**2.
    absE    =   np.sqrt(absE)
    """
    Res     =   (data['ext_shapeHSM_HsmPsfMoments_xx']+data['ext_shapeHSM_HsmPsfMoments_yy'])
    Res     =   Res/(data['ext_shapeHSM_HsmSourceMoments_xx']+data['ext_shapeHSM_HsmSourceMoments_yy'])
    Res     =   1.-Res
    """
    Res     =   data['matchResolution']
    mag     =   data['CModel_mag']
    snr     =   data['modelfit_CModel_flux']/data['modelfit_CModel_fluxSigma']
    ells    =   [Res,snr,mag,absE]
    mGrids  =   np.load('./%s/histGridsHSM.npy'%(sName))
    for i in range(4):
        dataOut[i,:]=np.histogram(ells[i],bins=mGrids[i])[0]
    return dataOut

if __name__ == '__main__':
    ndat    =   800
    files   =   range(ndat)
    p       =   Pool(24)
    ngrid   =   60
    sName   =   'sample4-matched'
    sNameA  =   itertools.repeat(sName)
    gFname  =   './%s/histGridsHSM.npy' %(sName)
    hgFname =   './%s/histGramHSMM600Mag25.npy'%(sName) 
    mGrids  =   []
    mGrids.append(np.linspace(0.1,1.,ngrid))
    mGrids.append(np.linspace(0.,80.,ngrid))
    mGrids.append(np.linspace(22.,25.2,ngrid))
    mGrids.append(np.linspace(0.,2.,ngrid))
    mGrids  =   np.vstack(mGrids)
    np.save(gFname,mGrids)

    outDat  =   p.map(histoGramWap,itertools.izip(files,sNameA))
    hg      =   np.zeros((4,ngrid-1))
    outDat  =   [x for x in outDat if x is not None]
    weights =   ascii.read('reweight_s16a_v4.txt')['weight']
    for i,dat in enumerate(outDat):
        hg+=dat*weights[i]
    for i in range(4):
        hg[i,:]=hg[i,:]/sum(hg[i,:])
    np.save(hgFname,hg)
