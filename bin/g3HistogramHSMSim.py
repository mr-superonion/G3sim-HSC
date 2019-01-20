#!/usr/bin/env python
import os
import logging
import itertools
import numpy as np
import astropy.io.fits as pyfits
import astropy.io.ascii as ascii
from multiprocessing import Pool
from astropy.table import Table

def histoGram(ifile):
    dataOut =   np.zeros((4,59))
    rescat =   './outcomeFPFS/src-%04d.fits'%(ifile)
    if not os.path.exists(rescat):
        logging.info('cannot find input file')
        return None 
    try:
        data    =   pyfits.getdata(rescat)
    except:
        logging.info('cannot open file %d' %(ifile))
        return None
        
    Res     =   data['ext_shapeHSM_HsmShapeRegauss_resolution']
    snr     =   data['modelfit_CModel_flux']/data['modelfit_CModel_fluxSigma']
    mask    =   np.ones(len(data)).astype(bool)
    mask    =   mask&(data['deblend_nChild']==0)
    #mask    =   mask&(np.all(~np.isnan(data['fpfs_moments']),axis=1))
    #mask    =   mask&(data['fpfs_moments'][:,0]>600)
    mask    =   mask&(data['Mag']<24.5)
    mask    =   mask&(data['base_ClassificationExtendedness_value']>0)
    mask    =   mask&(~np.isnan(data['ext_shapeHSM_HsmSourceMoments_xx']))
    mask    =   mask&(~np.isnan(data['ext_shapeHSM_HsmSourceMoments_yy']))
    mask    =   mask&(~np.isnan(data['ext_shapeHSM_HsmPsfMoments_xx']))
    mask    =   mask&(~np.isnan(data['ext_shapeHSM_HsmPsfMoments_yy']))
    mask    =   mask&(~np.isnan(data['modelfit_CModel_flux']))
    mask    =   mask&(~np.isnan(data['modelfit_CModel_fluxSigma']))
    mask    =   mask&(~np.isnan(data['Mag']))
    mask    =   mask&(~np.isnan(data['ext_shapeHSM_HsmShapeRegauss_e1']))
    mask    =   mask&(~np.isnan(data['ext_shapeHSM_HsmShapeRegauss_e2']))
    mask    =   mask&(~np.isnan(data['base_Variance_value']))
    mask    =   mask&(data['base_Variance_value']<0.01)
    mask    =   mask&(data['base_Blendedness_abs_flux']<10**(-0.375))
    mask    =   mask&(Res>=0.3)
    mask    =   mask&(snr>=10.)
    dx      =   abs(data['base_SdssCentroid_x']%64-32)
    dy      =   abs(data['base_SdssCentroid_y']%64-32)
    mask    =   mask&(dx**2.+dy**2.<100)
    mask    =   mask&(data['detected']!=0)
    data    =   data[mask]
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
    RA      =   (R1+R2)/2.
    mask2   =   (RA>0.)&(RA<5.)
    data    =   data[mask2]
    """
    absE    =   data['ext_shapeHSM_HsmShapeRegauss_e1']**2.+data['ext_shapeHSM_HsmShapeRegauss_e2']**2.
    absE    =   np.sqrt(absE)
    mag     =   data['Mag']
    ells    =   [Res,snr,mag,absE]
    mGrids  =   np.load('./histGridsHSM.npy')
    for i in range(4):
        dataOut[i,:]=np.histogram(ells[i],bins=mGrids[i])[0]
    return dataOut

if __name__ == '__main__':
    ndat    =   1600
    files   =   range(ndat)
    p       =   Pool(8)
    ngrid   =   60
    gFname  =   './histGridsHSM.npy' 
    hgFname =   './histGramHSMSim1.npy' 
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
    print(len(outDat))
    weiFname=   'noiVar_rPsf_s18a.csv'
    if os.path.exists(weiFname):
        info    =   ascii.read(weiFname)
        weights =   info['weightNoi']*info['weightPsf']
    else:
        weights  =   np.ones(len(outDat))

    for i,dat in enumerate(outDat):
        hg+=dat*weights[i]
    for i in range(4):
        hg[i,:]=hg[i,:]/sum(hg[i,:])
    np.save(hgFname,hg)
