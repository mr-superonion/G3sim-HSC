#!/usr/bin/env python
import os
import logging
import itertools
import numpy as np
import astropy.io.fits as pyfits
import astropy.io.ascii as ascii
import fitsio
from multiprocessing import Pool
from astropy.table import Table

def getHist(ells,mGrids,ngrid):
    dataOut =   np.zeros((4,ngrid-1))
    for i in range(4):
        dataOut[i,:]=np.histogram(ells[i],bins=mGrids[i])[0]
    return dataOut

def histoGram(ifield):
    if ifield=='AEGIS':
        return
    ngrid   =   60
    #Do mask
    directory= '/work/xiangchong.li/work/S16AFPFS/final/noBg/field'
    data    =   fitsio.read(os.path.join(directory,'%s.fits'%(ifield)))
    mask    =   np.ones(len(data)).astype(bool)
    mask    =   mask&(~np.isnan(data['reGaus_resolution']))
    mask    =   mask&(~np.isnan(data['CModel_snr']))
    mask    =   mask&(~np.isnan(data['CModel_mag']))
    mask    =   mask&(~np.isnan(data['reGaus_sigmaE']))
    mask    =   mask&((data['CModel_mag'])<24.5)
    mask    =   mask&(data['extendness']>0)
    mask    =   mask&(data['reGaus_resolution']>0.3)
    mask    =   mask&(data['CModel_snr']>10.)
    data    =   data[mask]
    absE    =   data['reGaus_sigmaE']
    Res     =   data['reGaus_resolution']
    mag     =   data['CModel_mag']
    snr     =   data['CModel_snr']
    ells    =   [Res,snr,mag,absE]
    mGrids  =   []
    mGrids.append(np.linspace(0.1,1.,ngrid))
    mGrids.append(np.linspace(0.,80.,ngrid))
    mGrids.append(np.linspace(22.,25.2,ngrid))
    mGrids.append(np.linspace(0.,2.,ngrid))
    mGrids  =   np.vstack(mGrids)
    return getHist(ells,mGrids,ngrid)

def main():
    ngrid   =   60
    pool    =   Pool(6)
    ngrid   =   60
    fields  =   np.load('./dr1FieldTract.npy').item().keys() 
    outDat  =   pool.map(histoGram,fields)
    outDat  =   [x for x in outDat if x is not None]
    if len(outDat)<2:
        return
    hg      =   np.zeros((4,ngrid-1))
    for dat in outDat:
        hg+=dat
    for i in range(4):
        hg[i,:]/= sum(hg[i,:])
    np.save('histGramHSMObs.npy',hg)
    return

if __name__ == '__main__':
    main()
