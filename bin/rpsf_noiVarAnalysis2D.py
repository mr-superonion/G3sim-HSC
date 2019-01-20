#!/usr/bin/env python
"""
Script to loop over simulated subfields and check the noise variance
estimated by the HSC pipeline, as compared to the true one.
"""
import os
import fitsio
import astropy.io.fits as pyfits
import numpy as np
import astropy.io.ascii as ascii
from multiprocessing import Pool
from astropy.table import Table
import matplotlib
matplotlib.use('agg')
matplotlib.rcParams['ps.useafm'] = True
matplotlib.rcParams['pdf.use14corefonts'] = True
matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt


def getHistSim(data):
    ngrid   =   20
    fwhms   =   data['rpsf_det']*0.168*2.355
    varNois =   data['var_meas']
    minVar  =   varNois.min()
    gridNoi =   np.linspace(minVar,0.01,ngrid)
    histNoi =   np.histogram(varNois,bins=gridNoi,density=True)[0]
    histNoi =   histNoi/sum(histNoi)
    gridPsf =   np.linspace(0.3,0.9,ngrid)
    histPsf =   np.histogram(fwhms,bins=gridPsf,density=True)[0]
    histPsf =   histPsf/sum(histPsf)
    outcome =   {'gridNoi':gridNoi,'gridPsf':gridPsf,'histNoiSim':histNoi,'histPsfSim':histPsf}
    np.save('noi_psf_histGram.npy',outcome)
    return

def getHistObs():
    pool    =   Pool(6)
    fields  =   np.load('./dr1FieldTract.npy').item().keys() 
    outDat  =   pool.map(getHistObsSub,fields)
    outDat  =   [x for x in outDat if x is not None]
    ngrid   =   20
    if len(outDat)<2:
        return
    hg      =   np.zeros((2,ngrid-1))
    for dat in outDat:
        hg+=dat
    for i in range(2):
        hg[i,:]/= sum(hg[i,:])
    histNoi =   hg[0]
    histPsf =   hg[1]
    histInfo=   np.load('noi_psf_histGram.npy').item()
    histInfo['histNoiObs']=histNoi
    histInfo['histPsfObs']=histPsf
    np.save('noi_psf_histGram.npy',histInfo)
    return histInfo

def getHistObsSub(ifield):
    if ifield=='AEGIS':
        return
    ngrid   =   20
    #Do mask
    directory= '/work/xiangchong.li/work/S16AFPFS/final/noBg/field'
    data    =   fitsio.read(os.path.join(directory,'%s.fits'%(ifield)))
    mask    =   np.ones(len(data)).astype(bool)
    mask    =   mask&(~np.isnan(data['noiVar']))
    mask    =   mask&(~np.isnan(data['psfMxx']))
    mask    =   mask&(~np.isnan(data['psfMyy']))
    mask    =   mask&(~np.isnan(data['psfMxy']))
    data    =   data[mask]
    fwhms   =   (data['psfMxx']*data['psfMyy'])**(0.25)*0.168*2.355
    varNois =   data['noiVar']
    gridInfo=   np.load('noi_psf_histGram.npy').item()
    gridNoi =   gridInfo['gridNoi']
    gridPsf =   gridInfo['gridPsf']
    histPsf =   np.histogram(fwhms,bins=gridPsf,density=True)[0]
    histNoi =   np.histogram(varNois,bins=gridNoi,density=True)[0]
    return [np.vstack([histNoi,histPsf]),hist2d]


if __name__ == '__main__':
    outFname=   'noiVar_rPsf_s18a.csv'
    t   =   Table.read(outFname)
    getHistSim(t)
    histInfo=   getHistObs()
     
