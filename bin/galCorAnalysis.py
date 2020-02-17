#!/usr/bin/env python
import os
import treecorr
import numpy as np
from astropy.table import Table, vstack
import astropy.io.fits as pyfits
import sim_analysis_utilities as utilities
import plot_catalog_utilities as pltUtil

import matplotlib
matplotlib.use('agg')
matplotlib.rcParams['ps.useafm'] = True
matplotlib.rcParams['pdf.use14corefonts'] = True
matplotlib.rcParams['text.usetex'] = True

import matplotlib.ticker as ticker
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import matplotlib.lines as mlines
from matplotlib.colors import LogNorm

pltDir  =   './s18-2/plot/galCor/'
if not os.path.exists(pltDir):
    os.system('mkdir -p %s' %pltDir)

version =   's18'
obsDir  =   './database/compare/'
if not os.path.isdir(obsDir):
    raise RuntimeError("Data catalog directory %s does not exist "
         "or is not a directory" % (obsDir))

def cor2PNK(catN,catK,pltDir,nbins=30,min_sep=0.05,max_sep=2.,unit='deg',label='base'):
    ddFname =   os.path.join(pltDir,'cor2PNK-%s.csv' %label)
    if not os.path.exists(ddFname):
        corDD   =   treecorr.NKCorrelation(nbins=nbins,min_sep=min_sep,max_sep=max_sep,sep_units=unit)
        corDD.process(catN,catK)
        corDD.write(ddFname)
    return

def getWLMask(data):
    mask    =   (data['npass']>=2)
    mask    &=  (data['logb']<=-0.375)
    mask    &=  (data['res']>0.3)
    mask    &=  (data['snr']>10.)
    mask    &=  (data['mag']<24.5)
    return mask

def getUnDetData():
    inFname =   os.path.join(obsDir,'scat-s16-2.fits')
    if not os.path.exists(inFname):
        preFname=   os.path.join(obsDir,'scat-s16.fits')
        data    =   Table.read(preFname)
        inS18   =   np.zeros(len(data),dtype=bool)
        inS19   =   np.zeros(len(data),dtype=bool)
        info    =   pyfits.getdata(os.path.join(obsDir,'match-s18.fits'))
        inS18[info['scat_id16']]=True
        info    =   pyfits.getdata(os.path.join(obsDir,'match-s19.fits'))
        inS19[info['scat_id16']]=True
        data['inS18']=inS18
        data['inS19']=inS19
        data.write(inFname)
    else:
        data=pyfits.getdata(inFname)
    return data

def getDetData():
    inFname =   os.path.join(obsDir,'scat-%s.fits' %version)
    data    =   pyfits.getdata(inFname)
    return data

def getMatData():
    inFname =   os.path.join(obsDir,'match-%s.fits' %version)
    data    =   pyfits.getdata(inFname)
    mask    =   (~np.isnan(data['mag']))
    mask   &=   (~np.isnan(data['mag16']))
    mask   &=   ((data['mag16']>-100))
    mask   &=   ((data['mag16']<100))
    mask   &=   ((data['mag']>-100))
    mask   &=   ((data['mag']<100))
    return data[mask]

def mainA():
    label   =   'wl-star-mag'
    if label.split('-')[2]=='undetect':
        data    =   getUnDetData()
        kappa   =   (~data['in%s' %version.upper()]).astype(int)
    elif label.split('-')[2]=='mag':
        data    =   getMatData()
        kappa   =   data['mag16']-data['mag']
    elif label.split('-')[2]=='detect':
        data    =   getDetData()
        kappa   =   np.ones(len(data))

    if label.split('-')[0]=='wl':
        mask    =   getWLMask(data)
        data    =   data[mask]
        kappa   =   kappa[mask]

    catK    =   treecorr.Catalog(ra=data['ra'],dec=data['dec'],k=kappa,ra_units='deg',dec_units='deg')

    if label.split('-')[1]=='star':
        bsAname =   '/work/xiangchong.li/data/brightObjMask/Arcturus/stars/stars_all.fits'
        bsCat   =   pyfits.getdata(bsAname,ext=1)
    elif label.split('-')[1]=='CAMIRA':
        bsAname =   '/work/xiangchong.li/work/lensCatalog/CAMIRA/CAMIRA_S16A.csv'
        bsCat   =   Table.read(bsAname)
        bsCat['coord_ra'].name='ra'
        bsCat['coord_dec'].name='dec'

    catN    =   treecorr.Catalog(ra=bsCat['ra'],dec=bsCat['dec'],ra_units='deg',dec_units='deg')
    cor2PNK(catN,catK,pltDir,min_sep=0.1,max_sep=20.,unit='arcmin',label=label+'-%s'%version)
    return
if __name__=='__main__':
    mainA()
