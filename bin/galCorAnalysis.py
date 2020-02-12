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

def cor2PNN(datA,datB,datR,pltDir,nbins=30,min_sep=0.05,max_sep=2.,unit='deg',label='base'):
    catA    =   treecorr.Catalog(ra=datA['ra'],dec=datA['dec'],ra_units='deg',dec_units='deg')
    catB    =   treecorr.Catalog(ra=datB['ra'],dec=datB['dec'],ra_units='deg',dec_units='deg')
    catR    =   treecorr.Catalog(ra=datR['ra'],dec=datR['dec'],ra_units='deg',dec_units='deg')

    ddFname =   os.path.join(pltDir,'cor2PNN-DD-%s.csv' %label)
    if not os.path.exists(ddFname):
        corDD   =   treecorr.NNCorrelation(nbins=nbins,min_sep=min_sep,max_sep=max_sep,sep_units=unit)
        corDD.process(catA,catB)
        corDD.write(ddFname)
    else:
        corDD   =   Table.read(ddFname)

    drFname =   os.path.join(pltDir,'cor2PNN-DR-%s.csv' %label)
    if not os.path.exists(drFname):
        corDR   =   treecorr.NNCorrelation(nbins=nbins,min_sep=min_sep,max_sep=max_sep,sep_units=unit)
        corDR.process(catB,catR)
        corDR.write(drFname)
    else:
        corDR   =   Table.read(drFname)

    rrFname =   os.path.join(pltDir,'cor2PNN-RR.csv')
    if not os.path.exists(rrFname):
        corRR   =   treecorr.NNCorrelation(nbins=nbins,min_sep=min_sep,max_sep=max_sep,sep_units=unit)
        corRR.process(catR,catR)
        corRR.write(rrFname)
    else:
        corRR   =   Table.read(rrFname)
    return

def prepareRan(version):
    obsDir  =   'database/%s-rand/fields/' %version
    if not os.path.isdir(obsDir):
        raise RuntimeError("Data catalog directory %s does not exist "
        "or is not a directory" % (obsDir))
    cmpDir =   'database/compare/'
    scatname=   os.path.join(cmpDir,'scat-%s-rand.fits' %version)
    if not os.path.exists(scatname):
        data  = []
        for fn in utilities.field_keys.values():
            print('Reading data from field %s' %fn)
            ranFname=  os.path.join(obsDir,'%s_rand.fits' %fn)
            tmp =   Table.read(ranFname)
            tmp['inS16FDFC'] =   utilities.get_FDFC_flag(tmp,'s16')
            data.append(tmp)
        data  =   vstack(data)
        data.write(scatname)
    else:
        print('Already have prepared catalog')
        data=Table.read(scatname)
    return data

def mainA():
    version =   's18'
    obsDir  =   './database/compare/'
    if not os.path.isdir(obsDir):
        raise RuntimeError("Data catalog directory %s does not exist "
             "or is not a directory" % (obsDir))
    #boAname =   os.path.join(obsDir,'scat-s16-brightGal.fits' )
    boAname =   '/work/xiangchong.li/data/brightObjMask/Arcturus/stars/stars_all.fits'
    inFname =   os.path.join(obsDir,'scat-s16-unmatch-%s.fits' %version)

    boCat   =   pyfits.getdata(boAname,ext=1)
    data    =   pyfits.getdata(inFname)

    pltDir  =   './s18-2/plot/galCor/'
    if not os.path.exists(pltDir):
        os.system('mkdir -p %s' %pltDir)
    ranCat  =   prepareRan(version)
    ranCat  =   ranCat[ranCat['inS16FDFC']]
    cor2PNN(data,boCat,ranCat,pltDir)
    return

def mainB():
    version =   's18'
    obsDir  =   './database/compare/'
    if not os.path.isdir(obsDir):
        raise RuntimeError("Data catalog directory %s does not exist "
             "or is not a directory" % (obsDir))
    #boAname =   os.path.join(obsDir,'scat-s16-brightGal.fits' )
    boAname =   '/work/xiangchong.li/data/brightObjMask/Arcturus/stars/stars_all.fits'
    inFname =   os.path.join(obsDir,'scat-s16-unmatch-%s.fits' %version)

    boCat   =   pyfits.getdata(boAname,ext=1)
    data    =   pyfits.getdata(inFname)

    pltDir  =   './s18-2/plot/galCor/'
    if not os.path.exists(pltDir):
        os.system('mkdir -p %s' %pltDir)
    ranCat  =   prepareRan(version)
    ranCat  =   ranCat[ranCat['inS16FDFC']]
    cor2PNN(data,boCat,ranCat,pltDir,min_sep=0.5,max_sep=30.,unit='arcmin')
    return
if __name__=='__main__':
    mainA()
