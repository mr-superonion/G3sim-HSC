#!/usr/bin/env python
import os
import smatch
import catSelect
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

maxMag  =   20
anaDir  =   'anaCat'
version =   's16'
obsDir  =   './database/%s/fields/' %version
if not os.path.isdir(obsDir):
    raise RuntimeError("Data catalog directory %s does not exist "
         "or is not a directory" % (obsDir))
obsDir0 =   './database/'
boAname =   os.path.join(obsDir0,'allBri-%s.fits' %version)
rAname  =   os.path.join(obsDir0,'allRan-%s.fits' %version)

pltDir  =   './s18-2/plot/galCor-%s/' %version
if not os.path.exists(pltDir):
    os.system('mkdir -p %s' %pltDir)

def main():
    allBri,allRan=getAllObj()
    return

def getAllObj():
    print('Achieving bright objects')
    if not os.path.exists(boAname):
        allBri  =   list(map(getBriObj,utilities.field_keys.values()))
        allBri  =   [x for x in allBri if x is not None]
        allBri  =   vstack(allBri)
        allBri.write(boAname)
    else:
        print('Already have prepared catalog')
        allBri=pyfits.getdata(boAname)

    print('Achieving random points')
    if not os.path.exists(rAname):
        allRan  =   list(map(getRanObj,utilities.field_keys.values()))
        allRan  =   [x for x in allRan if x is not None]
        allRan  =   vstack(allRan)
        allRan.write(rAname)
    else:
        print('Already have prepared catalog')
        allRan=pyfits.getdata(rAname)
    return allBri,allRan

def getBriObj(fieldname):
    datanames   =   ['ra','dec','mag','field','isStar']
    print('Processing field %s' %fieldname)
    print('Achieving bright galaxies')
    galFname=   os.path.join(obsDir,'%s.fits' %(fieldname))
    datG    =   pyfits.getdata(galFname)
    raG,decG=   utilities.get_radec(datG)
    magG    =   utilities.get_imag(datG)
    fnG     =   [fieldname]*len(datG)
    isStarG =   np.zeros(len(datG),dtype=bool)
    datG    =   Table(data=[raG,decG,magG,fnG,isStarG]
                ,names=datanames)
    datG    =   datG[magG<maxMag]

    print('Achieving bright stars')
    staFname=   os.path.join(obsDir,'%s_stars.fits' %(fieldname))
    datS    =   pyfits.getdata(staFname)
    raS,decS=   utilities.get_radec(datS)
    magS    =   utilities.get_imag_psf(datS)
    fnS     =   [fieldname]*len(datS)
    isStarS =   np.ones(len(datS),dtype=bool)
    datS    =   Table(data=[raS,decS,magS,fnS,isStarS]
                ,names=datanames)
    datS    =   datS[magS<maxMag]
    return vstack([datG,datS])

def getBriObj(fieldname):
if __name__=='__main__':
    main()
