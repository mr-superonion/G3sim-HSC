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


def getBrightObjCat():
    bsAname =   '/work/xiangchong.li/data/brightObjMask/Arcturus/stars/stars_all.fits'
    bsCat1  =   Table.read(bsAname)
    bsCat1  =   bsCat1[['ra','dec']]

    bsAname =   '/work/xiangchong.li/work/lensCatalog/CAMIRA/CAMIRA_S16A.csv'
    bsCat2  =   Table.read(bsAname)
    bsCat2['coord_ra'].name='ra'
    bsCat2['coord_dec'].name='dec'
    bsCat2  =   bsCat2[['ra','dec']]

    outDat  =   vstack([bsCat1,bsCat2])
    outName =   'database/brightObjs.fits'
    outDat.write(outName)
    return

if __name__=='__main__':
    getBrightObjCat()

