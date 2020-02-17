#!/usr/bin/env python
import os
import gc
import smatch
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
matplotlib.rcParams['font.size'] = 14

import matplotlib.ticker as ticker
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import matplotlib.lines as mlines
from matplotlib.colors import LogNorm

version =   's19'
cmpDir  =   'database/compare/'
catname =   os.path.join(cmpDir,'scat-%s.fits' %version)
matname =   os.path.join(cmpDir,'match-%s.fits' %version)
if not os.path.exists(catname):
    raise RuntimeError("Data catalog directory %s does not exist "
    "or is not a directory" % (catname))
cat16name= os.path.join(cmpDir,'scat-s16.fits' )

pltDir  =   's18-2/plot/compare/'
if not os.path.exists(pltDir):
    os.system('mkdir -p %s' %pltDir)

def matchT16():
    if not os.path.exists(matname):
        data    =   Table.read(catname)
        data16  =   Table.read(cat16name)
        match   =   smatch.match(data16['ra'],data16['dec'],0.4/3600.,\
            data['ra'],data['dec'],maxmatch=1)
        index   =   np.unique(match['i1'],return_index=True)[1]
        match   =   match[index]
        data    =   data[match['i2']]
        data16  =   data16[match['i1']]
        data['object_id16']=   data16['object_id']
        data['snr16']   =   data16['snr']
        data['res16']   =   data16['res']
        data['mag16']   =   data16['mag']
        data['absE16']  =   data16['absE']
        data['traceR16']=   data16['traceR']
        data['TPid16']  =   data16['TPid']
        data['scat_id16']=  match['i1']
        data['scat_id'] =   match['i2']
        data.write(matname)
    else:
        data    =   Table.read(matname)
    return data

def prepareColListB(pergroup):
    colns1=[]
    colns2=[]
    colns1.extend(['snr', 'snr','snr'])
    colns2.extend(['res', 'mag','absE'])

    col1List=   utilities.chunkPList(colns1,pergroup)
    col2List=   utilities.chunkPList(colns2,pergroup)
    return col1List,col2List

def prepareColListA(pergroup):
    colns1=[]
    colns2=[]
    colns1.extend(['dmag'])
    colns2.extend(['mag'])
    """
    colns1.extend(['obj','obj','obj'])
    colns2.extend(['snr','mag','magA10'])
    """

    col1List=   utilities.chunkPList(colns1,pergroup)
    col2List=   utilities.chunkPList(colns2,pergroup)
    return col1List,col2List

def onePlot(data,pergroup,pnd,col1List,col2List,cc='Greys'):
    label   =   version+pnd
    for i in range(len(col1List)):
        colns1=col1List[i]
        colns2=col2List[i]
        fig     =   plt.figure(figsize=(5*pergroup, 5))
        axes    =   []
        for j in range(pergroup):
            axes.append(fig.add_subplot(1,pergroup,j+1))
        pltUtil.cols2DPlot(data,colns1,colns2,axes,pltDir=pltDir,color=cc,label=label,ngridH=128)
        outFname=   os.path.join(pltDir, 'galDis2D-%s-%s_%s.png' %(colns1[0],colns2[0],label))
        fig.savefig(outFname)
    return

def mainA(matDat):
    # Prepare column names to plot
    pergroup=   1
    col1List,col2List=prepareColListA(pergroup)
    matDat['dTraceR'] =   matDat['traceR']-matDat['traceR16']
    matDat['dmag']    =   matDat['mag']-matDat['mag16']
    matDat['dsnr']    =   matDat['snr']-matDat['snr16']

    # Plot
    pndL    =   ['-all-all']
    for pnd in pndL:
        if pnd.split('-')[1]=='wl':
            mask    =   (matDat['npass']>=2)
            mask    &=  (matDat['logb']<=-0.375)
            mask    &=  (matDat['res']>0.3)
            mask    &=  (matDat['snr']>10.)
        if pnd.split('-')[1]=='all':
            mask    =   np.ones(len(matDat),dtype=bool)

        if pnd.split('-')[2]=='good':
            mask    &=  (matDat['mag16']<24.5)
            cc  =   'Blues'
        elif pnd.split('-')[2]=='bad':
            mask    &=  (matDat['mag16']>24.5)
            cc  =   'Reds'
        elif pnd.split('-')[2]=='all':
            cc   =  'Greys'

        dataU  =   matDat[mask]
        onePlot(dataU,pergroup,pnd,col1List,col2List,cc=cc)
    return

def mainB(matDat,sims=None):
    # Prepare column names to plot
    pergroup=   3
    col1List,col2List=prepareColListB(pergroup)
    matDat['dTraceR'] =   matDat['traceR']-matDat['traceR16']
    matDat['dmag']    =   matDat['mag']-matDat['mag16']
    matDat['dsnr']    =   matDat['snr']-matDat['snr16']
    # Plot
    pnd     =   '-wl'
    if pnd.split('-')[1]=='wl':
        mask    =   (matDat['npass']>=2)
        mask    &=  (matDat['logb']<=-0.375)
        mask    &=  (matDat['res']>0.3)
        mask    &=  (matDat['snr']>10.)
        matDat  =   matDat[mask]
        if sims is not None:
            maskS   =   (sims['logb']<=-0.375)
            sims    =   sims[maskS]
    cn  =   'magA10'
    for maglm in np.arange(26.0,25.0,-0.1):
        pndU    =   pnd+'-%s-%.1f' %(cn,maglm)
        dataU   =   matDat[matDat[cn]<maglm]
        simsU   =   sims[sims[cn]<maglm]
        onePlot(dataU,pergroup,'obs'+pndU,col1List,col2List)
        if sims is not None:
            onePlot(simsU,pergroup,'sim'+pndU,col1List,col2List)
    return

if __name__=='__main__':
    # Prepare data to plot
    matDat  =   matchT16()
    mainA(matDat)
    """
    sims    =   Table.read('s18-2/anaCat/allSim.fits')
    mainB(matDat,sims)
    """
