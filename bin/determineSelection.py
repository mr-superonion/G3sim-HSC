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

import matplotlib.ticker as ticker
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import matplotlib.lines as mlines
from matplotlib.colors import LogNorm

version =   's18'
obsDir  =   'database/%s/fields/' %version
if not os.path.isdir(obsDir):
    raise RuntimeError("Data catalog directory %s does not exist "
    "or is not a directory" % (obsDir))
outDir =   'database/selection'
if not os.path.exists(outDir):
    os.system('mkdir -p %s' %outDir)
obsAname=   os.path.join(outDir,'scat-%s.fits' %version)

pltDir  =   'plot/selection/'
if not os.path.exists(pltDir):
    os.system('mkdir -p %s' %pltDir)

def main():
    allObs=prepare()
    #onePlot(None,allObs)
    return

def prepare():
    print('Achieving galaxy data from observations')
    if not os.path.exists(obsAname):
        allObs  =   list(map(getSmallCat,utilities.field_keys.values()))
        allObs  =   [x for x in allObs if x is not None]
        allObs  =   vstack(allObs)
        allObs.write(obsAname)
    else:
        print('Already have prepared catalog')
        allObs=Table.read(obsAname)
    return allObs

def getSmallCat(fieldname):
    """
    #load the necessary data from catalog (after FDFC, WL cut)
    """
    obsFname=   os.path.join(obsDir,'%s.fits' %(fieldname))
    data    =   pyfits.getdata(obsFname)
    """
    # FDFC Column for S16A
    """
    data    =   data[utilities.get_FDFC_flag(data,'s16')]
    ra,dec  =   utilities.get_radec(data)

    """
    # Main Column for S18A
    """
    res     =   utilities.get_res(data)
    snr     =   utilities.get_snr(data)
    mag     =   utilities.get_imag(data)
    absE    =   utilities.get_abs_ellip(data)
    logb    =   utilities.get_logb(data)
    fwhm    =   utilities.get_psf_size(data,'fwhm')
    traceR  =   utilities.get_gal_size(data,type='trace')
    detR    =   utilities.get_gal_size(data,type='det')
    traceRP =   utilities.get_psf_size(data,type='trace')
    detRP   =   utilities.get_psf_size(data,type='det')
    absEP   =   utilities.get_abs_ellip_psf(data)
    snrA10,snrA15,snrA20    =   utilities.get_snr_apertures(data)
    snrlb   =   utilities.get_snr_localBG(data)
    logbA,logbR,logbO   =   utilities.get_logbAll(data)
    npass   =   utilities.get_npass(data)
    noivar  =   utilities.get_noi_var(data)
    fieldnameA=[fieldname]*len(data)
    magA10  =   utilities.get_imag_A10(data)
    magA15  =   utilities.get_imag_A15(data)
    magA20  =   utilities.get_imag_A20(data)
    maglb   =   utilities.get_imag_lb(data)
    magpsf  =   utilities.get_imag_psf(data)
    TPid    =   utilities.get_TPid(data)
    isIso   =   utilities.get_isIso(data)
    bsf     =   utilities.get_bs_factor(data)
    obj     =   utilities.get_cmodel_obj(data)
    bmsk18  =   utilities.get_mask_briObj(data)
    objId   =   data['object_id']
    del data
    gc.collect()

    datanames=  ['ra', 'dec', 'res', 'snr', 'mag', 'absE', 'fwhm', \
            'maglb', 'snrlb',  'isIso', 'bsf', 'obj', 'object_id', 'fieldname', \
            'magA10', 'magA15', 'magA20', 'noivar', 'magpsf', \
            'traceR', 'detR', 'traceRP', 'detRP', 'absEP', 'TPid', 'bmsk18', \
            'snrA10', 'snrA15', 'snrA20','logb', 'logbR', 'logbO', 'npass' ]
    dataTab =  [ra, dec, res, snr, mag, absE, fwhm ,\
            maglb, snrlb,  isIso, bsf, obj, objId, fieldnameA, \
            magA10, magA15, magA20, noivar, magpsf, \
            traceR, detR, traceRP, detRP, absEP, TPid, bmsk18, \
            snrA10, snrA15, snrA20,logb, logbR, logbO, npass]
    return Table(data=dataTab,names=datanames)

def prepareColListB(pergroup):
    colns1  =   ['snr','snr16','dTraceR']

    colns2  =   ['res','res16','snr']

    col1List=   utilities.chunkPList(colns1,pergroup)
    col2List=   utilities.chunkPList(colns2,pergroup)
    return col1List,col2List

def onePlot(sims,data):
    sims    =   None
    tryList =  ['std-s16-all','std-s16-good','std-s16-suspect','std-s16-dR','std-s16-magA10']

    data['dTraceR'] =   data['traceR']-data['traceR16']
    for itry in tryList:
        # Basic Mask for both data and sims
        mList   =   itry.split('-')
        if 'iso' in mList[0]:
            mask=   (data['isIso'])
            if sims:
                maskS=  (sims['isIso'])
        elif 'std' in mList[0]:
            mask=   (data['logb']<=-0.375)
            if sims:
                maskS=  (sims['logb']<=-0.375)

        if mList[1] in utilities.field_names.keys():
            mask&=  (data['fieldname'] == mList[2])
            if sims:
                maskS&= (sims['fieldname'] == mList[2])
        elif mList[1]=='s16':
            mask&=  (data['inS16FDFC'])

        mask    &=  (data['npass']>=2)

        if mList[2]=='good':
            mask&=  (data['inS16'])
            cc  =   'Blues'
        elif mList[2]=='suspect':
            mask&=  (~data['inS16'])
            cc  =   'Reds'
        elif mList[2]=='all':
            cc  =   'Greys'
        elif mList[2]=='dR':
            mask&=  (data['dTraceR']<0.12)
            cc  =   'Greys'
        elif mList[2]=='magA10':
            mask&=  (data['magA10']<25.5)
            if sims:
                maskS&=  (sims['magA10']<25.5)
            cc  =   'Greys'

        dataU   =   data[mask]
        if sims:
            simsU   =   sims[maskS]
        else:
            simsU   =   None

        # Do 2D Tests
        pergroup=   3
        col1List,col2List=prepareColListB(pergroup)
        for i in range(len(col1List)):
            pltUtil.cols2DPlot(dataU,col1List[i],col2List[i],pltDir=pltDir,color=cc,label=itry,sims=simsU,ngridH=128,ngridC=64,writeFits=False)
    return

if __name__=='__main__':
    main()
