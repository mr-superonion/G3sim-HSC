#!/usr/bin/env python
import os
import smatch
import fitsio
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


version =   's18-10-12-HECTOMAP'
anaDir  =   'anaCat'
obsAname=   os.path.join(anaDir,'allObs-%s-newMask.fits' %version)
if not os.path.isdir(anaDir):
    raise RuntimeError("Simulation catalog directory %s does not exist "
                 "or is not a directory" % (anaDir))
if not os.path.exists(obsAname):
    raise RuntimeError("required file does not exists")

def main():
    pergroup=   3
    pltDir  =   'plot/galDis-%s/' %version
    if not os.path.exists(pltDir):
        os.system('mkdir -p %s' %pltDir)
    datAname=   os.path.join(anaDir,'allMat16V1-%s-newmask.fits' %version)
    if not os.path.exists(datAname):
        allDat  =   prepareAllDat()
        allDat.write(datAname)
    else:
        print('Already have prepared catalog')
        allDat  =   Table.read(datAname)

    col1List,col2List=prepareColList(pergroup)
    itry    =   'S16v1Mat-all-HECTOMAP'
    mList   =   itry.split('-')
    mask    =   np.ones(len(allDat)).astype(bool)
    if 'suspect'in mList[1]:
        mask   &=   (~allDat['inS16'])
        cc      =   'Reds'
    elif 'all' in mList[1]:
        cc      =   'Blues'
    elif 'good' in mList[1]:
        cc      =   'Blues'
        mask   &=   (allDat['inS16'])

    allDat['dmagPC']=allDat['magpsf']-allDat['mag']
    allDat['dmagAlb']=allDat['maglb']-allDat['magA10']
    for i in range(len(col1List)):
        pltUtil.cols2DPlot(allDat[mask],col1List[i],col2List[i],pltDir=pltDir,color=cc,label=itry,nfig=pergroup)
    return

def prepareAllDat():
    allObs  =   Table.read(obsAname)
    mask    =   (allObs['logb']<=-0.375)
    mask    &=  (allObs['npass']>=2)
    mask    &=  (allObs['inS16FDFC'])
    #mask    &=  (~allObs['mask_briObj_s18'])
    #mask    &=  allObs['inS16M']
    allObs  =   allObs[mask]

    #fieldList=  utilities.field_names.keys()
    fieldList=  ['HECTOMAP']
    allDat =   []
    for field in fieldList:
        maskU   =  (allObs['fieldname'] == field)
        data    =   doMatch(field,allObs[maskU])
        allDat.append(data)
    allDat=    vstack(allDat)
    return allDat

def prepareColList(pergroup):
    colns1  =   []
    colns2  =   []
    colns1.extend(['snr','mag','res'])
    colns2.extend(['snr16','mag16','res16'])

    colns1.extend(['snr','snr','res'])
    colns2.extend(['res','mag','mag'])

    colns1.extend(['snr16','snr16','res16'])
    colns2.extend(['res16','mag16','mag16'])

    colns1.extend(['magA10','mag','snr'])
    colns2.extend(['dmagPC','dmagPC','snrlb'])

    colns1.extend(['magA10','mag','snr'])
    colns2.extend(['dmagAlb','dmagAlb','dmagAlb'])

    colns1.extend(['snrA10','snr','snrA10'])
    colns2.extend(['dmagPC','dmagPC','snrlb'])


    col1List=   utilities.chunkPList(colns1,pergroup)
    col2List=   utilities.chunkPList(colns2,pergroup)
    return col1List,col2List

def doMatch(field,data):
    """
    # Read data
    """
    refDir  =   '/work/xiangchong.li/work/S16AStandard/S16AStandardPre/field'
    refname =   os.path.join(refDir,'%s.fits'%(field))
    data16  =   Table.read(refname)

    """
    # Matching
    """
    ra16    =   data16['ira']
    dec16   =   data16['idec']
    ra      =   data['ra']
    dec     =   data['dec']
    match   =   smatch.match(ra16,dec16,\
            0.5/3600,ra,dec,maxmatch=1)

    """
    # Recording
    """
    data['absEP16'] =   np.zeros(len(data))
    data['res16'][match['i2']]  =   utilities.get_res(data16[match['i1']])
    data['snr16'][match['i2']]  =   utilities.get_snr(data16[match['i1']])
    data['mag16'][match['i2']]  =   utilities.get_imag(data16[match['i1']])
    data['logb16'][match['i2']] =   utilities.get_logb(data16[match['i1']])
    data['absE16'][match['i2']] =   utilities.get_abs_ellip(data16[match['i1']])
    data['absEP16'][match['i2']]=   utilities.get_abs_ellip_psf(data16[match['i1']])
    data['fwhm16'][match['i2']] =   utilities.get_psf_size(data16[match['i1']])
    print(len(match)/len(data))
    return data[match['i2']]

if __name__=='__main__':
    main()
