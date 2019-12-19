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

version =   's18-10-12'
anaDir  =   'anaCat'
simAname=   os.path.join(anaDir,'allSim-%s.fits' %version)
obsAname=   os.path.join(anaDir,'allObs-%s.fits' %version)
obsDir  =   'database/%s/fields/' %version
datanames   =   ['res','snr','mag','absE','logb','weight']
pltDir  =   'plot/galDis-%s/' %version


# Checks that the paths provided are valid
if not os.path.isdir(anaDir):
    raise RuntimeError("Simulation catalog directory %s does not exist "
                 "or is not a directory" % (anaDir))
if not os.path.isdir(obsDir):
    raise RuntimeError("Data catalog directory %s does not exist "
                 "or is not a directory" % (obsDir))
if not os.path.exists(pltDir):
    os.system('mkdir -p %s' %pltDir)

def main():
    if not os.path.exists(obsAname):
        allSim,allObs   =   prepareAllDat()
        allSim.write(simAname)
        allObs.write(obsAname)
    else:
        print('Already have prepared catalog')
        allSim=pyfits.getdata(simAname)
        allObs=pyfits.getdata(obsAname)

    """
        tryList  =  ['suspect_newmask-std-HECTOMAP','s16_newmask-std-HECTOMAP',\
                's16M_newmask-std-HECTOMAP','nomatch_newmask-std-HECTOMAP']
    for itry in tryList:
        # Basic Mask for both data and sims
        mList    =   itry.split('-')
        if 'iso' in mList[1]:
            mask0   =   (allObs['logb']<-5.)
        elif 'std' in mList[1]:
            mask0   =   (allObs['logb']<=-0.375)
        if mList[2] in utilities.field_names.keys():
            mask0   &=  (allObs['fieldname'] == mList[2])

        # Mask only for data
        mask    =  mask0&(allObs['npass']>=2)
        mask    &= (allObs['inS16FDFC'])
        if 'suspect_'in mList[0]:
            mask    &=   (~allObs['inS16'])
            mask    &=   allObs['inS16M']
            cc      =   'Reds'
        elif 's16_' in mList[0]:
            mask    &=   (allObs['inS16'])
            cc      =   'Blues'
        elif 's16M_' in mList[0]:
            mask    &=   (allObs['inS16M'])
            cc      =   'Blues'

        if 'newmask' in mList[0]:
            mask    &= (allObs['mask_briObj_s18']==0)

        # Do 2D Tests
        dataU       =   allObs[mask]
        pergroup    =   3
        col1List,col2List=prepareColListA(pergroup)
        for i in range(len(col1List)):
            pltUtil.cols2DPlot(dataU,col1List[i],col2List[i],pltDir=pltDir,color=cc,label=itry,nfig=pergroup)

    """

    finalPlot(allSim,allObs)
    return

def prepareColListA(pergroup):
    colns1  =   ['snr','mag','res']*3
    colns1.extend(['logb','snrA10','snrA10'])
    colns1.extend(['snr','snr','snr'])

    colns2  =   ['logb','logb','logb']
    colns2.extend(['snrA10','snrA10','snrA10'])
    colns2.extend(['snrlb','snrlb','snrlb'])
    colns2.extend(['snrlb','snrlb','logb'])
    colns2.extend(['res','mag','absE'])

    #colns2.extend(['logbR','logbR','logbR'])
    #colns2.extend(['logbO','logbO','logbO'])
    col1List=   utilities.chunkPList(colns1,pergroup)
    col2List=   utilities.chunkPList(colns2,pergroup)
    return col1List,col2List

def prepareAllDat():
    print('Achieving data')
    outcomes    =   list(map(process,utilities.field_keys.values()))
    outcomes    =   [x for x in outcomes if x is not None]
    allObsList  =   []
    allSimList  =   []
    for outcome in outcomes:
        allObsList.append(outcome[0])
        allSimList.append(outcome[1])
    allSim      =   vstack(allSimList)
    allObs      =   vstack(allObsList)
    del outcomes
    del allObsList
    del allSimList
    return allSim,allObs

def process(fieldname):
    print('Processing field: %s' %(fieldname))

    print('Achieving galaxy data from observations')
    obsFname    =   os.path.join(obsDir,'%s.fits' %(fieldname))
    obsDat      =   getObsGalData(obsFname,fieldname)

    print('Achieving galaxy data from simulations')
    anaFieldDir =   os.path.join(anaDir,'%s' %(fieldname))
    simFname    =   os.path.join(anaFieldDir,'%s.fits' %(fieldname))
    simDat      =   getSimGalData(simFname)
    return [obsDat,simDat]

def getObsGalData(obsFname,fieldname):
    """
    #load the necessary data from catalog (after FDFC, WL cut)
    """
    data    =   pyfits.getdata(obsFname)
    """
    # Match to S15B
    """
    obsFname16= os.path.join(os.environ['homeWrk'],\
            'S16AStandard/S16AStandardCalibrated/field',\
            '%s_calibrated.fits' %fieldname)

    data16  =   pyfits.getdata(obsFname16)
    ra      =   data['i_ra']
    dec     =   data['i_dec']
    """
    # Match to S16
    """
    matchL  =   smatch.match(data16['ira'],data16['idec'],20./3600.,\
            ra,dec,maxmatch=1)
    inS16L  =   np.zeros(len(data)).astype(bool)
    inS16L[matchL['i2']]  =   True
    mDist   =   np.arccos(matchL['cosdist'])*180./np.pi

    matchM  =  matchL[mDist<10./3600.]
    inS16M  =   np.zeros(len(data)).astype(bool)
    inS16M[matchM['i2']]  =   True

    matchS  =   matchL[mDist<5./3600.]
    inS16S  =   np.zeros(len(data)).astype(bool)
    inS16S[matchS['i2']]  =   True

    match   =   matchL[mDist<0.5/3600.]
    inS16   =   np.zeros(len(data)).astype(bool)
    inS16[match['i2']]  =   True

    """
    # Main Column for S15B
    """
    res16   =   np.zeros(len(data))
    res16[match['i2']]  =   utilities.get_res(data16[match['i1']])
    snr16   =   np.zeros(len(data))
    snr16[match['i2']]  =   utilities.get_snr(data16[match['i1']])
    mag16   =   np.zeros(len(data))
    mag16[match['i2']]  =   utilities.get_imag(data16[match['i1']])
    logb16   =   np.zeros(len(data))
    logb16[match['i2']]  =   utilities.get_logb(data16[match['i1']])
    absE16  =   np.zeros(len(data))
    absE16[match['i2']] =   utilities.get_abs_ellip(data16[match['i1']])
    absEP16  =   np.zeros(len(data))
    absEP16[match['i2']] =   utilities.get_abs_ellip_psf(data16[match['i1']])
    fwhm16   =   np.zeros(len(data))
    fwhm16[match['i2']]  =   utilities.get_psf_size(data16[match['i1']])

    """
    # FDFC Column for S15B
    """
    inS16FDFC   =   catSelect.getFDFCS16(data,'i_ra','i_dec')

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
    maglb   =   utilities.get_imag_lb(data)

    """
    # Flags for S18A
    """
    badPsf  =   data['i_pixelflags_inexact_psf']
    badPsfC =   data['i_pixelflags_inexact_psfcenter']
    objId   =   data['object_id']
    weight  =   np.ones(len(data))
    datanmA =   ['fwhm', 'badPsf', 'badPsfC', \
            'inS16', 'inS16L', 'inS16FDFC', 'absEP16', 'maglb', \
            'res16', 'snr16', 'mag16', 'absE16', 'logb16', 'fwhm16', \
            'ra', 'dec', 'object_id', 'fieldname', 'magA10', 'noivar', \
            'traceR', 'detR', 'traceRP', 'detRP', 'absEP', 'inS16M', 'inS16S',\
            'snrA10', 'snrA15', 'snrA20', 'snrlb', 'logbR', 'logbO', 'npass' ]
    dataTab =   [res, snr, mag, absE, logb, weight, fwhm, badPsf, badPsfC, \
            inS16, inS16L, inS16FDFC, absEP16, maglb, \
            res16, snr16, mag16, absE16, logb16, fwhm16, \
            ra, dec, objId, fieldnameA, magA10, noivar, \
            traceR, detR, traceRP, detRP, absEP, inS16M, inS16S, \
            snrA10, snrA15, snrA20, snrlb, logbR, logbO, npass]
    del data
    return Table(data=dataTab,names=datanames+datanmA)

def getSimGalData(simFname):
    data    =   pyfits.getdata(simFname)
    data    =   utilities.update_wl_flags_sim(data)
    mask    =   data['weak_lensing_flag']
    """
    mask    =   mask&(utilities.get_snr(data) < 80.)
    """
    data    =   data[mask]

    absE    =   utilities.get_abs_ellip(data)
    res     =   utilities.get_res(data)
    mag     =   utilities.get_imag(data)
    snr     =   utilities.get_snr(data)
    logb    =   utilities.get_logb(data)
    weight  =   data['reweight_s16a']
    snrA10,snrA15,snrA20    =   utilities.get_snr_apertures(data)
    snrlb   =   utilities.get_snr_localBG(data)
    logbA,logbR,logbO   =   utilities.get_logbAll(data)
    noivar  =   utilities.get_noi_var(data)
    magA10  =   utilities.get_imag_A10(data)
    maglb   =   utilities.get_imag_lb(data)

    dataTab =   [res, snr, mag, absE, logb, weight, \
                snrA10, snrA15, snrA20, snrlb, \
                logbR, logbO, noivar, magA10, maglb]
    datanmA =   ['snrA10', 'snrA15', 'snrA20', 'snrlb',\
                'logbR', 'logbO','noivar','magA10', 'maglb']
    return Table(data=dataTab,names=datanames+datanmA)

def compareS16(data):
    mask    =   data['bO_S18_r'].astype(bool)
    mask    &=  data['inS16']
    data    =   data[mask]
    fig2    =   plt.figure(figsize=(15, 10))

    """
    # SNR scatter
    """
    ax = fig2.add_subplot(2,3,1)
    ax.scatter(data['snr'],data['snr16'],s=0.1)
    ax.set_xlabel('SNR S18')
    ax.set_ylabel('SNR S16')

    """
    # Mag scatter
    """
    ax = fig2.add_subplot(2,3,2)
    ax.scatter(data['mag'],data['mag16'],s=0.1)
    ax.set_xlabel('mag S18')
    ax.set_ylabel('mag S16')

    """
    # logb scatter
    """
    ax = fig2.add_subplot(2,3,3)
    ax.scatter(data['logb'],data['logb16'],s=0.1)
    ax.set_xlabel('logb S18')
    ax.set_ylabel('logb S16')

    """
    # Resolution scatter
    """
    ax = fig2.add_subplot(2,3,4)
    ax.scatter(data['res'],data['res16'],s=0.1)
    ax.set_xlabel(r'$R_2$ S18')
    ax.set_ylabel(r'$R_2$ S16')

    """
    # FWHM scatter
    """
    ax = fig2.add_subplot(2,3,5)
    ax.scatter(data['fwhm'],data['fwhm16'],s=0.1)
    ax.set_xlabel('fwhm S18 (arcsec)')
    ax.set_ylabel('fwhm S16 (arcsec)')

    """
    # absE scatter
    """
    ax = fig2.add_subplot(2,3,6)
    ax.scatter(data['absE'],data['absE16'],s=0.1)
    ax.set_xlabel(r'$|e|$ S18')
    ax.set_ylabel(r'$|e|$ S16')
    """
    # Save figure
    """
    plt.tight_layout()
    fig2.savefig(os.path.join(pltDir,'compareS16_6Obs_newBriObjMsk.png'))
    return

def logbNumPlot(data):
    """
    # Plot the galaxy  number loss vs. logb
    """
    fig =   plt.figure(figsize=(6, 6))
    ax  =   fig.add_subplot(1,1,1)

    print('beginning logbNumPlot')
    data    =   data[data['bO_S18_r'].astype(bool)]
    """
    # Prepare
    """
    numAll  =   len(data)
    thres   =   np.arange(-1.2,-0.375,0.05)
    frac    =   np.array([np.sum(data['logb']<t) for t in thres])/numAll
    ax.plot(thres,frac)
    ax.set_xlabel('log10(b)')
    ax.set_ylabel('fraction')
    ax.grid()
    """
    # Plot
    """
    plt.tight_layout()
    fig.savefig(os.path.join(pltDir,'frac-logb_newBriObjMsk.png'))
    return

def finalPlot(sims,data):
    """
    # Plot the snr - (res,mag.|e|) Joint distributions
    """

    """
    # Do mask
    # Data
    """
    nbf     =   len(data)
    mask    =  (data['npass']>=2)
    print(np.sum(mask)/nbf)
    mask    &= (data['logb']<=-0.375)
    print(np.sum(mask)/nbf)
    mask    &= (data['snrlb']<0.04*(data['snrA10']-10)-0.1) #T1
    #mask    &= (data['snrlb']<0.01*(data['snrA10'])) #T2
    print(np.sum(mask)/nbf)
    data    =   data[mask]
    naf     =   len(data)
    print(naf/nbf)
    """
    # Simulations
    """
    maskS   =   (sims['logb']<=-0.375)
    maskS   &=  (sims['snrlb']<0.04*(sims['snrA10']-10)-0.1)
    sims    =   sims[maskS]

    fig2 = plt.figure(figsize=(15, 5))

    """
    # The First Panel
    """
    ax = fig2.add_subplot(1,3,1)
    counts1= np.histogram2d(data['snr'], data['res'], 32, weights=data['weight'], range=[[10,80],[0.3,1]], normed=True)[0]
    c = ax.contour(counts1.T,extent=[10,80,0.3,1],cmap='viridis', levels=[0.01, 0.025, 0.05, 0.075])
    ax.hist2d(data['snr'], data['res'], 128, weights=data['weight'], range=[[10,80],[0.3,1]], cmap='gray_r', normed=True,vmax=0.08,alpha=0.8);

    counts2= np.histogram2d(sims['snr'], sims['res'], 32, weights=sims['weight'], range=[[10,80],[0.3,1]], normed=True)[0]
    c = ax.contour(counts2.T,extent=[10,80,0.3,1],cmap='viridis',linestyles='dashed', levels=[0.01, 0.025, 0.05, 0.075])
    ax.set_xlabel('CModel SNR')
    ax.set_ylabel('reGaus resolution')
    #ax.set_xscale('log')
    ax.xaxis.set_major_formatter(mtick.FormatStrFormatter('%d'))
    ax.xaxis.set_minor_formatter(mtick.FormatStrFormatter('%d'))

    """
    # The Second Panel
    """
    ax = fig2.add_subplot(1,3,2)
    counts1= np.histogram2d(data['snr'], data['absE'], 32, weights=data['weight'], range=[[10,80],[0.,2]], normed=True)[0]
    c = ax.contour(counts1.T,extent=[10,80,0.,2],cmap='viridis', levels=[0.0065, 0.012, 0.025, 0.035])
    counts2= np.histogram2d(sims['snr'], sims['absE'], 32, weights=sims['weight'], range=[[10,80],[0.,2]], normed=True)[0]
    c = ax.contour(counts2.T,extent=[10,80,0.,2], cmap='viridis', linestyles='dashed', levels=[0.0065,0.012,0.025,0.035])
    ax.hist2d(data['snr'], data['absE'], 128, weights=data['weight'], range=[[10,80],[0.,2]], cmap='gray_r', normed=True,alpha=0.8,norm=LogNorm(vmin=0.0002,vmax=0.2));
    ax.set_xlabel('CModel SNR')
    ax.set_ylabel(r'reGaus $|e|$')
    #ax.set_xscale('log')
    ax.xaxis.set_major_formatter(mtick.FormatStrFormatter('%d'))
    ax.xaxis.set_minor_formatter(mtick.FormatStrFormatter('%d'))

    """
    #  Adding legend
    """
    data_line = mlines.Line2D([], [], color='black', ls='-', label='Data')
    sims_line = mlines.Line2D([], [], color='black', ls='--',label='Simulation')
    ax.legend(handles=[data_line,sims_line])

    """
    # The third panel
    """
    ax = fig2.add_subplot(1,3,3)
    counts1,xbins1,ybins1,image = ax.hist2d(data['snr'], data['mag'], 32, weights=data['weight'], range=[[10,80],[22,24.5]], cmap='gray_r', normed=True,vmax=0.1);
    c = ax.contour(counts1.T,extent=[10,80,22,24.5],cmap='viridis', levels=[0.01, 0.025, 0.05, 0.075])
    counts2,xbins2,ybins2,image = ax.hist2d(sims['snr'], sims['mag'], 32, weights=sims['weight'], range=[[10,80],[22,24.5]], cmap='gray_r', normed=True,vmax=0.1);
    c = ax.contour(counts2.T,extent=[10,80,22,24.5], cmap='viridis', linestyles='dashed', levels=[0.01,0.025,0.05,0.075])
    ax.hist2d(data['snr'], data['mag'], 128, weights=data['weight'], range=[[10,80],[22,24.5]], cmap='gray_r', normed=True,vmax=0.1,alpha=0.8,norm=LogNorm(vmin=0.005,vmax=0.1));
    ax.set_xlabel('CModel SNR')
    ax.set_ylabel('CModel mag')
    ax.set_ylim(24.5,22)
    #ax.set_xscale('log')
    ax.xaxis.set_major_formatter(mtick.FormatStrFormatter('%d'))
    ax.xaxis.set_minor_formatter(mtick.FormatStrFormatter('%d'))

    """
    # Save figure
    """
    plt.tight_layout()
    fig2.savefig(os.path.join(pltDir,'galDis2D_final_logbltm037_npass2_mag24.5_snrA10_snrlbT2.png'))
    return

def nullFun():
    return

if __name__=='__main__':
    main()
