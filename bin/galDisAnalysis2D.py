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
matplotlib.rcParams['font.size'] = 14

import matplotlib.ticker as ticker
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import matplotlib.lines as mlines
from matplotlib.colors import LogNorm

version =   's19'
matver  =   's16'
anaDir  =   'anaCat'
# Checks that the paths provided are valid
if not os.path.isdir(anaDir):
    raise RuntimeError("Simulation catalog directory %s does not exist "
         "or is not a directory" % (anaDir))
simAname=   os.path.join(anaDir,'allSim.fits')

obsDir  =   '../database/%s/fields/' %version
obsDir0 =   '../database/'
if not os.path.isdir(obsDir):
    raise RuntimeError("Data catalog directory %s does not exist "
         "or is not a directory" % (obsDir))
obsAname=   os.path.join(obsDir0,'allObs-%s.fits' %version)

pltDir  =   'plot/galDis-%s/' %version
if not os.path.exists(pltDir):
    os.system('mkdir -p %s' %pltDir)

pltDir2 =   'plot/galDis-%s-final/' %version
if not os.path.exists(pltDir2):
    os.system('mkdir -p %s' %pltDir2)

datanames   =   ['res','snr','mag','absE','logb','weight']

def main():
    allSim,allObs=getAllDat2()
    for maglim in np.arange(27,24.0,-0.1):
        finalPlot(allSim,allObs,'magA10',maglim,'scat')
    return

def getAllDat2():
    allSim=None
    allObs=pyfits.getdata('../database/compare/scat-%s.fits' %version)
    return allSim,allObs

def getAllDat():
    allSim  =   None
    print('Achieving galaxy data from observations')
    if not os.path.exists(obsAname):
        allObs  =   list(map(getObsGalData,utilities.field_keys.values()))
        allObs  =   [x for x in allObs if x is not None]
        allObs  =   vstack(allObs)
        allObs.write(obsAname)
    else:
        print('Already have prepared catalog')
        allObs=Table.read(obsAname)

    print('Achieving galaxy data from simulations')
    if not os.path.exists(simAname):
        allSim  =   list(map(getSimGalData,utilities.field_keys.values()))
        allSim  =   [x for x in allSim if x is not None]
        allSim  =   vstack(allSim)
        allSim.write(simAname)
    else:
        print('Already have prepared catalog')
        allSim=Table.read(simAname)
    return allSim,allObs

def getObsGalData(fieldname):
    """
    #load the necessary data from catalog (after FDFC, WL cut)
    """
    obsFname=   os.path.join(obsDir,'%s.fits' %(fieldname))
    data    =   pyfits.getdata(obsFname)
    mask    =   (data['weak_lensing_flag'])
    data    =   data[mask]
    ra,dec  =   utilities.get_radec(data)

    """
    # Match to S16A
    """
    obsFname16= '../database/%s/fields/%s.fits' %(matver,fieldname)
    data16  =   pyfits.getdata(obsFname16)

    """
    # Match to S16
    """
    match   =   smatch.match(data16['ira'],data16['idec'],0.5/3600.,\
            ra,dec,maxmatch=1)

    """
    # Main Column for S16A
    """
    TPid16  =   np.zeros(len(data),dtype=int)
    TPid16[match['i2']] =   utilities.get_TPid(data16[match['i1']])
    res16   =   np.zeros(len(data))
    res16[match['i2']]  =   utilities.get_res(data16[match['i1']])
    snr16   =   np.zeros(len(data))
    snr16[match['i2']]  =   utilities.get_snr(data16[match['i1']])
    mag16   =   np.zeros(len(data))
    mag16[match['i2']]  =   utilities.get_imag(data16[match['i1']])
    logb16  =   np.zeros(len(data))
    logb16[match['i2']] =   utilities.get_logb(data16[match['i1']])
    absE16  =   np.zeros(len(data))
    absE16[match['i2']] =   utilities.get_abs_ellip(data16[match['i1']])
    absEP16 =   np.zeros(len(data))
    absEP16[match['i2']]=   utilities.get_abs_ellip_psf(data16[match['i1']])
    fwhm16  =   np.zeros(len(data))
    fwhm16[match['i2']] =   utilities.get_psf_size(data16[match['i1']])

    traceR16=   np.zeros(len(data))
    traceR16[match['i2']]=  utilities.get_gal_size(data16[match['i1']],type='trace')
    detR16  =   np.zeros(len(data))
    detR16[match['i2']] =   utilities.get_gal_size(data16[match['i1']],type='det')
    object_id16 =   np.zeros(len(data),dtype=int)
    object_id16[match['i2']]=data16['object_id'][match['i1']]
    wlFlag16=   np.zeros(len(data),dtype=bool)
    wlFlag16[match['i2']]=  data16['weak_lensing_flag'][match['i1']]
    npass16 =   np.zeros(len(data),dtype=int)
    npass16[match['i2']] =  utilities.get_npass(data16[match['i1']])
    inS16   =   (wlFlag16)&(logb16<=-0.375)&(npass16>=2)

    """
    # FDFC Column for S16A
    """
    inS16FDFC   =   utilities.get_FDFC_flag(ra,dec,'s16')

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

    """
    # Flags for S18A
    """
    objId   =   data['object_id']
    weight  =   np.ones(len(data))
    datanmA =   ['fwhm', 'magA15', 'magA20', 'traceR16', 'detR16', 'object_id16', \
            'inS16', 'inS16FDFC', 'absEP16', 'maglb', 'isIso', 'bsf', 'obj', \
            'res16', 'snr16', 'mag16', 'absE16', 'logb16', 'fwhm16', 'TPid16', \
            'ra', 'dec', 'object_id', 'fieldname', 'magA10', 'noivar', 'magpsf', \
            'traceR', 'detR', 'traceRP', 'detRP', 'absEP', 'TPid', 'bmsk18', \
            'snrA10', 'snrA15', 'snrA20', 'snrlb', 'logbR', 'logbO', 'npass' ]
    dataTab =   [res, snr, mag, absE, logb, weight, \
            fwhm, magA15, magA20, traceR16, detR16, object_id16, \
            inS16, inS16FDFC, absEP16, maglb, isIso, bsf, obj, \
            res16, snr16, mag16, absE16, logb16, fwhm16, TPid16, \
            ra, dec, objId, fieldnameA, magA10, noivar, magpsf, \
            traceR, detR, traceRP, detRP, absEP, TPid, bmsk18, \
            snrA10, snrA15, snrA20, snrlb, logbR, logbO, npass]
    del data
    return Table(data=dataTab,names=datanames+datanmA)

def getSimGalData(fieldname):
    anaFieldDir =   os.path.join(anaDir,'%s' %(fieldname))
    simFname=   os.path.join(anaFieldDir,'%s.fits' %(fieldname))
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
    magA15  =   utilities.get_imag_A15(data)
    magA20  =   utilities.get_imag_A20(data)
    maglb   =   utilities.get_imag_lb(data)
    isIso   =   utilities.get_isIso(data)
    bsf     =   utilities.get_bs_factor(data)
    fwhm    =   utilities.get_psf_size(data,'fwhm')
    fieldname=  data['field_name']

    dataTab =   [res, snr, mag, absE, logb, weight, \
                snrA10, snrA15, snrA20, snrlb, fieldname, isIso, \
                bsf, logbR, logbO, noivar, magA10, magA15, magA20, maglb, \
                fwhm]
    datanmA =   ['snrA10', 'snrA15', 'snrA20', 'snrlb', 'fieldname', 'isIso', \
                'bsf', 'logbR', 'logbO','noivar','magA10', 'magA15', 'magA20', 'maglb', \
                'fwhm']
    return Table(data=dataTab,names=datanames+datanmA)

def finalPlot(sims,data,oname,olim,pend):
    """
    # Data Mask
    """
    nbf     =   len(data)+0.
    mask    =   np.ones(len(data),dtype=bool)
    mask    &=  (data['npass']>=2)
    mask    &=  (data['logb']<=-0.375)
    mask    &=  (data['res']>0.3)
    mask    &=  (data['snr']>10.)
    mask    &=  (data['absE']<2.)
    mask    &=  (data[oname]<olim)
    data    =   data[mask]
    if 'weight' in data.dtype.names:
        weight_data =   data['weight']
    else:
        weight_data =   np.ones(len(data))

    if sims is not None:
        """
        # Simulations Mask
        """
        maskS   =   (sims['logb']<=-0.375)
        maskS   &=  (sims[oname]<olim)
        sims    =   sims[maskS]
        if 'weight' in sims.dtype.names:
            weight_sims =   sims['weight']
        else:
            weight_sims =   np.ones(len(data))

    fig2 = plt.figure(figsize=(15, 5))

    """
    # The First Panel
    """
    ax = fig2.add_subplot(1,3,1)
    counts1= np.histogram2d(data['snr'], data['res'], 32, weights=weight_data, range=[[10,80],[0.3,1]], normed=True)[0]
    c = ax.contour(counts1.T,extent=[10,80,0.3,1],cmap='viridis', levels=[0.01, 0.025, 0.05, 0.075])
    ax.hist2d(data['snr'], data['res'], 128, weights=weight_data, range=[[10,80],[0.3,1]], cmap='gray_r', normed=True,vmax=0.08,alpha=0.8);
    if sims is not None:
        counts2= np.histogram2d(sims['snr'], sims['res'], 32, weights=weight_sims, range=[[10,80],[0.3,1]], normed=True)[0]
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
    counts1= np.histogram2d(data['snr'], data['absE'], 32, weights=weight_data, range=[[10,80],[0.,2]], normed=True)[0]
    c = ax.contour(counts1.T,extent=[10,80,0.,2],cmap='viridis', levels=[0.0065, 0.012, 0.025, 0.035])
    if sims is not None:
        counts2= np.histogram2d(sims['snr'], sims['absE'], 32, weights=weight_sims, range=[[10,80],[0.,2]], normed=True)[0]
        c = ax.contour(counts2.T,extent=[10,80,0.,2], cmap='viridis', linestyles='dashed', levels=[0.0065,0.012,0.025,0.035])
    ax.hist2d(data['snr'], data['absE'], 128, weights=weight_data, range=[[10,80],[0.,2]], cmap='gray_r', normed=True,alpha=0.8,norm=LogNorm(vmin=0.0002,vmax=0.2));
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
    counts1,xbins1,ybins1,image = ax.hist2d(data['snr'], data['mag'], 32, weights=weight_data, range=[[10,80],[22,24.5]], cmap='gray_r', normed=True,vmax=0.1);
    c = ax.contour(counts1.T,extent=[10,80,22,24.5],cmap='viridis', levels=[0.01, 0.025, 0.05, 0.075])
    if sims is not None:
        counts2,xbins2,ybins2,image = ax.hist2d(sims['snr'], sims['mag'], 32, weights=weight_sims, range=[[10,80],[22,24.5]], cmap='gray_r', normed=True,vmax=0.1);
        c = ax.contour(counts2.T,extent=[10,80,22,24.5], cmap='viridis', linestyles='dashed', levels=[0.01,0.025,0.05,0.075])
    ax.hist2d(data['snr'], data['mag'], 128, weights=weight_data, range=[[10,80],[22,24.5]], cmap='gray_r', normed=True,vmax=0.1,alpha=0.8,norm=LogNorm(vmin=0.005,vmax=0.1));
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
    fig2.savefig(os.path.join(pltDir2,'galDis2D-final-%s-%.1f-%s.png' %(oname,olim,pend)))
    return

if __name__=='__main__':
    main()
