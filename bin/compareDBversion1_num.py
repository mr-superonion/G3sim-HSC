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

pltDir  =   's18-2/plot/spacial/'
if not os.path.exists(pltDir):
    os.system('mkdir -p %s' %pltDir)

def getNumY(data,magX,doMask):
    if doMask:
        mask    =   (data['npass']>=2)
        mask    &=  (data['logb']<=-0.375)
        mask    &=  (data['res']>0.3)
        mask    &=  (data['snr']>10.)
        data    =   data[mask]
    numY    =   []
    for mag in magX:
        numY.append(sum(data['mag']<mag))
    return np.array(numY)

def numPlot(version,axes):
    outfname    =   os.path.join(pltDir,'num_mag_%s.csv' %version)
    if not os.path.exists(outfname):
        magX    =   np.arange(23,26.2,0.2)
        data    =   prepare(version)
        numY1   =   getNumY(data,magX,False)
        numY2   =   getNumY(data,magX,True)
        out     =   Table(data=[magX,numY1,numY2],names=['magX','numY1','numY2'])
        out.write(outfname)
    else:
        cat=Table.read(outfname)
        magX=cat['magX']
        numY1=cat['numY1']
        numY2=cat['numY2']
    ax1 =   axes[0]
    ax1.plot(magX,numY1,label=version)
    ax2 =   axes[1]
    ax2.plot(magX,numY2,label=version)
    return

def prepare(version):
    print('Achieving galaxy data from observations')
    obsDir  =   'database/%s/fields/' %version
    if not os.path.isdir(obsDir):
        raise RuntimeError("Data catalog directory %s does not exist "
        "or is not a directory" % (obsDir))
    cmpDir =   'database/compare/'
    if not os.path.exists(cmpDir):
        os.system('mkdir -p %s' %cmpDir)
    scatname=   os.path.join(cmpDir,'scat-%s.fits' %version)
    if not os.path.exists(scatname):
        data  = []
        for fn in utilities.field_keys.values():
            data.append(getSmallCat(obsDir,fieldname))
        data  =   vstack(data)
        data.write(scatname)
    else:
        print('Already have prepared catalog')
        data=Table.read(scatname)
    return data

def getSmallCat(obsDir,fieldname):
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
    fwhm    =   utilities.get_psf_size(data,'fwhm')

    maglb   =   utilities.get_imag_lb(data)
    snrlb   =   utilities.get_snr_localBG(data)
    isIso   =   utilities.get_isIso(data)
    bsf     =   utilities.get_bs_factor(data)
    obj     =   utilities.get_cmodel_obj(data)
    objId   =   data['object_id']
    fns     =   [fieldname]*len(data)

    magA10  =   utilities.get_imag_A10(data)
    magA15  =   utilities.get_imag_A15(data)
    magA20  =   utilities.get_imag_A20(data)
    noivar  =   utilities.get_noi_var(data)
    magpsf  =   utilities.get_imag_psf(data)

    traceR  =   utilities.get_gal_size(data,type='trace')
    detR    =   utilities.get_gal_size(data,type='det')
    traceRP =   utilities.get_psf_size(data,type='trace')
    detRP   =   utilities.get_psf_size(data,type='det')
    absEP   =   utilities.get_abs_ellip_psf(data)
    TPid    =   utilities.get_TPid(data)
    bmsk18  =   utilities.get_mask_briObj(data)

    snrA10,snrA15,snrA20    =   utilities.get_snr_apertures(data)
    logb    =   utilities.get_logb(data)
    logbA,logbR,logbO   =   utilities.get_logbAll(data)
    npass   =   utilities.get_npass(data)
    del data
    gc.collect()

    datanames=  ['ra', 'dec', 'res', 'snr', 'mag', 'absE', 'fwhm', \
            'maglb', 'snrlb',  'isIso', 'bsf', 'obj', 'object_id', 'fieldname', \
            'magA10', 'magA15', 'magA20', 'noivar', 'magpsf', \
            'traceR', 'detR', 'traceRP', 'detRP', 'absEP', 'TPid', 'bmsk18', \
            'snrA10', 'snrA15', 'snrA20','logb', 'logbR', 'logbO', 'npass' ]
    dataTab =  [ra, dec, res, snr, mag, absE, fwhm ,\
            maglb, snrlb,  isIso, bsf, obj, objId, fns, \
            magA10, magA15, magA20, noivar, magpsf, \
            traceR, detR, traceRP, detRP, absEP, TPid, bmsk18, \
            snrA10, snrA15, snrA20,logb, logbR, logbO, npass ]
    return Table(data=dataTab,names=datanames)

def main():
    fig =   plt.figure(figsize=(10,5))
    ax1 =   fig.add_subplot(1,2,1)
    ax2 =   fig.add_subplot(1,2,2)
    axes=[ax1,ax2]
    numPlot('s16',axes)
    numPlot('s18',axes)
    numPlot('s19',axes)
    ax1.set_xlim(26.2,23.8)
    ax1.set_ylabel('num')
    ax1.set_xlabel('cmodel mag')
    ax1.grid()
    ax2.set_xlim(26.2,23.8)
    ax2.set_xlabel('cmodel mag')
    ax2.legend()
    ax2.grid()
    plt.tight_layout()
    fig.savefig(os.path.join(pltDir,'num_mag.png'))
    return

if __name__=='__main__':
    main()
