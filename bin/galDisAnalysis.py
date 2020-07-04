#!/usr/bin/env python
import os
import sys
import glob
import fitsio
import optparse
import sim_analysis_utilities as utilities
import numpy as np
from astropy.table import Table, vstack

import matplotlib
matplotlib.use('agg')
matplotlib.rcParams['ps.useafm'] = True
matplotlib.rcParams['pdf.use14corefonts'] = True
matplotlib.rcParams['text.usetex'] = True

import matplotlib.pyplot as plt
import lsst.daf.base as dafBase
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase

from lsst.pipe.base import TaskRunner
from lsst.ctrl.pool.parallel import BatchPoolTask
from lsst.ctrl.pool.pool import Pool, abortOnError

anaDir  =   'anaCat'
obsDir  =   '../../database/s19-Hironao-2/fields/'
pltDir  =   'plot/galDis1D-2/'
## galDis1D is without magA10 cut
## galDis1D-2 is with magA10 cut

if not os.path.isdir(anaDir):
    raise RuntimeError("Simulation catalog directory %s does not exist "
         "or is not a directory" % (anaDir))
if not os.path.isdir(obsDir):
    raise RuntimeError("Data catalog directory %s does not exist "
         "or is not a directory" % (anaDir))
if not os.path.exists(pltDir):
    os.mkdir(pltDir)


nbin2      =   60
resolution_bin=np.linspace(0.1,1.,nbin2)
resolution_bin_0  =   resolution_bin[:-1]
resolution_bin_1  =   resolution_bin[1:]
resGrid    =   (resolution_bin_0+resolution_bin_1)/2.
snr_bin    =   np.linspace(0.,80.,nbin2)
snr_bin_0  =   snr_bin[:-1]
snr_bin_1  =   snr_bin[1:]
snrGrid    =   (snr_bin_0+snr_bin_1)/2.
mag_bin    =   np.linspace(22.,25.2,nbin2)
mag_bin_0  =   mag_bin[:-1]
mag_bin_1  =   mag_bin[1:]
magGrid    =   (mag_bin_0+mag_bin_1)/2.
absE_bin   =   np.linspace(0.,2.,nbin2)
absE_bin_0 =   absE_bin[:-1]
absE_bin_1 =   absE_bin[1:]
absEGrid   =   (absE_bin_0+absE_bin_1)/2.
datanames  =   ['gridRes','histRes','gridSnr','histSnr','gridMag','histMag','gridAbsE','histAbsE','num']

def main():
    #Prepare the pool
    for fname in utilities.field_keys.values():
        process(fname)
    finalPlot()
    return

def process(fieldname):
    print('Processing field: %s' %(fieldname))
    anaFieldDir =   os.path.join(anaDir,'%s' %(fieldname))

    print('Achieving the galax distributions from observations')
    obsFname =  os.path.join(obsDir,'%s.fits' %(fieldname))
    obsOFname=  os.path.join(anaFieldDir,'galDisObs.csv')
    if not os.path.exists(obsOFname):
        t_dis_obs   =   getObsGalDis(obsFname)
        t_dis_obs.write(obsOFname,overwrite=True)
    else:
        t_dis_obs   =   Table.read(obsOFname)

    print('Achieving the galaxy distributions from simulations')
    simFname    =   os.path.join(anaFieldDir,'%s.fits' %(fieldname))
    simOFname   =   os.path.join(anaFieldDir,'galDisSim.csv')
    if not os.path.exists(simOFname):
        t_dis_sim   =   getSimGalDis(simFname)
        t_dis_sim.write(simOFname,overwrite=True)
    else:
        t_dis_sim   =   Table.read(simOFname)

    plotGalDis(t_dis_obs,t_dis_sim,anaFieldDir)
    return

def getObsGalDis(obsFname):
    data    =   fitsio.read(obsFname)
    magA10  =   utilities.get_imag_A10(data)
    mask    =   np.ones(len(data)).astype(bool)
    mask    &=  (magA10<25.5)
    """
    mask    =   mask&(utilities.get_snr(data) < 80.)
    data    =   data[mask]
    """
    absE    =   utilities.get_abs_ellip(data)
    Res     =   utilities.get_res(data)
    mag     =   utilities.get_imag(data)
    snr     =   utilities.get_snr(data)
    ndata   =   len(data)*np.ones(nbin2-1)
    histRes =   np.histogram(Res,bins=resolution_bin,density=True)[0]
    histSnr =   np.histogram(snr,bins=snr_bin,density=True)[0]
    histMag =   np.histogram(mag,bins=mag_bin,density=True)[0]
    histAbsE=   np.histogram(absE,bins=absE_bin,density=True)[0]
    obscols = [resGrid,histRes,snrGrid,histSnr,magGrid,histMag,absEGrid,histAbsE,ndata]
    t_dis_obs=   Table(data=obscols,names=datanames)
    return t_dis_obs

def getSimGalDis(simFname):
    #Get the noise variance and psf information
    data    =   fitsio.read(simFname)
    data    =   utilities.update_wl_flags_sim(data)
    magA10  =   utilities.get_imag_A10(data)
    mask    =   data['weak_lensing_flag']
    mask    &=  (magA10<25.5)
    #mask    =   mask&(utilities.get_snr(data) < 80.)
    data    =   data[mask]
    ndata   =   len(data)*np.ones(nbin2-1)
    absE    =   utilities.get_abs_ellip(data)
    Res     =   utilities.get_res(data)
    mag     =   utilities.get_imag(data)
    snr     =   utilities.get_snr(data)
    weight  =   data['reweight_s16a']
    histRes =   np.histogram(Res,bins=resolution_bin,weights=weight,density=True)[0]
    histSnr =   np.histogram(snr,bins=snr_bin,weights=weight,density=True)[0]
    histMag =   np.histogram(mag,bins=mag_bin,weights=weight,density=True)[0]
    histAbsE=   np.histogram(absE,bins=absE_bin,weights=weight,density=True)[0]
    obscols = [resGrid,histRes,snrGrid,histSnr,magGrid,histMag,absEGrid,histAbsE,ndata]
    t_dis_sim=   Table(data=obscols,names=datanames)
    return t_dis_sim

def plotGalDis(t_dis_obs,t_dis_rw,anaFieldDir):
    fig =   plt.figure()
    ax  =   fig.add_subplot(111)
    ax.plot(t_dis_obs['gridRes'],t_dis_obs['histRes'],'-',label='obs',c='b')
    ax.plot(t_dis_rw['gridRes'],t_dis_rw['histRes'],'--',label='reweight',c='r')
    ax.legend()
    plt.xlabel('resolution')
    plt.ylabel('prob')
    plt.savefig(os.path.join(anaFieldDir,'resHist.eps'))
    plt.close()

    fig =   plt.figure()
    ax  =   fig.add_subplot(111)
    ax.plot(t_dis_obs['gridSnr'],t_dis_obs['histSnr'],'-',label='obs',c='b')
    ax.plot(t_dis_rw['gridSnr'],t_dis_rw['histSnr'],'--',label='reweight',c='r')
    ax.legend()
    plt.xlabel('CModel SNR')
    plt.ylabel('prob')
    plt.savefig(os.path.join(anaFieldDir,'snrHist.eps'))
    plt.close()

    fig =   plt.figure()
    ax  =   fig.add_subplot(111)
    ax.plot(t_dis_obs['gridMag'],t_dis_obs['histMag'],'-',label='obs',c='b')
    ax.plot(t_dis_rw['gridMag'],t_dis_rw['histMag'],'--',label='reweight',c='r')
    ax.legend()
    plt.xlabel('CModel mag')
    plt.ylabel('prob')
    plt.savefig(os.path.join(anaFieldDir,'magHist.eps'))
    plt.close()

    fig =   plt.figure()
    ax  =   fig.add_subplot(111)
    ax.plot(t_dis_obs['gridAbsE'],t_dis_obs['histAbsE'],'-',label='obs',c='b')
    ax.plot(t_dis_rw['gridAbsE'],t_dis_rw['histAbsE'],'--',label='reweight',c='r')

    ax.legend()
    plt.xlabel('reGaus |e|')
    plt.ylabel('prob')
    plt.savefig(os.path.join(anaFieldDir,'absEHist.eps'))
    plt.close()
    return

def finalPlot():
    obsDis_snr  =   np.zeros(nbin2-1)
    obsDis_res  =   np.zeros(nbin2-1)
    obsDis_mag  =   np.zeros(nbin2-1)
    obsDis_abs  =   np.zeros(nbin2-1)
    obsNum      =   np.zeros(nbin2-1)

    simRDis_snr =   np.zeros(nbin2-1)
    simRDis_res =   np.zeros(nbin2-1)
    simRDis_mag =   np.zeros(nbin2-1)
    simRDis_abs =   np.zeros(nbin2-1)
    simRNum     =   np.zeros(nbin2-1)

    figN=   plt.figure()
    figP=   plt.figure()
    figM=   plt.figure()
    figA=   plt.figure()

    for ii,fieldKey in enumerate(utilities.field_keys):
        fieldname   =   utilities.field_keys[fieldKey]
        anaFieldDir =   os.path.join(anaDir,'%s' %(fieldname))

        disnameObs  =   os.path.join(anaFieldDir,'galDisObs.csv')
        t_dis_obs   =   Table.read(disnameObs)
        obsDis_snr  +=  t_dis_obs['histSnr']*t_dis_obs['num']
        obsDis_res  +=  t_dis_obs['histRes']*t_dis_obs['num']
        obsDis_mag  +=  t_dis_obs['histMag']*t_dis_obs['num']
        obsDis_abs  +=  t_dis_obs['histAbsE']*t_dis_obs['num']
        obsNum      +=  t_dis_obs['num']

        disnameRw   =   os.path.join(anaFieldDir,'galDisSim.csv')
        t_dis_rw    =   Table.read(disnameRw)
        simRDis_snr +=  t_dis_rw['histSnr']*t_dis_rw['num']
        simRDis_res +=  t_dis_rw['histRes']*t_dis_rw['num']
        simRDis_mag +=  t_dis_rw['histMag']*t_dis_rw['num']
        simRDis_abs +=  t_dis_rw['histAbsE']*t_dis_rw['num']
        simRNum     +=  t_dis_rw['num']

        ax      =   figN.add_subplot(3,3,ii+1)
        gg1     =   t_dis_obs['gridSnr']
        hh1     =   t_dis_obs['histSnr']
        gg2     =   t_dis_rw['gridSnr']
        hh2     =   t_dis_rw['histSnr']
        dg      =   gg1[1]-gg1[0]
        ax.plot(gg1,hh1,'-',label='obs',c='b')
        ax.plot(gg2,hh2,'--',label='reweight',c='r')
        ax.axvline(x=sum(gg1*hh1)*dg,ls='-',c='blue')
        ax.axvline(x=sum(gg2*hh2)*dg,ls='--',c='red')
        ax.set_xlabel('CModel SNR')
        ax.set_title(fieldname)

        ax      =   figP.add_subplot(3,3,ii+1)
        ax.set_xlabel('reGaus resolution')
        gg1     =   t_dis_obs['gridRes']
        hh1     =   t_dis_obs['histRes']
        gg2     =   t_dis_rw['gridRes']
        hh2     =   t_dis_rw['histRes']
        dg      =   gg1[1]-gg1[0]
        ax.plot(gg1,hh1,'-',label='obs',c='b')
        ax.plot(gg2,hh2,'--',label='reweight',c='r')
        ax.axvline(x=sum(gg1*hh1)*dg,ls='-',c='blue')
        ax.axvline(x=sum(gg2*hh2)*dg,ls='--',c='red')
        ax.set_title(fieldname)

        ax      =   figM.add_subplot(3,3,ii+1)
        ax.set_xlabel('CModel mag')
        gg1     =   t_dis_obs['gridMag']
        hh1     =   t_dis_obs['histMag']
        gg2     =   t_dis_rw['gridMag']
        hh2     =   t_dis_rw['histMag']
        dg      =   gg1[1]-gg1[0]
        ax.plot(gg1,hh1,'-',label='obs',c='b')
        ax.plot(gg2,hh2,'--',label='reweight',c='r')
        ax.axvline(x=sum(gg1*hh1)*dg,ls='-',c='blue')
        ax.axvline(x=sum(gg2*hh2)*dg,ls='--',c='red')
        ax.set_title(fieldname)

        ax      =   figA.add_subplot(3,3,ii+1)
        ax.set_xlabel(r'reGaus $|e|$')
        gg1     =   t_dis_obs['gridAbsE']
        hh1     =   t_dis_obs['histAbsE']
        gg2     =   t_dis_rw['gridAbsE']
        hh2     =   t_dis_rw['histAbsE']
        dg      =   gg1[1]-gg1[0]
        ax.plot(gg1,hh1,'-',label='obs',c='b')
        ax.plot(gg2,hh2,'--',label='reweight',c='r')
        ax.axvline(x=sum(gg1*hh1)*dg,ls='-',c='blue')
        ax.axvline(x=sum(gg2*hh2)*dg,ls='--',c='red')
        ax.set_title(fieldname)

    obsDis_snr /=   obsNum
    obsDis_res /=   obsNum
    obsDis_mag /=   obsNum
    obsDis_abs /=   obsNum

    simRDis_snr /=  simRNum
    simRDis_res /=  simRNum
    simRDis_mag /=  simRNum
    simRDis_abs /=  simRNum

    ax      =   figN.add_subplot(3,3,7)
    gg      =   snrGrid
    hh1     =   obsDis_snr
    hh2     =   simRDis_snr
    dg      =   gg[1]-gg[0]
    ax.plot(gg,hh1,'-',label='obs',c='b')
    ax.plot(gg,hh2,'--',label='reweight',c='r')
    ax.axvline(x=sum(gg*hh1)*dg,ls='-',c='blue')
    ax.axvline(x=sum(gg*hh2)*dg,ls='--',c='red')
    ax.set_xlabel('CModel SNR')
    ax.set_title('All')
    figN.tight_layout()
    ax.legend(loc='center', bbox_to_anchor=(1.5, 0.5),fontsize=15)
    figN.savefig(os.path.join(pltDir,'snrDis.eps'))

    ax      =   figP.add_subplot(3,3,7)
    gg      =   resGrid
    hh1     =   obsDis_res
    hh2     =   simRDis_res
    dg      =   gg[1]-gg[0]
    ax.plot(gg,hh1,'-',label='obs',c='b')
    ax.plot(gg,hh2,'--',label='reweight',c='r')
    ax.axvline(x=sum(gg*hh1)*dg,ls='-',c='blue')
    ax.axvline(x=sum(gg*hh2)*dg,ls='--',c='red')
    ax.set_xlabel('reGaus resolution')
    ax.set_title('All')
    figP.tight_layout()
    ax.legend(loc='center', bbox_to_anchor=(1.5, 0.5),fontsize=15)
    figP.savefig(os.path.join(pltDir,'resolutionDis.eps'))

    ax      =   figM.add_subplot(3,3,7)
    gg      =   magGrid
    hh1     =   obsDis_mag
    hh2     =   simRDis_mag
    dg      =   gg[1]-gg[0]
    ax.plot(gg,hh1,'-',label='obs',c='b')
    ax.plot(gg,hh2,'--',label='reweight',c='r')
    ax.axvline(x=sum(gg*hh1)*dg,ls='-',c='blue')
    ax.axvline(x=sum(gg*hh2)*dg,ls='--',c='red')
    ax.set_xlabel('CModel mag')
    ax.set_title('All')
    figM.tight_layout()
    ax.legend(loc='center', bbox_to_anchor=(1.5, 0.5),fontsize=15)
    figM.savefig(os.path.join(pltDir,'magDis.eps'))

    ax      =   figA.add_subplot(3,3,7)
    gg      =   absEGrid
    hh1     =   obsDis_abs
    hh2     =   simRDis_abs
    dg      =   gg[1]-gg[0]
    ax.plot(gg,hh1,'-',label='obs',c='b')
    ax.plot(gg,hh2,'--',label='reweight',c='r')
    ax.axvline(x=sum(gg*hh1)*dg,ls='-',c='blue')
    ax.axvline(x=sum(gg*hh2)*dg,ls='--',c='red')
    ax.set_xlabel(r'reGaus $|e|$')
    ax.set_title('All')
    figA.tight_layout()
    ax.legend(loc='center', bbox_to_anchor=(1.5, 0.5),fontsize=15)
    figA.savefig(os.path.join(pltDir,'absEDis.eps'))
    return

if __name__=='__main__':
    main()
