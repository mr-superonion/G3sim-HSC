#!/usr/bin/env python
import os
import numpy as np
from astropy.table import Table, vstack
import pandas as pd
import sim_analysis_utilities as utilities

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
anaDir  =   'anaCat_newS19Mask'
# Checks that the paths provided are valid
if not os.path.isdir(anaDir):
    raise RuntimeError("Simulation catalog directory %s does not exist "
        "or is not a directory" % (anaDir))
simAname=   os.path.join(anaDir,'allSim.h5')

obsDir =   '../../database/'
if not os.path.isdir(obsDir):
    raise RuntimeError("Data catalog directory %s does not exist "
         "or is not a directory" % (obsDir))
obsAname=   os.path.join(obsDir,'smallcat-wl-s19-2-2.fits')

pltDir2 =   'plot/galDis2D/'
if not os.path.exists(pltDir2):
    os.system('mkdir -p %s' %pltDir2)

datanames   =   ['res','snr','mag','absE','logb','weight']

def main():
    allSim,allObs=getAllDat()
    for maglim in np.arange(25.5,25.4,-0.1):
        finalPlot(allSim,allObs,'magA10',maglim,'wl')
    return

def getAllDat():
    allSim  =   None
    print('Achieving galaxy data from observations')
    allObs  =   Table.read(obsAname)

    print('Achieving galaxy data from simulations')
    allSim  =   pd.read_hdf(simAname,mode='r').to_records(index=True)

    return allSim,allObs

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
    if oname!='mag':
        mask    &=  ((data['mag']-data['a_i'])<24.5)
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
        maskS   &=  (sims['res']>0.3)
        maskS   &=  (sims['snr']>10.)
        maskS   &=  (sims['absE']<2.)
        if oname!='mag':
            maskS   &=  ((sims['mag']-sims['a_i'])<24.5)
        maskS   &=  (sims[oname]<olim)
        sims    =   sims[maskS]
        if 'weight' in sims.dtype.names:
            weight_sims =   sims['weight']
        else:
            weight_sims =   np.ones(len(data))

    fig2,axes = plt.subplots(1,3,figsize=(15, 5),constrained_layout=True)
    fig2.suptitle('%s: %s %.1f' %(version,oname,olim,), fontsize=18)

    """
    # The First Panel
    """
    ax = axes[0]
    counts1= np.histogram2d(data['snr'], data['res'], 32, weights=weight_data, range=[[10,80],[0.3,1]], normed=True)[0]
    ax.contour(counts1.T,extent=[10,80,0.3,1],cmap='viridis', levels=[0.01, 0.025, 0.05, 0.075])
    ax.hist2d(data['snr'], data['res'], 128, weights=weight_data, range=[[10,80],[0.3,1]], cmap='gray_r', normed=True,vmax=0.08,alpha=0.8);
    if sims is not None:
        counts2= np.histogram2d(sims['snr'], sims['res'], 32, weights=weight_sims, range=[[10,80],[0.3,1]], normed=True)[0]
        ax.contour(counts2.T,extent=[10,80,0.3,1],cmap='viridis',linestyles='dashed', levels=[0.01, 0.025, 0.05, 0.075])
    ax.xaxis.set_major_formatter(mtick.FormatStrFormatter('%d'))
    ax.set_xlabel('CModel SNR')
    ax.set_ylabel('reGaus resolution')
    #ax.set_xscale('log')

    """
    # The Second Panel
    """
    ax = axes[1]
    counts1= np.histogram2d(data['snr'], data['absE'], 32, weights=weight_data, range=[[10,80],[0.,2]], normed=True)[0]
    ax.contour(counts1.T,extent=[10,80,0.,2],cmap='viridis', levels=[0.0065, 0.012, 0.025, 0.035])
    if sims is not None:
        counts2= np.histogram2d(sims['snr'], sims['absE'], 32, weights=weight_sims, range=[[10,80],[0.,2]], normed=True)[0]
        ax.contour(counts2.T,extent=[10,80,0.,2], cmap='viridis', linestyles='dashed', levels=[0.0065,0.012,0.025,0.035])
    ax.hist2d(data['snr'], data['absE'], 128, weights=weight_data, range=[[10,80],[0.,2]], cmap='gray_r', normed=True,alpha=0.8,norm=LogNorm(vmin=0.0002,vmax=0.2));
    ax.xaxis.set_major_formatter(mtick.FormatStrFormatter('%d'))
    ax.set_xlabel('CModel SNR')
    ax.set_ylabel(r'reGaus $|e|$')
    #ax.set_xscale('log')

    """
    #  Adding legend
    """
    data_line = mlines.Line2D([], [], color='black', ls='-', label='Data')
    sims_line = mlines.Line2D([], [], color='black', ls='--',label='Simulation')
    ax.legend(handles=[data_line,sims_line])

    """
    # The third panel
    """
    ax = axes[2]
    counts1,xbins1,ybins1,image = ax.hist2d(data['snr'], data['mag'], 32, weights=weight_data, range=[[10,80],[22,24.5]], cmap='gray_r', normed=True,vmax=0.1);
    ax.contour(counts1.T,extent=[10,80,22,24.5],cmap='viridis', levels=[0.01, 0.025, 0.05, 0.075])
    if sims is not None:
        counts2,xbins2,ybins2,image = ax.hist2d(sims['snr'], sims['mag'], 32, weights=weight_sims, range=[[10,80],[22,24.5]], cmap='gray_r', normed=True,vmax=0.1);
        ax.contour(counts2.T,extent=[10,80,22,24.5], cmap='viridis', linestyles='dashed', levels=[0.01,0.025,0.05,0.075])
    ax.hist2d(data['snr'], data['mag'], 128, weights=weight_data, range=[[10,80],[22,24.5]], cmap='gray_r', normed=True,vmax=0.1,alpha=0.8,norm=LogNorm(vmin=0.005,vmax=0.1));
    ax.xaxis.set_major_formatter(mtick.FormatStrFormatter('%d'))
    ax.set_xlabel('CModel SNR')
    ax.set_ylabel('CModel mag')
    ax.set_ylim(24.5,22)
    #ax.set_xscale('log')

    """
    # Save figure
    """
    fig2.savefig(os.path.join(pltDir2,'galDis2D-final-%s-%.1f-%s.eps' %(oname,olim,pend)))
    return

if __name__=='__main__':
    main()
