#!/usr/bin/env python
import os
import fitsio
import treecorr
import numpy as np
from astropy.table import Table

import matplotlib as mpl
mpl.use('TkAgg')
mpl.rcParams['ps.useafm'] = True
mpl.rcParams['pdf.use14corefonts'] = True
mpl.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt

"""
# Prepare the magnitude bins
"""
nMBin       =   4
magEdge     =   np.linspace(2,18,nMBin+1)
magBin      =   (magEdge[:-1]+magEdge[1:])//2

def main():
    """
    # Prepare for Plot
    """
    pltDir      =   'plot/galDis/'
    fig         =   plt.figure()
    ax          =   fig.add_subplot(111)
    cor0        =   treecorr.NNCorrelation(nbins=10,min_sep=0.01,max_sep=1.5,sep_units='deg')
    x           =   cor0.rnom
    solidAng    =   ((x*np.exp(cor0.bin_size/2))**2.-(x/np.exp(cor0.bin_size/2)))*np.pi
    print(solidAng)

    """
    # Load the observed data
    """
    anaDir      =   'anaCat'
    data        =   fitsio.read(os.path.join(anaDir,'allObs4.fits'))
    """
    # Select the extended faint galaxies
    """
    res         =   data['res']
    snr         =   data['snr']
    mask        =   (res>0.8)&(res<1.)&(snr>10.)&(snr<15.)
    data        =   data[mask]
    raG         =   data['ra']
    decG        =   data['dec']
    print('number of extended faint galaxies: %s' %len(data))
    catG        =   treecorr.Catalog(ra=raG,dec=decG,ra_units='deg',dec_units='deg')

    """
    # Load the star catalog
    """
    starName    =   'brightObjMask/HSC-SSP_brightStarMask_Arcturus/stars/stars_all.fits'
    stars       =   fitsio.read(starName,ext=1)
    mags        =   stars['G_Gaia']
    print('number of bright stars: %d' %len(stars))

    for i in range(nMBin):
        mag     =   magBin[i] 
        maskU   =   (mags>magEdge[i])& (mags<magEdge[i+1])
        starU   =   stars[maskU]
        numS    =   len(starU)
        raS     =   starU['ra']
        decS    =   starU['dec']
        cor     =   treecorr.NNCorrelation(nbins=10,min_sep=0.01,max_sep=1.5,sep_units='deg')
        catS    =   treecorr.Catalog(ra=raS,dec=decS,ra_units='deg',dec_units='deg')
        cor.process(catS,catG)
        y       =   cor.npairs/numS/solidAng
        ax.plot(x,y,label='%.1f' %mag)

    plt.xlabel(r'$\theta$ (deg)')
    plt.ylabel(r'num/${\rm{deg}}^2$')
    plt.tight_layout()
    plt.legend()
    plt.savefig(os.path.join(pltDir,'WL_requirement.png'))
    return

if __name__=='__main__':
    main()
