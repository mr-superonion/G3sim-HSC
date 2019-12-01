#!/usr/bin/env python
import os
import smatch
import fitsio
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


def makeRequirement():

    outTab      =   Table()
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
    ra          =   data['ra']
    dec         =   data['dec']
    print('number of extended faint galaxies: %s' %len(data))

    """
    # Load the star catalog
    """
    starName    =   'brightObjMask/HSC-SSP_brightStarMask_Arcturus/stars/stars_all.fits'
    stars       =   fitsio.read(starName,ext=1)
    mags        =   stars['G_Gaia']
    indMBinS    =   np.digitize(mags,magEdge) # len(stars)
    print('number of bright stars: %d' %len(stars))
    
    """
    # Match between bright Stars and galaxies
    """
    match   =   smatch.match(ra,dec,1.2,\
            stars['ra'],stars['dec'],nside=512,maxmatch=-1)


    """
    # Separation List
    """
    nRbin       =   9
    distEdge    =   np.linspace(0,1.2,nRbin+1) #in degree
    distBin     =   (distEdge[:-1]+distEdge[1:])/2.
    print('distBin: %s (deg)' %distBin)
    cosEdge     =   np.cos(distEdge*np.pi/180)
    solidAngBin =   (cosEdge[:-1]-cosEdge[1:])*2.*np.pi
    print('solidAngBin: %s' %solidAngBin)

    cosDist     =   match['cosdist']
    indRBinM    =   np.digitize(cosDist,cosEdge)    # len(match)
    indMBinM    =   indMBinS[match['i2']]           # len(match) 

    """
    # Make columns
    """
    outTab['dist']      =   distBin
    outTab['solidAng']  =   solidAngBin

    for i in range(nMBin):
        """
        # Compute perBriObjBin
        """
        perBriObjBin=   np.zeros(nRbin)
        for j in range(nRbin):
            matU        =   match[(indMBinM==(i+1))&(indRBinM==(j+1))]
            numBriObj   =   len(set(matU['i2']))
            perBriObjBin[j]=len(matU)/numBriObj
            print('number of matches %d' %len(matU))
            print('number of bright Objs %d' %numBriObj)
        outTab['%s_BadperBriObj' %i]  =   perBriObjBin
    outTab.write('./brightObjMask/WL_requirements.fits')
    return outTab

def main():
    """
    # requirement
    """
    if os.path.exists('./brightObjMask/WL_requirements.fits'):
        print('Already have the dens-radius file')
        req =   Table.read('./brightObjMask/WL_requirements.fits')
    else:
        print('Generating the dens-radius file')
        req =   makeRequirement()

    """
    # Make Plot
    """
    pltDir  =   'plot/galDis/'
    fig =   plt.figure()
    ax  =   fig.add_subplot(111)
    x   =   req['dist']
    sA  =   req['solidAng']*(180/np.pi)**2.
    for im,mag in enumerate(magBin):
        y   =   req['%s_BadperBriObj' %im]/sA
        ax.plot(x,y,label='%.1f' %mag)
    plt.xlabel(r'$\theta$ (deg)')
    plt.ylabel(r'num/${\rm{deg}}^2$')
    plt.tight_layout()
    plt.legend()
    plt.savefig(os.path.join(pltDir,'WL_requirement.png'))
    return

if __name__=='__main__':
    main()
