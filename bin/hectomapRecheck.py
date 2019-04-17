#!/usr/bin/env python
import os
import numpy as np
import matplotlib
matplotlib.use('agg')
matplotlib.rcParams['ps.useafm'] = True
matplotlib.rcParams['pdf.use14corefonts'] = True
matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt
import astropy.io.fits as pyfits
from astropy.table import Table

def main():
    hiroFname   =   '/work/xiangchong.li/work/S16AStandard/S16AStandardCalibrated/field/HECTOMAP_calibrated.fits'
    a           =   pyfits.getdata(hiroFname)
    b           =   a['ishape_sdss_psf_ixx']*a['ishape_sdss_psf_iyy']-a['ishape_sdss_psf_ixy']**2.
    b           =   b**(0.25)
    b           =   b*2.35
    #s16a-4.5
    nbin1       =   20
    fwhm_bin    =   np.linspace(0.3,0.9,nbin1)
    fwhm_bin_0  =   fwhm_bin[:-1]
    fwhm_bin_1  =   fwhm_bin[1:]
    gridPsf     =   (fwhm_bin_0+fwhm_bin_1)/2.
    histPsf     =   np.histogram(b,bins=fwhm_bin,density=True)[0]
    anaFieldDir =   '/work/xiangchong.li/work/S18ASimHSC/s16/anaCat/HECTOMAP/'
    fnameRw     =   os.path.join(anaFieldDir,'varPSF_dis_rw.csv')
    fnameOb     =   os.path.join(anaFieldDir,'varPSF_dis_obs.csv')
    t_dis_rw    =   Table.read(fnameRw)
    t_dis_ob    =   Table.read(fnameOb)
    
    fig =   plt.figure()
    ax  =   fig.add_subplot(111)
    ax.set_xlabel('PSF FWHM [arcsec]')
    ax.set_ylabel('prob')
    ax.plot(gridPsf,histPsf,'--',label='s16a-calibrated',c='b')
    ax.plot(t_dis_rw['gridPsf'],t_dis_rw['histPsf'],'--',label='sim',c='r')
    ax.plot(t_dis_ob['gridPsf'],t_dis_ob['histPsf'],'--',label='s16a-6.7',c='black')
    ax.legend()
    
    plt.savefig('problem_hectomap_fwhm.png')
    
    return
    
if __name__ == '__main__':
    main()
