#!/usr/bin/env python
"""
Script to loop over simulated subfields and check the noise variance
estimated by the HSC pipeline, as compared to the true one.
"""
import os
import itertools
import astropy.io.fits as pyfits
import numpy as np
import astropy.io.ascii as ascii
from multiprocessing import Pool
from astropy.table import Table
import matplotlib
matplotlib.use('agg')
matplotlib.rcParams['ps.useafm'] = True
matplotlib.rcParams['pdf.use14corefonts'] = True
matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt

def histoGram(ifile):
    srcFname=   './outcomeFPFS/meas%04d.fits' %(ifile)
    if not os.path.exists(srcFname):
        return
    data    =   pyfits.getdata(srcFname)
    mask    =   np.ones(len(data)).astype(bool)
    mask    =   mask&(~np.isnan(data['base_SdssShape_psf_xx']))
    mask    =   mask&(~np.isnan(data['base_SdssShape_psf_yy']))
    mask    =   mask&(~np.isnan(data['base_SdssShape_psf_xy']))
    ss      =   data[mask][0]
    rpsf_trace  =   np.sqrt((ss['base_SdssShape_psf_xx']+ss['base_SdssShape_psf_yy'])/2.)
    rpsf_det    =   (ss['base_SdssShape_psf_xx']*ss['base_SdssShape_psf_yy']-ss['base_SdssShape_psf_xy']**2.)**(0.25)
    preFname=   './psfPre/psf%04d.fits'%(ifile) 
    pre     =   pyfits.getheader(preFname)
    var_target= pre['var']
    field   =   pre['field']
    expFname=   './expSim/exposure%04d.fits'%(ifile)
    img     =   pyfits.getdata(expFname,1)
    msk     =   pyfits.getdata(expFname,2)
    img     =   img[msk==0]
    img     =   img[~np.isnan(img)]
    var_meas=   np.var(img)
    return (ifile,field,var_target,var_meas,rpsf_trace,rpsf_det)

def plotHist(n_bins,ratio,dyn_range,out_pref):
    # Histogram of (estimated/input) noise variance.
    fig = plt.figure()
    ax = fig.add_subplot(111)
    if dyn_range > 5:
        my_bins = np.logspace(np.log10(np.min(ratio)), np.log10(np.max(ratio)), n_bins+1)
        h, ed = np.histogram(ratio, bins=my_bins)
    else:
        h, ed = np.histogram(ratio, bins=n_bins)
    x = 0.5*(ed[1:]+ed[:-1])
    ax.plot(x, h)
    if dyn_range > 5:
        ax.set_xscale('log')
    plt.xlabel('Estimated/input variance')
    plt.savefig('%s-hist.png'%out_pref)
    return


def plotScatt(dyn_range,ratio,true_noise_var_list,psf_fwhm,out_pref):
    # Scatter plot for that quantity vs. seeing.
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.scatter(psf_fwhm, ratio)
    plt.xlabel('PSF FWHM [arcsec]')
    plt.ylabel('Estimated/input variance')
    if dyn_range > 5:
        ax.set_yscale('log')
    plt.savefig('%s-scatter1.png'%out_pref)

    # Scatter plot for that quantity vs. input variance.
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.scatter(np.array(true_noise_var_list), ratio)
    if dyn_range > 5:
        ax.set_yscale('log')
    plt.xlabel('Input variance')
    plt.ylabel('Estimated/input variance')
    plt.savefig('%s-scatter2.png'%out_pref)
    return

if __name__ == '__main__':
    outFname=   'noiVar_rPsf_s18a.csv'
    if not os.path.exists(outFname):
        ndat    =   1600
        files   =   range(ndat)
        p       =   Pool(40)

        datarows=   p.map(histoGram,files)
        datarows=   [x for x in datarows if x is not None]
        names   =   ('num','field','var_target','var_meas','rpsf_trace','rpsf_det')
        t   =   Table(rows=datarows, names=names)
        t.write(outFname)
    else:
        t   =   Table.read(outFname)

    ratio   =   np.array(t['var_meas']/t['var_target'])
    psf_fwhm=   np.array(t['rpsf_det']*0.168*2.355)
    var_target= np.array(t['var_target'])
    dyn_range=  np.max(ratio)/np.min(ratio)
    out_pref=  'checkNoiseVar'
    n_bins  =   20
    plotHist(n_bins,ratio,dyn_range,out_pref)
    plotScatt(dyn_range,ratio,var_target,psf_fwhm,out_pref)
