#!/usr/bin/env python
"""
Script to loop over simulated subfields and check the noise variance
estimated by the HSC pipeline, as compared to the true one.
"""
import os
import fitsio
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

def getData(ifile):
    srcFname=   './outcomeFPFS/src-%04d.fits' %(ifile)
    if not os.path.exists(srcFname):
        return
    data    =   pyfits.getdata(srcFname)
    mask    =   np.ones(len(data)).astype(bool)
    mask    =   mask&(~np.isnan(data['base_SdssShape_psf_xx']))
    mask    =   mask&(~np.isnan(data['base_SdssShape_psf_yy']))
    mask    =   mask&(~np.isnan(data['base_SdssShape_psf_xy']))
    mask    =   mask&(~np.isnan(data['base_Variance_value']))
    ss      =   data[mask][0]
    rpsf_trace  =   np.sqrt((ss['base_SdssShape_psf_xx']+ss['base_SdssShape_psf_yy'])/2.)
    rpsf_det    =   (ss['base_SdssShape_psf_xx']*ss['base_SdssShape_psf_yy']-ss['base_SdssShape_psf_xy']**2.)**(0.25)
    var_meas=   ss['base_Variance_value']
    preFname=   './psfPre/psf%04d.fits'%(ifile) 
    pre     =   pyfits.getheader(preFname)
    var_target= pre['var']
    field   =   pre['field']
    return (ifile,field,var_target,var_meas,rpsf_trace,rpsf_det,1.,1.)


def plotVarRatio(t):
    ratio   =   np.array(t['var_meas']/t['var_target'])
    psf_fwhm=   np.array(t['rpsf_det']*0.168*2.355)
    var_target= np.array(t['var_target'])
    dyn_range=  np.max(ratio)/np.min(ratio)
    out_pref=  'checkNoiseVar'
    n_bins  =   20
    plotHistVarRatio(n_bins,ratio,dyn_range,out_pref)
    plotScattVarRatio(dyn_range,ratio,var_target,psf_fwhm,out_pref)
    return

def plotHistVarRatio(n_bins,ratio,dyn_range,out_pref):
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


def plotScattVarRatio(dyn_range,ratio,true_noise_var_list,psf_fwhm,out_pref):
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

def prepareHistInfo(data):
    ngrid   =   20
    fwhms   =   data['rpsf_det']*0.168*2.355
    varNois =   data['var_meas']
    minVar  =   varNois.min()
    gridNoi =   np.linspace(minVar,0.01,ngrid)
    gridPsf =   np.linspace(0.35,0.85,ngrid)
    histInfo =   {'gridNoi':gridNoi,'gridPsf':gridPsf}
    np.save('noi_psf_histGram.npy',histInfo)
    return histInfo

def getHistSim(t,histInfo):
    ngrid   =   20
    fwhms   =   t['rpsf_det']*0.168*2.355
    varNois =   t['var_meas']
    gridNoi =   histInfo['gridNoi']
    gridPsf =   histInfo['gridPsf']
    histNoi =   np.histogram(varNois,bins=gridNoi,density=True)[0]
    histNoi =   histNoi/sum(histNoi)
    histPsf =   np.histogram(fwhms,bins=gridPsf,density=True)[0]
    histPsf =   histPsf/sum(histPsf)
    histInfo['histNoiSim0']=histNoi
    histInfo['histPsfSim0']=histPsf
    #do reweight for noise
    minVars =   histInfo['gridNoi'][:-1]
    maxVars =   histInfo['gridNoi'][1:]
    ratios  =   np.array(histInfo['histNoiObs']/histInfo['histNoiSim0'])
    for ts in t:
        if ts['var_meas']<=minVars[0] or ts['var_meas']>maxVars[-1]:
            ts['weightNoi']    =   0.
        else:
            ind         =   np.where((ts['var_meas']>minVars) & (ts['var_meas']<maxVars))[0][0]
            ts['weightNoi']=   ratios[ind]
    #apply the noise weight
    weight  =   t['weightNoi']
    histNoi =   np.histogram(varNois,bins=gridNoi,weights=weight,density=True)[0]
    histNoi =   histNoi/sum(histNoi)
    histPsf =   np.histogram(fwhms,bins=gridPsf,weights=weight,density=True)[0]
    histPsf =   histPsf/sum(histPsf)
    histInfo['histNoiSim1']=histNoi
    histInfo['histPsfSim1']=histPsf
    #do reweight for PSF
    minPsfr =   histInfo['gridPsf'][:-1]
    maxPsfr =   histInfo['gridPsf'][1:]
    ratios  =   np.array(histInfo['histPsfObs']/histInfo['histPsfSim1'])
    

    for ts in t:
        fwhm = ts['rpsf_det']*0.168*2.355
        if fwhm<=minPsfr[0] or fwhm>maxPsfr[-1]:
            ts['weightPsf']    =   0.
        else:
            ind         =   np.where((fwhm>minPsfr) & (fwhm<maxPsfr))[0][0]
            ratio   =   ratios[ind]
            if np.isnan(ratio) or abs(ratio)==np.inf:
                ts['weightPsf']=10
            else:
                ts['weightPsf']=ratio
    #apply the noise weight
    weight  =   t['weightNoi']*t['weightPsf'] 
    histNoi =   np.histogram(varNois,bins=gridNoi,weights=weight,density=True)[0]
    histNoi =   histNoi/sum(histNoi)
    histPsf =   np.histogram(fwhms,bins=gridPsf,weights=weight,density=True)[0]
    histPsf =   histPsf/sum(histPsf)
    histInfo['histNoiSim2']=histNoi
    histInfo['histPsfSim2']=histPsf
    return histInfo

def getHistObs(histInfo):
    pool    =   Pool(6)
    fields  =   np.load('./dr1FieldTract.npy').item().keys() 
    outDat  =   pool.map(getHistObsSub,fields)
    outDat  =   [x for x in outDat if x is not None]
    ngrid   =   20
    if len(outDat)<2:
        return
    hg      =   np.zeros((2,ngrid-1))
    for dat in outDat:
        hg+=dat
    for i in range(2):
        hg[i,:]/= sum(hg[i,:])
    histNoi =   hg[0]
    histPsf =   hg[1]
    histInfo=   np.load('noi_psf_histGram.npy').item()
    histInfo['histNoiObs']=histNoi
    histInfo['histPsfObs']=histPsf
    return histInfo

def getHistObsSub(ifield):
    if ifield=='AEGIS':
        return
    ngrid   =   20
    #Do mask
    directory= '/work/xiangchong.li/work/S16AFPFS/final/noBg/field'
    data    =   fitsio.read(os.path.join(directory,'%s.fits'%(ifield)))
    mask    =   np.ones(len(data)).astype(bool)
    mask    =   mask&(~np.isnan(data['noiVar']))
    mask    =   mask&(~np.isnan(data['psfMxx']))
    mask    =   mask&(~np.isnan(data['psfMyy']))
    mask    =   mask&(~np.isnan(data['psfMxy']))
    data    =   data[mask]
    fwhms   =   (data['psfMxx']*data['psfMyy'])**(0.25)*0.168*2.355
    varNois =   data['noiVar']
    histInfo=   np.load('noi_psf_histGram.npy').item()
    gridNoi =   histInfo['gridNoi']
    gridPsf =   histInfo['gridPsf']
    histPsf =   np.histogram(fwhms,bins=gridPsf,density=True)[0]
    histNoi =   np.histogram(varNois,bins=gridNoi,density=True)[0]
    return np.vstack([histNoi,histPsf])

def plotNoiPsf(histInfo):
    gridNoi =   histInfo['gridNoi']
    histNoiReWei =   histInfo['histNoiSim2']
    histNoiSim =   histInfo['histNoiSim0']
    histNoiObs =   histInfo['histNoiObs']
    nn      =   len(gridNoi)-1
    gridNoiPlot =   np.zeros(nn)
    for i in range(nn):
        gridNoiPlot[i]=(gridNoi[i]+gridNoi[i+1])/2.
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(gridNoiPlot,histNoiReWei,'-',c='r')
    ax.plot(gridNoiPlot,histNoiSim,'--',label='sim',c='r')
    ax.plot(gridNoiPlot,histNoiObs,'--',label='obs',c='b')
    ax.legend()
    plt.xlabel('noise variance')
    plt.ylabel('prob')
    plt.savefig('noiHist.png')
    
    
    gridPsf =   histInfo['gridPsf']
    histPsfReWei =   histInfo['histPsfSim2']
    histPsfSim =   histInfo['histPsfSim0']
    histPsfObs =   histInfo['histPsfObs']
    nn      =   len(gridPsf)-1
    gridPsfPlot =   np.zeros(nn)
    for i in range(nn):
        gridPsfPlot[i]=(gridPsf[i]+gridPsf[i+1])/2.
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(gridPsfPlot, histPsfReWei,'-',c='r')
    ax.plot(gridPsfPlot, histPsfSim,'--',label='sim',c='r')
    ax.plot(gridPsfPlot, histPsfObs,'--',label='obs',c='b')
    ax.legend()
    plt.xlabel('PSF FWHM [arcsec]')
    plt.ylabel('prob')
    plt.savefig('psfHist.png')
    return

if __name__ == '__main__':
    outFname=   'noiVar_rPsf_s18a.csv'
    if not os.path.exists(outFname):
        ndat    =   3200
        files   =   range(ndat)
        p       =   Pool(20)
        datarows=   p.map(getData,files)
        datarows=   [x for x in datarows if x is not None]
        names   =   ('num','field','var_target','var_meas','rpsf_trace','rpsf_det','weightNoi','weightPsf')
        t   =   Table(rows=datarows, names=names)
        t.write(outFname)
    else:
        t   =   Table.read(outFname)
    plotVarRatio(t)
    histOut =   'noi_psf_histGram.npy'
    if not os.path.exists(histOut):
        histInfo=   prepareHistInfo(t)
        histInfo=   getHistObs(histInfo)
        histInfo=   getHistSim(t,histInfo)
    else:
        histInfo=   np.load(histOut).item()
    plotNoiPsf(histInfo)
    t.write(outFname,overwrite=True)
    np.save('noi_psf_histGram.npy',histInfo)
