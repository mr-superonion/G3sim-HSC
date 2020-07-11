#!/usr/bin/env python
import os
import gc
import numpy as np
import healpy as hp
import astropy.io.fits as pyfits
from astropy.table import Table
import astropy.io.ascii as ascii
import sim_analysis_utilities as utilities

import matplotlib as mpl
mpl.use('agg')
mpl.rcParams['ps.useafm'] = True
mpl.rcParams['pdf.use14corefonts'] = True
mpl.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt

anaDir  =   'anaCat_newS19Mask/'
obsDir  =   '../../database/s19-Hironao-2/fields/'
pltDir  =   'plot/psfNoi/'
if not os.path.isdir(anaDir):
    raise RuntimeError("Simulation catalog directory %s does not exist "
                 "or is not a directory" % (anaDir))
if not os.path.isdir(obsDir):
    raise RuntimeError("Data catalog directory %s does not exist "
                 "or is not a directory" % (anaDir))
if not os.path.exists(pltDir):
    os.system('mkdir -p %s'%pltDir)

nbin1       =   20
disnames    =   ('gridNoi','histNoi','gridPsf','histPsf','num')
fwhm_bin    =   np.linspace(0.3,0.9,nbin1)
fwhm_bin_0  =   fwhm_bin[:-1]
fwhm_bin_1  =   fwhm_bin[1:]
gridPsf     =   (fwhm_bin_0+fwhm_bin_1)/2.
var_bin     =   np.linspace(0.001,0.01,nbin1)
var_bin_0   =   var_bin[:-1]
var_bin_1   =   var_bin[1:]
gridNoi     =   (var_bin_0+var_bin_1)/2.
colors      =   ['blue','red']

def main():
    for fname in utilities.field_keys.values():
        process(fname)
    finalPlot()
    return

def process(fieldname):
    anaFieldDir =   os.path.join(anaDir,'%s' %(fieldname))
    outDnameRw  =   os.path.join(anaFieldDir,'varPSF_dis_rw.csv')
    if os.path.exists(outDnameRw):
        return

    simFname    =   os.path.join(anaFieldDir,'%s.fits' %(fieldname))
    dataSim     =   pyfits.getdata(simFname)
    dataSim['reweight_s16a']    =   0.

    print('Achieving the distriutions for observations')
    outDnameObs =   os.path.join(anaFieldDir,'varPSF_dis_obs.csv')
    if not os.path.exists(outDnameObs):
        obsFname    =   os.path.join(obsDir,'%s.fits' %(fieldname))
        dataObs     =   pyfits.getdata(obsFname)
        t_dis_obs   =   getObsHist(dataObs)
        t_dis_obs.write(outDnameObs,overwrite=True)
    else:
        t_dis_obs   =   Table.read(outDnameObs)

    print('Achieving the distriutions for simulations')
    outDname    =   os.path.join(anaFieldDir,'varPSF_dis_sim.csv')
    if not os.path.exists(outDname):
        t_dis_sim   =   getSimHist(dataSim)
        t_dis_sim.write(outDname,overwrite=True)
    else:
        t_dis_sim   =   Table.read(outDname)

    print('Reweighting the simulations')
    if not os.path.exists(outDnameRw):
        t_dis_rw    =   reWeight(dataSim,t_dis_sim,t_dis_obs)
        t_dis_rw.write(outDnameRw,overwrite=True)
        pyfits.writeto(simFname,dataSim,overwrite=True)

    plotNoiPsf(t_dis_sim,t_dis_obs,t_dis_rw,anaFieldDir)
    return

def getObsHist(data):
    mask    =   np.ones(len(data)).astype(bool)
    mask    =   mask&(~np.isnan(data['i_variance_value']))
    mask    =   mask&(~np.isnan(data['i_sdssshape_psf_shape11']))
    mask    =   mask&(~np.isnan(data['i_sdssshape_psf_shape12']))
    mask    =   mask&(~np.isnan(data['i_sdssshape_psf_shape22']))
    data    =   data[mask]
    fwhms   =   utilities.get_psf_size(data)
    varNois =   utilities.get_noi_var(data)
    histNoi =   np.histogram(varNois,bins=var_bin,density=True)[0]
    histPsf =   np.histogram(fwhms,bins=fwhm_bin,density=True)[0]
    return Table(data=[gridNoi,histNoi,gridPsf,histPsf,len(data)*np.ones(len(gridNoi))],\
                names=disnames)

def getSimHist(data):
    #Get the noise variance and psf information
    data    =   utilities.update_wl_flags_sim(data)
    mask    =   data['weak_lensing_flag']
    data    =   data[mask]
    varNois =   utilities.get_noi_var(data)
    fwhms   =   utilities.get_psf_size(data)
    histNoi =   np.histogram(varNois,bins=var_bin,density=True)[0]
    histPsf =   np.histogram(fwhms,bins=fwhm_bin,density=True)[0]
    return Table(data=[gridNoi,histNoi,gridPsf,histPsf,len(data)*np.ones(len(gridNoi))],\
    names=disnames)

def reWeight(data,t_dis,t_dis_obs):
    varNois =   utilities.get_noi_var(data)
    fwhms   =   utilities.get_psf_size(data)
    #reweight for noise
    ratio_noi=   np.array(t_dis_obs['histNoi']/t_dis['histNoi'])
    weightNoi=  np.zeros(len(data))
    for ii,vn in enumerate(varNois):
        if vn<=var_bin[0] or vn>var_bin[-1]:
            weightNoi[ii]   =   0.
        else:
            ind             =   np.where((vn>var_bin_0) & (vn<var_bin_1))[0][0]
            weightNoi[ii]   =   ratio_noi[ind]
    #apply the noise weight
    histNoi_rn  =   np.histogram(varNois,bins=var_bin,weights=weightNoi,density=True)[0]
    histPsf_rn  =   np.histogram(fwhms,bins=fwhm_bin,weights=weightNoi,density=True)[0]

    #reweight for PSF
    ratios  =   np.array(t_dis_obs['histPsf']/histPsf_rn)
    weightPsf=  np.zeros(len(data))
    for ii,rp in enumerate(fwhms):
        if rp<=fwhm_bin[0] or rp>fwhm_bin[-1]:
            weightPsf[ii] =   0.
        else:
            ind     =   np.where((rp>fwhm_bin_0) & (rp<fwhm_bin_1))[0][0]
            ratio   =   ratios[ind]
            if np.isnan(ratio) or abs(ratio)==np.inf:
                weightPsf[ii]   =   0.
            else:
                weightPsf[ii]   =   ratio
    reweight    =   weightPsf*weightNoi
    data['reweight_s16a']   =   reweight
    histNoi_rp  =   np.histogram(varNois,bins=var_bin,weights=data['reweight_s16a'],density=True)[0]
    histPsf_rp  =   np.histogram(fwhms,bins=fwhm_bin,weights=data['reweight_s16a'],density=True)[0]
    return Table(data=[gridNoi,histNoi_rp,gridPsf,histPsf_rp,len(data)*np.ones(len(gridNoi))],\
    names=disnames)

def plotNoiPsf(t_dis,t_dis_obs,t_dis_rw,anaFieldDir):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(t_dis['gridNoi'],t_dis['histNoi'],'--',label='sim',c='y')
    ax.plot(t_dis_obs['gridNoi'],t_dis_obs['histNoi'],'-',label='obs',c='b')
    ax.plot(t_dis_rw['gridNoi'],t_dis_rw['histNoi'],'--',label='sim',c='r')
    ax.legend()
    plt.xlabel('noise variance')
    plt.ylabel('prob')
    plt.savefig(os.path.join(anaFieldDir,'noiHist.png'))

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(t_dis['gridPsf'],t_dis['histPsf'],'--',label='sim',c='y')
    ax.plot(t_dis_obs['gridPsf'],t_dis_obs['histPsf'],'-',label='obs',c='b')
    ax.plot(t_dis_rw['gridPsf'],t_dis_rw['histPsf'],'--',label='sim',c='r')
    ax.legend()
    plt.xlabel('PSF FWHM [arcsec]')
    plt.ylabel('prob')
    plt.savefig(os.path.join(anaFieldDir,'psfHist.png'))
    return

def finalPlot():
    obsDis_var  =   np.zeros(nbin1-1)
    obsDis_fwhm =   np.zeros(nbin1-1)
    obsNum      =   np.zeros(nbin1-1)

    simDis_var  =   np.zeros(nbin1-1)
    simDis_fwhm =   np.zeros(nbin1-1)
    simNum      =   np.zeros(nbin1-1)

    simRDis_var  =   np.zeros(nbin1-1)
    simRDis_fwhm =   np.zeros(nbin1-1)
    simRNum      =   np.zeros(nbin1-1)

    figN=   plt.figure()
    figP=   plt.figure()
    """
    Plot for 6 fields
    """
    for ii,fieldKey in enumerate(utilities.field_keys):
        fieldname   =   utilities.field_keys[fieldKey]
        anaFieldDir =   os.path.join(anaDir,'%s' %(fieldname))

        disnameOri  =   os.path.join(anaFieldDir,'varPSF_dis_sim.csv')
        t_dis       =   Table.read(disnameOri)
        simDis_var  +=  t_dis['histNoi']*t_dis['num']
        simDis_fwhm +=  t_dis['histPsf']*t_dis['num']
        simNum      +=  t_dis['num']

        disnameObs  =    os.path.join(anaFieldDir,'varPSF_dis_obs.csv')
        t_dis_obs   =   Table.read(disnameObs)
        obsDis_var  +=  t_dis_obs['histNoi']*t_dis_obs['num']
        obsDis_fwhm +=  t_dis_obs['histPsf']*t_dis_obs['num']
        obsNum      +=  t_dis_obs['num']

        disnameRw   =   os.path.join(anaFieldDir,'varPSF_dis_rw.csv')
        t_dis_rw    =   Table.read(disnameRw)
        simRDis_var +=  t_dis_rw['histNoi']*t_dis_rw['num']
        simRDis_fwhm+=  t_dis_rw['histPsf']*t_dis_rw['num']
        simRNum     +=  t_dis_rw['num']

        axN =   figN.add_subplot(3,3,ii+1)
        gg1 =   t_dis_obs['gridNoi']*100
        hh1 =   t_dis_obs['histNoi']
        gg2 =   t_dis_rw['gridNoi']*100
        hh2 =   t_dis_rw['histNoi']
        dg  =   gg1[1]-gg1[0]
        axN.plot(gg1,hh1,'-',label='obs',c='b')
        axN.plot(gg2,hh2,'--',label='sim',c='r')
        ave1=   np.sum(gg1*hh1)*dg/100.
        ave2=   np.sum(gg2*hh2)*dg/100.
        axN.axvline(x=ave1,ls='-',c='blue')
        axN.axvline(x=ave2,ls='--',c='red')
        ymin,ymax=axN.get_ylim()
        axN.text(0.8,0.8*ymax+0.2*ymin, '%.2f' %ave1,
            color='blue')
        axN.text(0.8,0.6*ymax+0.4*ymin , '%.2f' %ave2,
            color='red')
        axN.set_xlabel(r'noise variance [$\times 10^2$]')
        axN.set_title(fieldname)

        axP =   figP.add_subplot(3,3,ii+1)
        gg1 =   t_dis_obs['gridPsf']
        hh1 =   t_dis_obs['histPsf']
        gg2 =   t_dis_rw['gridPsf']
        hh2 =   t_dis_rw['histPsf']
        dg  =   gg1[1]-gg1[0]
        axP.plot(gg1,hh1,'-',label='obs',c='b')
        axP.plot(gg2,hh2,'--',label='sim',c='r')
        ave1=   np.sum(gg1*hh1)*dg
        ave2=   np.sum(gg2*hh2)*dg
        axP.axvline(x=ave1,ls='-',c='blue')
        axP.axvline(x=ave2,ls='--',c='red')

        ymin,ymax=axP.get_ylim()
        axP.text(0.8,0.8*ymax+0.2*ymin, '%.2f' %ave1,
            color='blue')
        axP.text(0.8,0.6*ymax+0.4*ymin , '%.2f' %ave2,
            color='red')
        axP.set_xlabel('PSF FWHM [arcsec]')
        axP.set_title(fieldname)
    """
    Plot for all fields
    """
    obsDis_var /=obsNum
    obsDis_fwhm/=obsNum
    simDis_var /=simNum
    simDis_fwhm/=simNum
    simRDis_var /=simRNum
    simRDis_fwhm/=simRNum

    axN =   figN.add_subplot(3,3,7)
    axN.plot(gridNoi*100,obsDis_var,'-',label='obs',c='b')
    axN.plot(gridNoi*100,simRDis_var,'--',label='sim',c='r')
    dg  =   (gridNoi[1]-gridNoi[0])*100.
    ave1=   sum(gridNoi*obsDis_var)*dg
    ave2=   sum(gridNoi*simRDis_var)*dg
    axN.axvline(x=ave1,ls='-',c='blue')
    axN.axvline(x=ave2,ls='--',c='red')
    ymin,ymax=axN.get_ylim()
    axN.text(0.8,0.8*ymax+0.2*ymin, '%.2f' %ave1,
        color='blue')
    axN.text(0.8,0.6*ymax+0.4*ymin , '%.2f' %ave2,
        color='red')

    axN.set_xlabel(r'noise variance [$\times 10^2$]')
    axN.set_title('All')
    figN.tight_layout()
    axN.legend(loc='center', bbox_to_anchor=(1.5, 0.5),fontsize=15)
    figN.savefig(os.path.join(pltDir,'noiVar2.eps'))

    axP =   figP.add_subplot(3,3,7)
    axP.plot(gridPsf,obsDis_fwhm,'-',label='obs',c='b')
    axP.plot(gridPsf,simRDis_fwhm,'--',label='sim',c='r')
    dg  =   gridPsf[1]-gridPsf[0]
    ave1=sum(gridPsf*obsDis_fwhm)*dg
    ave2=sum(gridPsf*simRDis_fwhm)*dg
    axP.axvline(x=ave1,ls='-',c='blue')
    axP.axvline(x=ave2,ls='--',c='red')
    ymin,ymax=axP.get_ylim()
    axP.text(0.8,0.8*ymax+0.2*ymin, '%.2f' %ave1,
        color='blue')
    axP.text(0.8,0.6*ymax+0.4*ymin , '%.2f' %ave2,
        color='red')
    axP.set_xlabel('PSF FWHM [arcsec]')
    axP.set_title('All')
    figP.tight_layout()
    axP.legend(loc='center', bbox_to_anchor=(1.5, 0.5),fontsize=15)
    figP.savefig(os.path.join(pltDir,'psfFWHM2.eps'))
    return

if __name__=='__main__':
    main()
