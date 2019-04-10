#!/usr/bin/env python

import os
import sys
import glob
import fitsio
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


class psfNoiDisAnalysisConfig(pexConfig.Config):
    anaDir  =   pexConfig.Field(dtype=str, default='anaCat', doc = 'simulation catalogs directory')
    obsDir  =   pexConfig.Field(dtype=str, default='/work/xiangchong.li/work/S16AFPFS/final/noBg/field-matched/', doc = 'observed catalogs directory')
    pltDir  =   pexConfig.Field(dtype=str, default='plot', doc = 'plot directory')
    def setDefaults(self):
        pexConfig.Config.setDefaults(self)

    def validate(self):
        pexConfig.Config.validate(self)
       # Checks that the paths provided are valid
        if not os.path.isdir(self.anaDir):
            raise RuntimeError("Simulation catalog directory %s does not exist "
                         "or is not a directory" % (self.anaDir))
        if not os.path.isdir(self.obsDir):
            raise RuntimeError("Data catalog directory %s does not exist "
                         "or is not a directory" % (self.anaDir))
        if not os.path.exists(self.pltDir):
            os.mkdir(self.plotDir)

class psfNoiDisAnalysisRunner(TaskRunner):
    @staticmethod
    def getTargetList(parsedCmd, **kwargs):
        return [(ref, kwargs) for ref in range(1)] 

def unpickle(factory, args, kwargs):
    """Unpickle something by calling a factory"""
    return factory(*args, **kwargs)

class psfNoiDisAnalysisTask(BatchPoolTask):
    ConfigClass     =   psfNoiDisAnalysisConfig
    RunnerClass     =   psfNoiDisAnalysisRunner
    _DefaultName    =   "psfNoiDisAnalysis"
    
    def __init__(self,**kwargs):
        BatchPoolTask.__init__(self, **kwargs)
        self.nbin1      =   20
        self.varPsfnames=   ('ifield','rpsf_trace','rpsf_det','epsf1','epsf2','var_meas')
        self.disnames   =   ('gridNoi','histNoi','gridPsf','histPsf','num')
        self.fwhm_bin   =   np.linspace(0.3,0.9,self.nbin1)
        self.fwhm_bin_0 =   self.fwhm_bin[:-1]
        self.fwhm_bin_1 =   self.fwhm_bin[1:]
        self.gridPsf    =   (self.fwhm_bin_0+self.fwhm_bin_1)/2.
        self.var_bin    =   np.linspace(0.001,0.01,self.nbin1)
        self.var_bin_0  =   self.var_bin[:-1]
        self.var_bin_1  =   self.var_bin[1:]
        self.gridNoi    =   (self.var_bin_0+self.var_bin_1)/2.
    def __reduce__(self):
        """Pickler"""
        return unpickle, (self.__class__, [], dict(config=self.config, name=self._name,
                parentTask=self._parentTask, log=self.log))
    
    @abortOnError
    
    def run(self,Id):
        #Prepare the pool
        pool        =   Pool("psfNoiDisAnalysis")
        pool.cacheClear()
        pool.map(self.process,utilities.field_keys.keys())
        return
    
    def process(self,cache,fieldKey):
        if fieldKey=='aegis':
            return
        fieldname   =   utilities.field_keys[fieldKey]
        self.log.info('Processing field: %s' %(fieldname))
        anaFieldDir =   os.path.join(self.config.anaDir,'%s' %(fieldname))
        
        self.log.info('Achieving the distriutions for observations')
        outDnameObs =   os.path.join(anaFieldDir,'varPSF_dis_obs.csv')
        if not os.path.exists(outDnameObs):
            obsFname    =   os.path.join(self.config.obsDir,'%s.fits' %(fieldname))
            dataObs     =   fitsio.read(obsFname)
            t_dis_obs   =   self.getObsHist(dataObs)
            t_dis_obs.write(outDnameObs,overwrite=True)
        else:
            t_dis_obs   =   Table.read(outDnameObs)
            
        self.log.info('Achieving the distriutions for simulations')
        outDname    =   os.path.join(anaFieldDir,'varPSF_dis_sim.csv')
        simFname    =   os.path.join(anaFieldDir,'%s.fits' %(fieldname))
        dataSim     =   fitsio.read(simFname)
        dataSim     =   self.maskSim(dataSim)
        if not os.path.exists(outDname):
            t_dis_sim   =   self.getSimHist(dataSim)
            t_dis_sim.write(outDname,overwrite=True)
        else:
            t_dis_sim   =   Table.read(outDname)
         
        self.log.info('Reweighting the simulations')
        outDnameRw      =   os.path.join(anaFieldDir,'varPSF_dis_rw.csv')
        if not os.path.exists(outDnameRw):
            t_dis_rw    =   self.reWeight(dataSim,t_dis_sim,t_dis_obs)
            t_dis_rw.write(outDnameRw,overwrite=True)
        else:
            t_dis_rw    =   Table.read(outDnameRw)   
        self.plotNoiPsf(t_dis_sim,t_dis_obs,t_dis_rw,anaFieldDir)
        return
    
    def getObsHist(self,data):
        mask    =   np.ones(len(data)).astype(bool)
        mask    =   mask&(~np.isnan(data['noiVar']))
        mask    =   mask&(~np.isnan(data['psfMxx']))
        mask    =   mask&(~np.isnan(data['psfMyy']))
        mask    =   mask&(~np.isnan(data['psfMxy']))
        mask    =   mask&((data['CModel_mag'])<24.5)
        mask    =   mask&(data['extendness']>0)
        mask    =   mask&(data['ishape_hsm_regauss_sigma']>0.)
        mask    =   mask&(data['ishape_hsm_regauss_sigma']<0.4)
        mask    =   mask&(data['reGaus_resolution']>0.3)
        mask    =   mask&(data['CModel_snr']>10.)
        mask    =   mask&(data['reGaus_sigmaE']<2.)
        data    =   data[mask]
        fwhms   =   (data['psfMxx']*data['psfMyy']-data['psfMxy']**2)**(0.25)*0.168*2.355
        varNois =   data['noiVar']
        histNoi =   np.histogram(varNois,bins=self.var_bin,density=True)[0]
        histPsf =   np.histogram(fwhms,bins=self.fwhm_bin,density=True)[0] 
        return Table(data=[self.gridNoi,histNoi,self.gridPsf,histPsf,len(data)*np.ones(len(self.gridNoi))],\
                    names=self.disnames)
    
    def maskSim(self,data): 
        mask    =   np.ones(len(data)).astype(bool)
        mask    =   mask&(data['weak_lensing_flag'])
        data    =   data[mask]
        return data
    
    def getSimHist(self,data):
        #Get the noise variance and psf information
        varNois =   data['base_Variance_value']
        rpsf_det=   (data['ext_shapeHSM_HsmPsfMoments_xx']*data['ext_shapeHSM_HsmPsfMoments_yy']-data['ext_shapeHSM_HsmPsfMoments_xy']**2.)**(0.25)
        fwhms   =   rpsf_det*0.168*2.355
        epsf1   =   (data['ext_shapeHSM_HsmPsfMoments_xx']-data['ext_shapeHSM_HsmPsfMoments_yy'])/(data['ext_shapeHSM_HsmPsfMoments_xx']+data['ext_shapeHSM_HsmPsfMoments_yy'])
        epsf2   =   2.*data['ext_shapeHSM_HsmPsfMoments_xy']/(data['ext_shapeHSM_HsmPsfMoments_xx']+data['ext_shapeHSM_HsmPsfMoments_yy'])
        rpsf_trace  =   np.sqrt((data['ext_shapeHSM_HsmPsfMoments_xx']+data['ext_shapeHSM_HsmPsfMoments_yy'])/2.)
        histNoi =   np.histogram(varNois,bins=self.var_bin,density=True)[0]
        histPsf =   np.histogram(fwhms,bins=self.fwhm_bin,density=True)[0]
        return Table(data=[self.gridNoi,histNoi,self.gridPsf,histPsf,len(data)*np.ones(len(self.gridNoi))],\
                    names=self.disnames)
      
    
    def reWeight(self,data,t_dis,t_dis_obs):
        varNois =   data['base_Variance_value']
        rpsf_det=   (data['ext_shapeHSM_HsmPsfMoments_xx']*data['ext_shapeHSM_HsmPsfMoments_yy']-data['ext_shapeHSM_HsmPsfMoments_xy']**2.)**(0.25)
        fwhms   =   rpsf_det*0.168*2.355
        #reweight for noise
        ratio_noi=   np.array(t_dis_obs['histNoi']/t_dis['histNoi'])
        weightNoi=  np.zeros(len(data))
        for ii,vn in enumerate(varNois):
            if vn<=self.var_bin[0] or vn>self.var_bin[-1]:
                weightNoi[ii]   =   0.
            else:
                ind             =   np.where((vn>self.var_bin_0) & (vn<self.var_bin_1))[0][0]
                weightNoi[ii]   =   ratio_noi[ind]
        #apply the noise weight
        histNoi_rn  =   np.histogram(varNois,bins=self.var_bin,weights=weightNoi,density=True)[0]
        histPsf_rn  =   np.histogram(fwhms,bins=self.fwhm_bin,weights=weightNoi,density=True)[0]
        #reweight for PSF
        ratios  =   np.array(t_dis_obs['histPsf']/histPsf_rn)
        weightPsf=  np.zeros(len(data))
        for ii,rp in enumerate(fwhms):
            if rp<=self.fwhm_bin[0] or rp>self.fwhm_bin[-1]:
                weightPsf[ii] =   0.
            else:
                ind     =   np.where((rp>self.fwhm_bin_0) & (rp<self.fwhm_bin_1))[0][0]
                ratio   =   ratios[ind]
                if np.isnan(ratio) or abs(ratio)==np.inf:
                    weightPsf[ii]   =   10.
                else:
                    weightPsf[ii]   =   ratio
        data['reweight_s16a']  =   weightPsf*weightNoi
        histNoi_rp =   np.histogram(varNois,bins=self.var_bin,weights=data['reweight_s16a'],density=True)[0]
        histPsf_rp =   np.histogram(fwhms,bins=self.fwhm_bin,weights=data['reweight_s16a'],density=True)[0]
        return Table(data=[self.gridNoi,histNoi_rp,self.gridPsf,histPsf_rp,len(data)*np.ones(len(self.gridNoi))],\
                    names=self.disnames)
    
    def plotNoiPsf(self,t_dis,t_dis_obs,t_dis_rw,anaFieldDir):
        fig = plt.figure()
        ax = fig.add_subplot(111)
        #ax.plot(t_dis['gridNoi'],t_dis['histNoi'],'--',label='sim',c='r')
        ax.plot(t_dis_obs['gridNoi'],t_dis_obs['histNoi'],'-',label='obs',c='b')
        ax.plot(t_dis_rw['gridNoi'],t_dis_rw['histNoi'],'--',label='reweight',c='r')
        ax.legend()
        plt.xlabel('noise variance')
        plt.ylabel('prob')
        plt.savefig(os.path.join(anaFieldDir,'noiHist.png'))

        fig = plt.figure()
        ax = fig.add_subplot(111)
        #ax.plot(t_dis['gridPsf'],t_dis['histPsf'],'--',label='sim',c='r')
        ax.plot(t_dis_obs['gridPsf'],t_dis_obs['histPsf'],'-',label='obs',c='b')
        ax.plot(t_dis_rw['gridPsf'],t_dis_rw['histPsf'],'--',label='reweight',c='r')
        ax.legend()
        plt.xlabel('PSF FWHM [arcsec]')
        plt.ylabel('prob')
        plt.savefig(os.path.join(anaFieldDir,'psfHist.png'))
        return
    
    @classmethod
    def _makeArgumentParser(cls, *args, **kwargs):
        kwargs.pop("doBatch", False)
        parser = pipeBase.ArgumentParser(name=cls._DefaultName)
        return parser
    
    @classmethod
    def batchWallTime(cls, time, parsedCmd, numCpus):
        return None

    def writeConfig(self, butler, clobber=False, doBackup=False):
        pass

    def writeSchemas(self, butler, clobber=False, doBackup=False):
        pass
    
    def writeMetadata(self, ifield):
        pass
    
    def writeEupsVersions(self, butler, clobber=False, doBackup=False):
        pass
    
    def _getConfigName(self):
        return None
   
    def _getEupsVersionsName(self):
        return None
    
    def _getMetadataName(self):
        return None
