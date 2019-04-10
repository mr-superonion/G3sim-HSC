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
        self.nbin2      =   60
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
            obsFname    =   os.path.join(anaFieldDir,'%s.fits' %(fieldname))
            dataObs     =   fitsio.read(obsFname)
            t_dis_obs   =   self.getObsHist(dataObs)
            t_dis_obs.write(outDnameObs,overwrite=True)
        else:
            t_dis_obs   =   Table.read(outDnameObs)
            
        self.log.info('Achieving the distriutions for simulations')
        outDname    =   os.path.join(anaFieldDir,'varPSF_dis_sim.csv')
        dataSim     =   fitsio.read(simFname)
        if not os.path.exists(outDname):
            simFname    =   os.path.join(anaFieldDir,'%s.fits' %(fieldname))
            t_dis_sim   =   self.getObsHist(dataSim)
            t_dis_sim.write(outDname,overwrite=True)
        else:
            t_dis_sim   =   Table.read(outDname)
        
        ratios_noi  =   np.array(t_dis_obs['histNoi']/t_dis_sim['histNoi'])
        """            
        outDnameRw  =   os.path.join(anaFieldDir,'varPSF_dis_rw.csv')
        t_dis_rw    =   self.reWeight(dataSim,t_dis_sim,t_dis_obs)
        t_dis_rw.write(outDnameRw,overwrite=True)
        self.plotNoiPsf(t_dis_sim,t_dis_obs,t_dis_rw,anaFieldDir)
        """
        return
    
    def getObsHist(self,data):
        mask    =   np.ones(len(data)).astype(bool)
        mask    =   mask&(~np.isnan(data['noiVar']))
        mask    =   mask&(~np.isnan(data['psfMxx']))
        mask    =   mask&(~np.isnan(data['psfMyy']))
        mask    =   mask&(~np.isnan(data['psfMxy']))
        data    =   data[mask]
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
    
    
    
    def getSimHist(self,data):
        #Get the noise variance and psf information
        mask    =   np.ones(len(data)).astype(bool)
        mask    =   mask&(~np.isnan(data['ext_shapeHSM_HsmPsfMoments_xx']))
        mask    =   mask&(~np.isnan(data['ext_shapeHSM_HsmPsfMoments_yy']))
        mask    =   mask&(~np.isnan(data['ext_shapeHSM_HsmPsfMoments_xy']))
        mask    =   mask&(~np.isnan(data['base_Variance_value']))
        mask    =   mask&(data['weak_lensing_flag'])
        data    =   data[mask]
        rpsf_trace  =   np.sqrt((data['ext_shapeHSM_HsmPsfMoments_xx']+data['ext_shapeHSM_HsmPsfMoments_yy'])/2.)
        rpsf_det    =   (data['ext_shapeHSM_HsmPsfMoments_xx']*data['ext_shapeHSM_HsmPsfMoments_yy']-data['ext_shapeHSM_HsmPsfMoments_xy']**2.)**(0.25)
        epsf1   =   (data['ext_shapeHSM_HsmPsfMoments_xx']-data['ext_shapeHSM_HsmPsfMoments_yy'])/(data['ext_shapeHSM_HsmPsfMoments_xx']+data['ext_shapeHSM_HsmPsfMoments_yy'])
        epsf2   =   2.*data['ext_shapeHSM_HsmPsfMoments_xy']/(data['ext_shapeHSM_HsmPsfMoments_xx']+data['ext_shapeHSM_HsmPsfMoments_yy'])
        varNois =   data['base_Variance_value']
        fwhms   =   rpsf_det*0.168*2.355
        histNoi =   np.histogram(varNois,bins=self.var_bin,density=True)[0]
        histPsf =   np.histogram(fwhms,bins=self.fwhm_bin,density=True)[0]
        return Table(data=[self.gridNoi,histNoi,self.gridPsf,histPsf,len(data)*np.ones(len(self.gridNoi))],\
                    names=self.disnames)
      
    
    def reWeight(self,t,t_dis,t_dis_obs):
        #reweight for noise
       
        weightNoi=  np.zeros(len(t))
        for ii,ts in enumerate(t):
            if ts['var_meas']<=self.var_bin[0] or ts['var_meas']>self.var_bin[-1]:
                weightNoi[ii]   =   0.
            else:
                ind             =   np.where((ts['var_meas']>self.var_bin_0) & (ts['var_meas']<self.var_bin_1))[0][0]
                weightNoi[ii]   =   ratios[ind]
        t['weightNoi']  =   weightNoi
        #apply the noise weight
        fwhms       =   t['rpsf_det']*0.168*2.355
        varNois     =   t['var_meas']
        histNoi_rn  =   np.histogram(varNois,bins=self.var_bin,weights=weightNoi,density=True)[0]
        histPsf_rn  =   np.histogram(fwhms,bins=self.fwhm_bin,weights=weightNoi,density=True)[0]
        #reweight for PSF
        ratios  =   np.array(t_dis_obs['histPsf']/histPsf_rn)
        weightPsf=  np.zeros(len(t))
        for ii,ts in enumerate(t):
            fwhm = ts['rpsf_det']*0.168*2.355
            if fwhm<=self.fwhm_bin[0] or fwhm>self.fwhm_bin[-1]:
                weightPsf[ii] =   0.
            else:
                ind     =   np.where((fwhm>self.fwhm_bin_0) & (fwhm<self.fwhm_bin_1))[0][0]
                ratio   =   ratios[ind]
                if np.isnan(ratio) or abs(ratio)==np.inf:
                    weightPsf[ii]   =   0.
                else:
                    weightPsf[ii]   =   ratio
        t['weightPsf']  =   weightPsf
        #apply the noise weight
        weight  =   weightNoi*weightPsf 
        t['weight'] =   weight
        histNoi_rp =   np.histogram(varNois,bins=self.var_bin,weights=weight,density=True)[0]
        histPsf_rp =   np.histogram(fwhms,bins=self.fwhm_bin,weights=weight,density=True)[0]
        return Table(data=[self.gridNoi,histNoi_rp,self.gridPsf,histPsf_rp,len(t)*np.ones(len(self.gridNoi))],\
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
