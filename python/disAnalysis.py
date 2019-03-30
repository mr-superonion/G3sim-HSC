#!/usr/bin/env python

import os
import sys
import glob
import optparse
import astropy.io.fits as pyfits
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


class disAnalysisConfig(pexConfig.Config):
    anaDir  =   pexConfig.Field(dtype=str, default='anaCat', doc = 'simulation catalogs directory')
    obsDir  =   pexConfig.Field(dtype=str, default='/work/xiangchong.li/work/S16AFPFS/final/noBg/field-matched/', doc = 'observed catalogs directory')
    pltDir  =   pexConfig.Field(dtype=str, default='plot', doc = 'plot directory')
    pltType =   pexConfig.Field(dtype=str, default='eps', doc = 'plot directory')
    version =   pexConfig.Field(dtype=str, default='v1', doc = 'input catalogs')
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

class disAnalysisRunner(TaskRunner):
    @staticmethod
    def getTargetList(parsedCmd, **kwargs):
        return [(ref, kwargs) for ref in utilities.field_keys.keys()] 

def unpickle(factory, args, kwargs):
    """Unpickle something by calling a factory"""
    return factory(*args, **kwargs)

class disAnalysisTask(BatchPoolTask):
    ConfigClass     =   disAnalysisConfig
    RunnerClass     =   disAnalysisRunner
    _DefaultName    =   "disAnalysis"
    
    def __init__(self,**kwargs):
        BatchPoolTask.__init__(self, **kwargs)
        self.nbin1      =   20
        self.nbin2      =   60
        self.varPsfnames=   ('ifield','rpsf_trace','rpsf_det','epsf1','epsf2','var_meas')
        self.disnames   =   ('gridNoi','histNoi','gridPsf','histPsf')
        self.fwhm_bin   =   np.linspace(0.3,0.9,self.nbin1)
        self.fwhm_bin_0 =   self.fwhm_bin[:-1]
        self.fwhm_bin_1 =   self.fwhm_bin[1:]
        self.var_bin    =   np.linspace(0.001,0.01,self.nbin1)
        self.var_bin_0  =   self.var_bin[:-1]
        self.var_bin_1  =   self.var_bin[1:]
        self.resolution_bin=np.linspace(0.1,1.,self.nbin2)
        self.resolution_bin_0  =   self.resolution_bin[:-1]
        self.resolution_bin_1  =   self.resolution_bin[1:]
        self.snr_bin    =   np.linspace(0.,80.,self.nbin2)
        self.snr_bin_0  =   self.snr_bin[:-1]
        self.snr_bin_1  =   self.snr_bin[1:]
        self.mag_bin    =   np.linspace(22.,25.2,self.nbin2)
        self.absE_bin   =   np.linspace(0.,2.,self.nbin2)
    def __reduce__(self):
        """Pickler"""
        return unpickle, (self.__class__, [], dict(config=self.config, name=self._name,
                parentTask=self._parentTask, log=self.log))
    
    @abortOnError
    def run(self,fieldKey):
        if fieldKey=='aegis':
            return
        fieldname   =   utilities.field_keys[fieldKey]
        self.log.info('Processing field: %s' %(fieldname))
        anaFieldDir =   os.path.join(self.config.anaDir,'%s' %(fieldname))
        outFname    =   os.path.join(anaFieldDir,'varPSF.csv')
        if not os.path.exists(outFname):
            anafnames   =   glob.glob('%s/*.fits' %(anaFieldDir)) 
            #Prepare the pool
            pool        =   Pool("disAnalysis")
            pool.cacheClear()
            rows        =   pool.map(self.getSimVarPsf,anafnames)
            rows        =   [x for x in rows if x is not None]
            t           =   Table(rows=rows, names=self.varPsfnames)
            t.sort('ifield')
            t.write(outFname,overwrite=True)
        else:
            t   =   Table.read(outFname)
        fwhms   =   t['rpsf_det']*0.168*2.355
        varNois =   t['var_meas']
        histNoi =   np.histogram(varNois,bins=self.var_bin,density=True)[0]
        histPsf =   np.histogram(fwhms,bins=self.fwhm_bin,density=True)[0]
        gridNoi =   (self.var_bin_0+self.var_bin_1)/2.
        gridPsf =   (self.fwhm_bin_0+self.fwhm_bin_1)/2.
        outDname=   os.path.join(anaFieldDir,'varPSF_dis_sim.csv')
        t_dis   =   Table(data=[gridNoi,histNoi,gridPsf,histPsf],\
                    names=self.disnames)
        t_dis.write(outDname,overwrite=True)
        outDnameObs=    os.path.join(anaFieldDir,'varPSF_dis_obs.csv')
        obsFame =   os.path.join(self.config.obsDir,'%s.fits' %(fieldname))
        t_dis_obs=  self.getObsHist(obsFname)
        t_dis_obs.write(outDnameObs,overwrite=True)
        self.plotNoiPsf(t_dis,t_dis_obs,anaFieldDir)
        return
    
    def getSimVarPsf(self,cache,anafname):
        #Get the noise variance and psf information
        ifield  =   anafname.split('-')[1].split('.')[0]
        data    =   pyfits.getdata(anafname)
        mask    =   np.ones(len(data)).astype(bool)
        mask    =   mask&(~np.isnan(data['base_SdssShape_psf_xx']))
        mask    =   mask&(~np.isnan(data['base_SdssShape_psf_yy']))
        mask    =   mask&(~np.isnan(data['base_SdssShape_psf_xy']))
        mask    =   mask&(~np.isnan(data['base_Variance_value']))
        ss      =   data[mask][0]
        rpsf_trace  =   np.sqrt((ss['base_SdssShape_psf_xx']+ss['base_SdssShape_psf_yy'])/2.)
        rpsf_det    =   (ss['base_SdssShape_psf_xx']*ss['base_SdssShape_psf_yy']-ss['base_SdssShape_psf_xy']**2.)**(0.25)
        epsf1   =   (ss['base_SdssShape_psf_xx']-ss['base_SdssShape_psf_yy'])/(ss['base_SdssShape_psf_xx']+ss['base_SdssShape_psf_yy'])
        epsf2   =   2.*ss['base_SdssShape_psf_xy']/(ss['base_SdssShape_psf_xx']+ss['base_SdssShape_psf_yy'])
        var_meas=   ss['base_Variance_value']
        return (ifield,rpsf_trace,rpsf_det,epsf1,epsf2,var_meas)
      
    def getObsHist(self,obsFname):
        data    =   pyfits.getdata(obsFname)
        mask    =   np.ones(len(data)).astype(bool)
        mask    =   mask&(~np.isnan(data['noiVar']))
        mask    =   mask&(~np.isnan(data['psfMxx']))
        mask    =   mask&(~np.isnan(data['psfMyy']))
        mask    =   mask&(~np.isnan(data['psfMxy']))
        data    =   data[mask]
        fwhms   =   (data['psfMxx']*data['psfMyy']-data['psfMxy']**2)**(0.25)*0.168*2.355
        varNois =   data['noiVar']
        gridNoi =   (self.var_bin_0+self.var_bin_1)/2.
        gridPsf =   (self.fwhm_bin_0+self.fwhm_bin_1)/2.
        histNoi =   np.histogram(varNois,bins=self.var_bin,density=True)[0]
        histPsf =   np.histogram(fwhms,bins=self.fwhm_bin,density=True)[0] 
        return Table(data=[gridNoi,histNoi,gridPsf,histPsf],\
                    names=self.disnames)
    
    
    def plotNoiPsf(self,t_dis,t_dis_obs,anaFieldDir):
        
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(t_dis['gridNoi'],t_dis['histNoi'],'--',label='sim',c='r')
        ax.plot(t_dis_obs['gridNoi'],t_dis_obs['histNoi'],'--',label='obs',c='b')
        ax.legend()
        plt.xlabel('noise variance')
        plt.ylabel('prob')
        plt.savefig(os.path.join(anaFieldDir,'noiHist.png'))

        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(t_dis['gridPsf'],t_dis['histPsf'],'--',label='sim',c='r')
        ax.plot(t_dis_obs['gridPsf'],t_dis_obs['histPsf'],'--',label='obs',c='b')
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
