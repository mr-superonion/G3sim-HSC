#!/usr/bin/env python

import os
import sys
import glob
import fitsio
import optparse
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


class galDisAnalysisConfig(pexConfig.Config):
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

class galDisAnalysisRunner(TaskRunner):
    @staticmethod
    def getTargetList(parsedCmd, **kwargs):
        return [(ref, kwargs) for ref in range(1)] 

def unpickle(factory, args, kwargs):
    """Unpickle something by calling a factory"""
    return factory(*args, **kwargs)

class galDisAnalysisTask(BatchPoolTask):
    ConfigClass     =   galDisAnalysisConfig
    RunnerClass     =   galDisAnalysisRunner
    _DefaultName    =   "galDisAnalysis"
    
    def __init__(self,**kwargs):
        BatchPoolTask.__init__(self, **kwargs)
        self.nbin2      =   60
        self.resolution_bin=np.linspace(0.1,1.,self.nbin2)
        self.resolution_bin_0  =   self.resolution_bin[:-1]
        self.resolution_bin_1  =   self.resolution_bin[1:]
        self.resGrid    =   (self.resolution_bin_0+self.resolution_bin_1)/2.
        self.snr_bin    =   np.linspace(0.,80.,self.nbin2)
        self.snr_bin_0  =   self.snr_bin[:-1]
        self.snr_bin_1  =   self.snr_bin[1:]
        self.snrGrid    =   (self.snr_bin_0+self.snr_bin_1)/2.
        self.mag_bin    =   np.linspace(22.,25.2,self.nbin2)
        self.mag_bin_0  =   self.mag_bin[:-1]
        self.mag_bin_1  =   self.mag_bin[1:]
        self.magGrid    =   (self.mag_bin_0+self.mag_bin_1)/2.
        self.absE_bin   =   np.linspace(0.,2.,self.nbin2)
        self.absE_bin_0 =   self.absE_bin[:-1]
        self.absE_bin_1 =   self.absE_bin[1:]
        self.absEGrid   =   (self.absE_bin_0+self.absE_bin_1)/2.
        self.datanames   =   ['gridRes','histRes','gridSnr','histSnr','gridMag','histMag','gridAbsE','histAbsE','num']
    def __reduce__(self):
        """Pickler"""
        return unpickle, (self.__class__, [], dict(config=self.config, name=self._name,
                parentTask=self._parentTask, log=self.log))
    
    @abortOnError
    def run(self,Id)
        #Prepare the pool
        pool        =   Pool("galDisAnalysis")
        pool.cacheClear()
        pool.map(self.getSimGalDis,utilities.field_keys.keys())
        return
    
    def process(self,fieldKey):
        if fieldKey=='aegis':
            return
        fieldname   =   utilities.field_keys[fieldKey]
        self.log.info('Processing field: %s' %(fieldname))
        anaFieldDir =   os.path.join(self.config.anaDir,'%s' %(fieldname))
        
        self.log.info('Achieving the distriutions from observations')
        obsFname =  os.path.join(self.config.obsDir,'%s.fits' %(fieldname))
        obsOFname=  os.path.join(anaFieldDir,'galDisObs.csv')
        if not os.path.exists(obsOFname):
            t_dis_obs   =   self.getObsGalDis(obsFname)
            t_dis_obs.write(obsOFname)
        else:
            t_dis_obs   =   Table.read(obsOFname)
        
        self.log.info('Achieving the distriutions from simulations')
        simFname    =   os.path.join(anaFieldDir,'%s.fits' %(fieldname))
        simOFname   =   os.path.join(anaFieldDir,'galDisSim.csv')
        if not os.path.exists(simOFname):
            t_dis_sim.write(simOFname)
        else:
            t_dis_sim   =   Table.read(simOFname)
            
        self.plotGalDis(t_dis_obs,t_dis_sim,anaFieldDir)
        return
    
    def getObsGalDis(self,obsFname):
        data    =   fitsio.read(obsFname)
        mask    =   np.ones(len(data)).astype(bool)
        mask    =   mask&(~np.isnan(data['reGaus_resolution']))
        mask    =   mask&(~np.isnan(data['CModel_snr']))
        mask    =   mask&(~np.isnan(data['CModel_mag']))
        mask    =   mask&(~np.isnan(data['reGaus_sigmaE']))
        mask    =   mask&((data['CModel_mag'])<24.5)
        mask    =   mask&(data['extendness']>0)
        mask    =   mask&(data['ishape_hsm_regauss_sigma']>0.)
        mask    =   mask&(data['ishape_hsm_regauss_sigma']<0.4)
        mask    =   mask&(data['reGaus_resolution']>0.3)
        mask    =   mask&(data['CModel_snr']>10.)
        mask    =   mask&(data['reGaus_sigmaE']<2.)
        data    =   data[mask]
        absE    =   data['reGaus_sigmaE']
        Res     =   data['reGaus_resolution']
        mag     =   data['CModel_mag']
        snr     =   data['CModel_snr']
        ndata   =   len(data)*np.ones(len(self.nbin2-1))
        histRes =   np.histogram(Res,bins=self.resolution_bin,density=True)[0]
        histSnr =   np.histogram(snr,bins=self.snr_bin,density=True)[0]
        histMag =   np.histogram(mag,bins=self.mag_bin,density=True)[0]
        histAbsE=   np.histogram(absE,bins=self.absE_bin,density=True)[0]
        obscols = [self.resGrid,histRes,self.snrGrid,histSnr,self.magGrid,histMag,self.absEGrid,histAbsE,ndata]
        t_dis_obs=   Table(data=obscols,names=self.datanames)
        
        return t_dis_obs
        
    def getSimGalDis(self,cache,simFname):
        #Get the noise variance and psf information
        ifield  =   anafname.split('-')[1].split('.')[0]
        data    =   fitsio.read(simFname)
        #do mask
        data    =   data[data['weak_lensing_flag']]
        ndata   =   len(data)*np.ones(len(self.nbin2-1))
        Res     =   data['ext_shapeHSM_HsmShapeRegauss_resolution']
        snr     =   data['modelfit_CModel_flux']/data['modelfit_CModel_fluxSigma']
        absE    =   data['ext_shapeHSM_HsmShapeRegauss_e1']**2.+data['ext_shapeHSM_HsmShapeRegauss_e2']**2.
        absE    =   np.sqrt(absE)
        mag     =   data['Mag']
        weight  =   data['reweight_s16a']
        histRes =   np.histogram(Res,bins=self.resolution_bin,weights=weight,density=True)[0]
        histSnr =   np.histogram(snr,bins=self.snr_bin,weights=weight,density=True)[0]
        histMag =   np.histogram(mag,bins=self.mag_bin,weights=weight,density=True)[0]
        histAbsE=   np.histogram(absE,bins=self.absE_bin,weights=weight,density=True)[0]
        obscols = [self.resGrid,histRes,self.snrGrid,histSnr,self.magGrid,histMag,self.absEGrid,histAbsE,ndata]
        t_dis_sim=   Table(data=obscols,names=self.datanames)
        return t_dis_sim
    
    
    def plotGalDis(self,t_dis_obs,t_dis_rw,anaFieldDir):
        fig =   plt.figure()
        ax  =   fig.add_subplot(111)
        ax.plot(t_dis_obs['gridRes'],t_dis_obs['histRes'],'-',label='obs',c='b')
        ax.plot(t_dis_rw['gridRes'],t_dis_rw['histRes'],'--',label='reweight',c='r')

        ax.legend()
        plt.xlabel('resolution')
        plt.ylabel('prob')
        plt.savefig(os.path.join(anaFieldDir,'resHist.png'))

        fig =   plt.figure()
        ax  =   fig.add_subplot(111)
        ax.plot(t_dis_obs['gridSnr'],t_dis_obs['histSnr'],'-',label='obs',c='b')
        ax.plot(t_dis_rw['gridSnr'],t_dis_rw['histSnr'],'--',label='reweight',c='r')

        ax.legend()
        plt.xlabel('CModel SNR')
        plt.ylabel('prob')
        plt.savefig(os.path.join(anaFieldDir,'snrHist.png'))
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
