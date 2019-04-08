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
        return [(ref, kwargs) for ref in utilities.field_keys.keys()] 

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
        simFname    =   os.path.join(anaFieldDir,'galDisSim.csv')
        self.log.info('Achieving the distriutions from simulations')
        if not os.path.exists(simFname):
            anafnames   =   glob.glob('%s/*.fits' %(anaFieldDir)) 
            #Prepare the pool
            pool        =   Pool("galDisAnalysis")
            pool.cacheClear()
            outcomes    =   pool.map(self.getSimGalDis,anafnames)
            outcomes    =   [x for x in outcomes if x is not None]
            outcomes    =   np.stack(outcomes)
            histAll     =   np.sum(outcomes,axis=0)
            datacols    =   [self.resGrid,histAll[0],self.snrGrid,histAll[1],self.magGrid,histAll[2],self.absEGrid,histAll[3]]
            datanames   =   ['gridRes','histRes','gridSnr','histSnr','gridMag','histMag','gridAbsE','histAbsE']
            t_dis_obs   =   Table(data=datacols,names=datanames)
            t_dis_obs.write(simFname)
        else:
            t_dis_obs   =   Table.read(simFname)
        
        self.log.info('Achieving the distriutions from observations')
        obsFname =   os.path.join(self.config.obsDir,'%s.fits' %(fieldname))
        """
        t_dis_obs=  self.getObsHist(obsFname)
        t_dis_obs.write(outDnameObs,overwrite=True)
            
        outDnameRw  =   os.path.join(anaFieldDir,'varPSF_dis_rw.csv')
        t_dis_rw    =   self.reWeight(t,t_dis,t_dis_obs)
        t_dis_rw.write(outDnameRw,overwrite=True)
        t.write(outFname,overwrite=True)
        self.plotNoiPsf(t_dis,t_dis_obs,t_dis_rw,anaFieldDir)
        """
        return
    
    def getSimGalDis(self,cache,anafname):
        #Get the noise variance and psf information
        ifield  =   anafname.split('-')[1].split('.')[0]
        data    =   pyfits.getdata(anafname)
        #do mask
        data    =   data[data['weak_lensing_flag']]
        ndata   =   len(data)
        if ndata==0:
            return None
        Res     =   data['ext_shapeHSM_HsmShapeRegauss_resolution']
        snr     =   data['modelfit_CModel_flux']/data['modelfit_CModel_fluxSigma']
        absE    =   data['ext_shapeHSM_HsmShapeRegauss_e1']**2.+data['ext_shapeHSM_HsmShapeRegauss_e2']**2.
        absE    =   np.sqrt(absE)
        mag     =   data['Mag']
        rw      =   data[0]['reweight_s16a']
        histRes =   np.histogram(Res,bins=self.resolution_bin,density=True)[0]*ndata*rw
        histSnr =   np.histogram(snr,bins=self.snr_bin,density=True)[0]*ndata*rw
        histMag =   np.histogram(mag,bins=self.mag_bin,density=True)[0]*ndata*rw
        histAbsE=   np.histogram(absE,bins=self.absE_bin,density=True)[0]*ndata*rw
        return  np.vstack([histRes,histSnr,histMag,histAbsE])
    
    
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
