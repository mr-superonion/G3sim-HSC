#!/usr/bin/env python

import os
import sys
import optparse
import sim_analysis_utilities as utilities
import numpy as np
from reGaus_sigmaE import reGaus_sigmaE_Task
from astropy.table import Table, vstack

import lsst.daf.base as dafBase
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase

from lsst.pipe.base import TaskRunner
from lsst.ctrl.pool.parallel import BatchPoolTask
from lsst.ctrl.pool.pool import Pool, abortOnError


class reGausAnalysisConfig(pexConfig.Config):
    anaDir  =   pexConfig.Field(dtype=str, default='anaCat', doc = 'simulation catalogs')
    obsDir  =   pexConfig.Field(dtype=str, default='./', doc = 'observed catalogs')
    pltDir  =   pexConfig.Field(dtype=str, default='plot', doc = 'plot directory')
    pltType =   pexConfig.Field(dtype=str, default='eps', doc = 'plot directory')
    nSims   =   pexConfig.Field(dtype=int, default=3200, doc = 'number of simulations')
    rescale_noise  =   pexConfig.Field(dtype=bool, default=True, doc = 'rescale measured noise to target noise?')
    useCalib=   pexConfig.Field(dtype=bool, default=False, doc = 'apply calibration?')
    reGaus_sigmaE= pexConfig.ConfigurableField(
        target  = reGaus_sigmaE_Task,
        doc     = "Subtask to determine sigmae"
    )
    
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

class reGausAnalysisRunner(TaskRunner):
    @staticmethod
    def getTargetList(parsedCmd, **kwargs):
        return [(ref, kwargs) for ref in range(1)] 

def unpickle(factory, args, kwargs):
    """Unpickle something by calling a factory"""
    return factory(*args, **kwargs)

class reGausAnalysisTask(BatchPoolTask):
    ConfigClass     =   reGausAnalysisConfig
    RunnerClass     =   reGausAnalysisRunner
    _DefaultName    =   "reGausAnalysis"
    
    def __init__(self,**kwargs):
        BatchPoolTask.__init__(self, **kwargs)
        self.makeSubtask('reGaus_sigmaE')
    
    def __reduce__(self):
        """Pickler"""
        return unpickle, (self.__class__, [], dict(config=self.config, name=self._name,
                parentTask=self._parentTask, log=self.log))
    
    @abortOnError
    def run(self,Id):
        #Prepare the pool
        rows    =   Table.read('noiVar_rPsf_s18a.csv')
        pool    =   Pool("reGausAnalysis")
        pool.cacheClear()
        outcomes=   pool.map(self.process,rows)
        outtable=   vstack(outcomes)
        self.reGaus_sigmaE.finalise(outtable,self.config.pltDir,self.config.pltType)
        return
        
    def process(self,cache,row):
        ifield  =   row['num']
        self.log.info('Processing subfield %s' %(ifield))
        fieldname=  utilities.field_keys[row['field']]
        fname   =   os.path.join(self.config.anaDir,fieldname,'src-%04d.fits'%(ifield))
        try:
            catalog =   Table.read(fname)   
        except:
            self.log.info(fname)
            return
        return self.reGaus_sigmaE.run(catalog)

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
