#!/usr/bin/env python
# This script converts the HSC fits data files into one HDF-5 file for faster
# access where only columns of interest can be selectively loaded
# Plus, it now also aggregate the noise rescaling factors for each (tract, patch)
import sys
import os
import optparse
from astropy.table import Table, vstack
from sim_analysis_utilities import *


class prepareSimAnalysisConfig(pexConfig.Config):
    anaDir  =   pexConfig.Field(dtype=str, default='anaCat', doc = 'output catalogs')
    simDir  =   pexConfig.Field(dtype=str, default='outcomeFPFS', doc = 'input catalogs')
    version =   pexConfig.Field(dtype=str, default='v1', doc = 'input catalogs')
    nSims   =   pexConfig.Field(dtype=int, default=3200, doc = 'number of simulations')
    rescale_noise  =   pexConfig.Field(dtype=bool, default=False, doc = 'output catalogs')
    def setDefaults(self):
        pexConfig.Config.setDefaults(self)

class processSimRunner(TaskRunner):
    @staticmethod
    def getTargetList(parsedCmd, **kwargs):
        return [(ref, kwargs) for ref in range(1)] 

def unpickle(factory, args, kwargs):
    """Unpickle something by calling a factory"""
    return factory(*args, **kwargs)

class prepareSimAnalysisTask(BatchPoolTask):
    ConfigClass     =   prepareSimAnalysisConfig
    RunnerClass     =   processSimRunner
    _DefaultName    =   "prepareSimAnalysis"
    
    def __reduce__(self):
        """Pickler"""
        return unpickle, (self.__class__, [], dict(config=self.config, name=self._name,
                parentTask=self._parentTask, log=self.log))

    def __init__(self,**kwargs):
        BatchPoolTask.__init__(self, **kwargs)
    
    @abortOnError
    def run(self,Id):
        # Get the  galaxy generator      
        # Load data
        #Prepare the pool
        pool    =   Pool("processSim")
        pool.cacheClear()
        fieldList=  range(self.config.nSims)
        pool.map(self.process,fieldList)
        return
        
    def process(self,cache,ifield):
         try:
            catalog  = prepare_sim_catalog(self.config.simDir,ifield,rescale_est_variance=self.config.rescale_noise)
        except IOError:
            # Technically there could be some other IOError here,
            # e.g., missing variance file rather than this subfield
            # being missing.  But for now let's keep it like this.
            self.log.info("Subfield %d missing, skipping to the next one."%(ifield))
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


if __name__ == "__main__":
    prepareSimAnalysis.parseAndSubmit()
