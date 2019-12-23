#!/usr/bin/env python
import os
import sys
import glob
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


class addMeasobsConfig(pexConfig.Config):
    obsDir  =   pexConfig.Field(dtype=str, default='', doc = 'observed catalogs directory')
    def setDefaults(self):
        pexConfig.Config.setDefaults(self)
    def validate(self):
        pexConfig.Config.validate(self)
        # Checks that the paths provided are valid
        if not os.path.isdir(self.obsDir):
            raise RuntimeError("Data catalog directory %s does not exist or is not a directory" % (self.anaDir))

class addMeasobsRunner(TaskRunner):
    @staticmethod
    def getTargetList(parsedCmd, **kwargs):
        #only 1 group
        return [(ref, kwargs) for ref in range(1)]
def unpickle(factory, args, kwargs):
    """Unpickle something by calling a factory"""
    return factory(*args, **kwargs)

class addMeasobsTask(BatchPoolTask):
    ConfigClass     =   addMeasobsConfig
    RunnerClass     =   addMeasobsRunner
    _DefaultName    =   "addMeasobs"

    def __init__(self,**kwargs):
        BatchPoolTask.__init__(self, **kwargs)
    def __reduce__(self):
        """Pickler"""
        return unpickle, (self.__class__, [], dict(config=self.config, name=self._name,
                parentTask=self._parentTask, log=self.log))

    @abortOnError
    def run(self,fname):
        #Prepare the pool
        pool        =   Pool("addMeasobs")
        pool.cacheClear()
        pool.map(self.process,tractList)
        return

    def process(self,cache,fname):
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
