#!/usr/bin/env python
import os
import smatch
import fitsio
import numpy as np
import catSelect
from astropy.table import Table, vstack
import sim_analysis_utilities as utilities

import matplotlib
matplotlib.use('agg')
matplotlib.rcParams['ps.useafm'] = True
matplotlib.rcParams['pdf.use14corefonts'] = True
matplotlib.rcParams['text.usetex'] = True

import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import matplotlib.lines as mlines
from matplotlib.colors import LogNorm
import lsst.daf.base as dafBase
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase

from lsst.pipe.base import TaskRunner
from lsst.ctrl.pool.parallel import BatchPoolTask
from lsst.ctrl.pool.pool import Pool, abortOnError

class galDisAnalysis2DSizeConfig(pexConfig.Config):
    anaDir  =   pexConfig.Field(dtype=str, default='anaCat', doc = 'simulation catalogs directory')
    obsDir  =   pexConfig.Field(dtype=str, default='obsCat', doc = 'observed catalogs directory')
    pltDir  =   pexConfig.Field(dtype=str, default='plot/galDis/', doc = 'plot directory')
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
            os.mkdir(self.pltDir)

class galDisAnalysis2DSizeRunner(TaskRunner):
    @staticmethod
    def getTargetList(parsedCmd, **kwargs):
        return [(ref, kwargs) for ref in range(1)]
def unpickle(factory, args, kwargs):
    """Unpickle something by calling a factory"""
    return factory(*args, **kwargs)
class galDisAnalysis2DSizeTask(BatchPoolTask):
    ConfigClass     =   galDisAnalysis2DSizeConfig
    RunnerClass     =   galDisAnalysis2DSizeRunner
    _DefaultName    =   "galDisAnalysis2DSize"
    def __init__(self,**kwargs):
        BatchPoolTask.__init__(self, **kwargs)
        self.datanames   =   ['res','snr','mag','traceR','logb','weight']
    def __reduce__(self):
        """Pickler"""
        return unpickle, (self.__class__, [], dict(config=self.config, name=self._name,
                parentTask=self._parentTask, log=self.log))

    @abortOnError
    def run(self,Id):
        if not os.path.exists(os.path.join(self.config.anaDir,'allObs4.fits')):
            self.log.info('canot find input files')
            return
        else:
            self.log.info('already have prepared catalog')
            allSim=Table.read(os.path.join(self.config.anaDir,'allSim4.fits'))
            allObs=Table.read(os.path.join(self.config.anaDir,'allObs4.fits'))
        #self.finalPlot_gal(allSim,allObs)
        self.finalPlot_psf2(allSim,allObs)
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
