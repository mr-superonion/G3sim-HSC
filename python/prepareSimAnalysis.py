#!/usr/bin/env python
#
# LSST Data Management System
# Copyright 20082014 LSST Corporation.
#
# This product includes software developed by the
# LSST Project (http://www.lsst.org/).
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.    See the
# GNU General Public License for more details.
#
# You should have received a copy of the LSST License Statement and
# the GNU General Public License along with this program.  If not,
# see <http://www.lsstcorp.org/LegalNotices/>.
"""
 This script converts the HSC fits data files into one HDF-5 file for faster
 access where only columns of interest can be selectively loaded
 Plus, it now also aggregate the noise rescaling factors for each (tract, patch)
"""

import sys
import os
from astropy.table import Table, vstack
import sim_analysis_utilities as utilities

import lsst.daf.base as dafBase
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase

from lsst.pipe.base import TaskRunner
from lsst.ctrl.pool.parallel import BatchPoolTask
from lsst.ctrl.pool.pool import Pool, abortOnError

class prepareSimAnalysisConfig(pexConfig.Config):
    anaDir  =   pexConfig.Field(dtype=str, default='anaCat', doc = 'output catalogs')
    simDir  =   pexConfig.Field(dtype=str, default='outcomeFPFS', doc = 'input catalogs')
    def setDefaults(self):
        pexConfig.Config.setDefaults(self)

    def validate(self):
        pexConfig.Config.validate(self)
        if not os.path.exists(self.anaDir):
            os.mkdir(self.anaDir)
        for key in utilities.field_keys:
            if key=='aegis':
                continue
            fN  =   utilities.field_keys[key]
            tmN =   os.path.join(self.anaDir,fN)
            if not os.path.exists(tmN):
                os.mkdir(tmN)

class processSimRunner(TaskRunner):
    @staticmethod
    def getTargetList(parsedCmd, **kwargs):
        return [(ref, kwargs) for ref in utilities.field_keys.keys()] 

def unpickle(factory, args, kwargs):
    """Unpickle something by calling a factory"""
    return factory(*args, **kwargs)

class prepareSimAnalysisTask(BatchPoolTask):
    ConfigClass     =   prepareSimAnalysisConfig
    RunnerClass     =   processSimRunner
    _DefaultName    =   "prepareSimAnalysis"
    
    def __init__(self,**kwargs):
        BatchPoolTask.__init__(self, **kwargs)
    
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
        #Prepare the pool
        pool    =   Pool("processSim")
        pool.cacheClear()
        simList =   Table.read('psfPre/info.csv')
        simList =   simList[simList['field']==fieldKey]['num']
        outcomes    =   pool.map(self.process,simList)
        outcomes    =   [x for x in outcomes if x is not None]
        catalogAll  =   vstack(outcomes)
        outFname    =   os.path.join(self.config.anaDir,fieldname,'%s.fits'%(fieldname))
        catalogAll.write(outFname)
        return
        
    def process(self,cache,ifield):
        self.log.info('Processing subfield %s' %(ifield))
        # Load the catalog from fits file
        filename    =   os.path.join(self.config.simDir,'src-%04d.fits' % ifield)
        # If the file does not exist, an IOError will be raised
        catalog     =   utilities.prepare_catalog(filename,isim=ifield)
        if catalog is None:
            self.log.info('failed to process subfield %d' %(ifield))
        return catalog 

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
