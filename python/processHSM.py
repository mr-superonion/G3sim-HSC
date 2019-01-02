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
#
# python lib
import os
import logging
import galsim
import numpy as np
import astropy.io.fits as pyfits

# lsst Tasks
import lsst.daf.base as dafBase
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
import lsst.afw.math as afwMath
import lsst.afw.table as afwTable
import lsst.afw.image as afwImg
import lsst.afw.detection as afwDet
import lsst.afw.geom as afwGeom
import lsst.afw.coord as afwCoord
import lsst.meas.algorithms as meaAlg

from lsst.pipe.base import TaskRunner
from lsst.ctrl.pool.parallel import BatchPoolTask
from lsst.ctrl.pool.pool import Pool, abortOnError


class processHSMConfig(pexConfig.Config):
    "config"
    detection = pexConfig.ConfigurableField(
        target = SourceDetectionTask,
        doc = "Detect sources"
    )
    deblend = pexConfig.ConfigurableField(
        target = SourceDeblendTask,
        doc = "Split blended source into their components"
    )
    measurement = pexConfig.ConfigurableField(
        target = SingleFrameMeasurementTask,
        doc = "Measure sources"
    )
    catalogCalculation = pexConfig.ConfigurableField(
        target = CatalogCalculationTask,
        doc = "Subtask to run catalogCalculation plugins on catalog"
    )
    def setDefaults(self):
        pexConfig.Config.setDefaults(self)
        self.detection.thresholdValue = 5.0
        self.detection.isotropicGrow  = True
        self.detection.reEstimateBackground=False
        self.deblend.propagateAllPeaks  =   True 
        self.deblend.maxFootprintArea   =   64*64
        self.deblend.maxFootprintSize   =   64
        self.measurement.load(os.path.join(getPackageDir("obs_subaru"), "config", "hsm.py"))
        self.measurement.slots.apFlux       =   None
        self.measurement.slots.instFlux     =   None
        self.measurement.slots.calibFlux    =   None

class processHSMTask(pipeBase.CmdLineTask):
    _DefaultName = "processHSM"
    ConfigClass = processHSMConfig
    def __init__(self,**kwargs):
        pipeBase.CmdLineTask.__init__(self, **kwargs)
        schema     =   afwTable.SourceTable.makeMinimalSchema()
        self.schemaMapper = afwTable.SchemaMapper(schema)
        self.schemaMapper.addMinimalSchema(schema)
        self.schema= self.schemaMapper.getOutputSchema()
        self.algMetadata = dafBase.PropertyList()
        self.makeSubtask('detection', schema=self.schema)
        self.makeSubtask("deblend", schema=self.schema)
        self.makeSubtask('measurement', schema=self.schema, algMetadata=self.algMetadata)
        self.makeSubtask('catalogCalculation', schema=self.schema)
        
    @pipeBase.timeMethod
    def run(self,ifield,rgc):
        self.log.info('begining for field %04d' %(ifield))
        outputdir   =  './outcome/' %(sName)
        outFname    =   'measHSM%04d.fits'  %(ifield)
        outFname    =   os.path.join(outputdir,outFname)
        table = afwTable.SourceTable.make(self.schema)
        sources = afwTable.SourceCatalog(table)
        table.setMetadata(self.algMetadata)
        detRes  = self.detection.run(table=table, exposure=exposure, doSmooth=True)
        sources = detRes.sources
        # do deblending
        self.deblend.run(exposure=exposure,sources=sources)
        # do measurement
        self.measurement.run(measCat=sources,exposure=exposure)
        # measurement on the catalog level
        self.catalogCalculation.run(sources)
        sources.writeFits(outFname)
        return
    
    
    @classmethod
    def _makeArgumentParser(cls):
        """Create an argument parser
        """
        doBatch = kwargs.pop("doBatch", False)
        parser = pipeBase.ArgumentParser(name=cls._DefaultName)
        return parser

    def writeConfig(self, butler, clobber=False, doBackup=False):
        pass

    def writeSchemas(self, butler, clobber=False, doBackup=False):
        pass

    def writeMetadata(self, dataRef):
        pass

    def writeEupsVersions(self, butler, clobber=False, doBackup=False):
        pass

class processHSMBatchConfig(pexConfig.Config):
    minField =   pexConfig.Field(dtype=int, default=0, doc = 'minField')
    maxField =   pexConfig.Field(dtype=int, default=800, doc = 'minField')
    processHSM = pexConfig.ConfigurableField(
        target = processHSMTask,
        doc = "processHSM task to run on multiple cores"
    )
    
class processHSMRunner(TaskRunner):
    @staticmethod
    def getTargetList(parsedCmd, **kwargs):
        return [(ref, kwargs) for ref in range(1)] 

def unpickle(factory, args, kwargs):
    """Unpickle something by calling a factory"""
    return factory(*args, **kwargs)

class processHSMBatchTask(BatchPoolTask):
    ConfigClass = processHSMBatchConfig
    RunnerClass = processHSMRunner
    _DefaultName = "processHSMBatch"

    def __reduce__(self):
        """Pickler"""
        return unpickle, (self.__class__, [], dict(config=self.config, name=self._name,
                parentTask=self._parentTask, log=self.log))

    def __init__(self,**kwargs):
        BatchPoolTask.__init__(self, **kwargs)
        self.makeSubtask("processHSM")
    
    @abortOnError
    def run(self,Id):
        # Get the  galaxy generator      
        # Load data
        Xnames      =   ['best','median','worst']
        rgcs        =   []
        for iXn in range(3):
            Xname       =   Xnames[iXn]
            catName     =   'real_galaxy_catalog_%s.fits' %(Xname)
            directory   =   '/work/xiangchong.li/work/galsim_train/parent_%s_processed' %(Xname)
            rgcs.append(galsim.RealGalaxyCatalog(catName,dir=directory))
        minField=   self.config.minField
        maxField=   self.config.maxField
        fMin    =   minField
        fMax    =   maxField
        #Prepare the pool
        pool    =   Pool("processHSM")
        pool.cacheClear()
        pool.storeSet(rgcs=rgcs)
        fieldList=  range(fMin,fMax)
        pool.map(self.process,fieldList)
        return
        
    def process(self,cache,ifield):
        rgc =   cache.rgcs[ifield%3]
        try:
            self.processHSM.run(ifield,rgc)
            self.log.info('finish field %04d' %(ifield))
        except:
            self.log.info('fail to finish field %04d' %(ifield))
        return 

    @classmethod
    def _makeArgumentParser(cls, *args, **kwargs):
        kwargs.pop("doBatch", False)
        parser = pipeBase.ArgumentParser(name=cls._DefaultName)
        return parser
    
    @classmethod
    def writeConfig(self, butler, clobber=False, doBackup=False):
        pass

    def writeSchemas(self, butler, clobber=False, doBackup=False):
        pass
    
    def writeMetadata(self, dataRef):
        pass
    
    def writeEupsVersions(self, butler, clobber=False, doBackup=False):
        pass
    
    def _getConfigName(self):
        return None
   
    def _getEupsVersionsName(self):
        return None
    
    def _getMetadataName(self):
        return None
