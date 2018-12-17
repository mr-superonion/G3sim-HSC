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

from lsst.pipe.base import TaskRunner
from lsst.ctrl.pool.parallel import BatchPoolTask
from lsst.ctrl.pool.pool import Pool, abortOnError


class rgcSimConfig(pexConfig.Config):
    "config"
    def setDefaults(self):
        pexConfig.Config.setDefaults(self)

class rgcSimTask(pipeBase.CmdLineTask):
    _DefaultName = "rgcSim"
    ConfigClass = rgcSimConfig
    def __init__(self,**kwargs):
        pipeBase.CmdLineTask.__init__(self, **kwargs)
    
    def SimulateExposure(self,ifield,outputdir):
        # Basic parameters
        ngrid       =   64
        nx          =   100 
        ny          =   nx
        ndata       =   nx*ny
        scale       =   0.168
        ngridTot    =   ngrid*nx
        bigfft      =   galsim.GSParams(maximum_fft_size=10240)
        gal_image   =   galsim.ImageF((nx+2)*ngrid,(ny+2)*ngrid,scale=scale)
        g1          =   0.02
        g2          =   0.

        # Get the noise image
        ud          =   galsim.UniformDeviate(ifield)
        varNoi      =   0.008769248618965289
        corData     =   pyfits.getdata('correlation.fits')
        corImg      =   galsim.ImageF(corData,scale=scale)
        #corNoise    =   galsim.CorrelatedNoise(corImg,ud,scale=scale,corImg=True)
        corNoise    =   galsim.GaussianNoise(ud)
        corNoise    =   corNoise.withVariance(varNoi)
        
        # Get the  galaxy image      
        flux_scaling=   2.587
        # Load data
        catName     =   'real_galaxy_catalog_25.2.fits'
        dir         =   '../galsim_train/COSMOS_25.2_training_sample/'
        cosmos_cat  =   galsim.COSMOSCatalog(catName, dir=dir)
        # index
        index_use   =   cosmos_cat.orig_index
        # parametric catalog
        param_cat   =   cosmos_cat.param_cat[index_use]
        index       =   ifield%81400
        record      =   param_cat[index]
        # prepare the galaxies        
        gal0        =   cosmos_cat.makeGalaxy(gal_type='parametric',index=index,gsparams=bigfft)
        gal0        *=  flux_scaling
        use_bulgefit=   record['use_bulgefit']
        if use_bulgefit:
            bparams = record['bulgefit']
            gal_q   = bparams[3]
            gal_beta= bparams[7]*galsim.radians
            hlr     = record['hlr'][2]
            bulge_flux= record['flux'][1]
            disk_flux= record['flux'][2]
            flux    = (bulge_flux+disk_flux)*flux_scaling
        else:
            sparams =   record['sersicfit']
            gal_q   =   sparams[3]
            gal_beta=   sparams[7]*galsim.radians
            hlr     =   record['hlr'][0]
            flux    =   record['flux'][0]*flux_scaling
        # Get the SNR
        fluxRadius  =   record['flux_radius']
        snr         =   flux/fluxRadius/np.sqrt(np.pi)/np.sqrt(varNoi)
        for i in range(ndata):
            # Prepare the subimage
            ix      =   i%nx
            iy      =   i/nx
            b       =   galsim.BoundsI(ix*ngrid+1, (ix+1)*ngrid,iy*ngrid+1,(iy+1)*ngrid)
            sub_gal_image = gal_image[b]
            # rotate the galaxy
            ud      =   galsim.UniformDeviate(ifield*1000+i)
            ang     =   ud()*2.*np.pi * galsim.radians
            gal     =   gal.rotate(ang)
            # Shear the galaxy
            gal     =   gal.shear(g1=g1,g2=g2)
            gal     =   galsim.Convolve([psf,gal],gsparams=bigfft)
            # Draw the galaxy image
            gal.drawImage(sub_gal_image)
        corNoise.applyTo(gal_image)
        exposure    =   afwImg.ExposureF((nx+2)*ngrid,(ny+2)*ngrid)
        exposure.getMaskedImage().getImage().getArray()[:,:]=gal_image.array
        exposure.getMaskedImage().getVariance().getArray()[:,:]=varNoi
        logging.info('finish making exposure')
        psfImg      =   afwImg.ImageF(45,45)
        psfImg.getArray()[:,:]= psfGalsim.array
        psfImg      =   psfImg.convertD()
        kernel      =   afwMath.FixedKernel(psfImg)
        kernelPSF   =   meaAlg.KernelPsf(kernel)
        exposure.setPsf(kernelPSF)
        return exposure

    
    @pipeBase.timeMethod
    def run(self,ifield):
        #Output 
        outputdir=  './expSim/' 
        if not os.path.exists(outputdir):
            os.popen('mkdir %s' %(outputdir) )
        logging.info('simulating field %s' %(ifield))
        #Get PSF
        self.SimulateExposure(ifield,outputdir)
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

class rgcSimBatchConfig(pexConfig.Config):
    minField =   pexConfig.Field(dtype=int, default=0, doc = 'minField')
    maxField =   pexConfig.Field(dtype=int, default=800, doc = 'minField')
    rgcSim = pexConfig.ConfigurableField(
        target = rgcSimTask,
        doc = "rgcSim task to run on multiple cores"
    )
    
class rgcSimRunner(TaskRunner):
    @staticmethod
    def getTargetList(parsedCmd, **kwargs):
        return [(ref, kwargs) for ref in range(1)] 

def unpickle(factory, args, kwargs):
    """Unpickle something by calling a factory"""
    return factory(*args, **kwargs)

class rgcSimBatchTask(BatchPoolTask):
    ConfigClass = rgcSimBatchConfig
    RunnerClass = rgcSimRunner
    _DefaultName = "rgcSimBatch"

    def __reduce__(self):
        """Pickler"""
        return unpickle, (self.__class__, [], dict(config=self.config, name=self._name,
                parentTask=self._parentTask, log=self.log))

    def __init__(self,**kwargs):
        BatchPoolTask.__init__(self, **kwargs)
        self.makeSubtask("rgcSim")
    
    @abortOnError
    def run(self,Id):
        minField=   self.config.minField
        maxField=   self.config.maxField
        fMin    =   minField
        fMax    =   maxField
        #Prepare the pool
        pool    =   Pool("rgcSim")
        pool.cacheClear()
        fieldList=  range(fMin,fMax)
        pool.map(self.process,fieldList)
        return
        
    def process(self,cache,ifield):
        return self.rgcSim.run(ifield)

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
