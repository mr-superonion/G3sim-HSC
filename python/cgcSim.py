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
import galsim
import numpy as np
from astropy.table import Table
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


class cgcSimConfig(pexConfig.Config):
    preDir      =   's19-pre'
    outDir      =   's19-1'
    expDir      =   pexConfig.Field(dtype=str, default='%s/expDir' %outDir, doc = 'directory to store exposures')
    catPre      =   pexConfig.Field(dtype=str, default='%s/catPre' %outDir, doc = 'directory to store pre-catalog')
    psfPre      =   pexConfig.Field(dtype=str, default='%s/psfPre' %preDir, doc = 'directory to store pre-PSF')
    corPre      =   pexConfig.Field(dtype=str, default='%s/corPre' %preDir, doc = 'directory to store pre-PSF')
    def setDefaults(self):
        pexConfig.Config.setDefaults(self)
    def validate(self):
        pexConfig.Config.validate(self)
        if not os.path.exists(self.expDir):
            os.mkdir(self.expDir)
        if not os.path.exists(self.catPre):
            os.mkdir(self.catPre)
        if not os.path.exists(self.psfPre):
            raise RuntimeError("Cannot find PSF images")
        if not os.path.exists(self.corPre):
            raise RuntimeError("Cannot find noise correlations images")

class cgcSimTask(pipeBase.CmdLineTask):
    _DefaultName = "cgcSim"
    ConfigClass = cgcSimConfig
    def __init__(self,**kwargs):
        pipeBase.CmdLineTask.__init__(self, **kwargs)

    @pipeBase.timeMethod
    def runDataRef(self,ihpix):
        self.log.info('begining for healPIX %05d' %(ihpix))
        outFname    =   os.path.join(self.config.expDir,'image-%05d.fits' %(ihpix))
        if os.path.exists(outFname):
            self.log.info('Already have the outcome')
            return
        self.log.info('simulating field %s' %(ihpix))
        # Basic parameters
        scale       =   0.168
        bigfft      =   galsim.GSParams(maximum_fft_size=10240)
        flux_scaling=   2.587
        varRatio    =   1.25

        # Get the shear information
        g1          =   0.02
        g2          =   0.
        fwhm        =   noiseHead['fwhm']*0.168
        variance    =   noiseHead['varPlane']/varRatio#1.78

        # Get the PSF
        psfFname    =   os.path.join(self.config.psfPre,'psf-%05d.fits'%(ifield))
        psfImg      =   galsim.fits.read(psfFname)
        noiseHead   =   pyfits.getheader(psfFname)
        psfInt      =   galsim.InterpolatedImage(psfImg,scale=scale,flux = 1.)

        self.log.info('Adding correlated noise')
        """
        Simulate the original noise on image uniform
        """
        rng = galsim.BaseDeviate(ihpix)
        max_variance=   np.max(var_image.array)
        var_image   =   max_variance - var_image
        vn          =   galsim.VariableGaussianNoise(rng,var_image)
        gal_image.addNoise(vn)

        """
        Make noise correlated and match to the input variance
        """
        corFname    =   os.path.join(self.config.corPre,'corGroup-%02d.fits' %group)
        corNoise    =   galsim.getCOSMOSNoise(file_name=corFname,rng=rng,cosmos_scale=scale,variance=variance)
        corNoise.applyTo(gal_image)

        """
        Generate the image
        """
        exposure    =   afwImg.ExposureF(nx*ngrid,ny*ngrid)
        exposure.getMaskedImage().getImage().getArray()[:,:]=gal_image.array
        exposure.getMaskedImage().getVariance().getArray()[:,:]=variance*varRatio#1.78
        del gal_image
        del var_image
        #Set the PSF
        psfArray    =   psfImg.array
        ngridPsf    =   psfArray.shape[0]
        psfLsst     =   afwImg.ImageF(ngridPsf,ngridPsf)
        psfLsst.getArray()[:,:]= psfArray
        psfLsst     =   psfLsst.convertD()
        kernel      =   afwMath.FixedKernel(psfLsst)
        kernelPSF   =   meaAlg.KernelPsf(kernel)
        exposure.setPsf(kernelPSF)
        #prepare the wcs
        #Rotation
        cdelt   =   (0.168*afwGeom.arcseconds)
        CD      =   afwGeom.makeCdMatrix(cdelt, afwGeom.Angle(0.))#no rotation
        #wcs
        crval   =   afwGeom.SpherePoint(afwGeom.Angle(0.,afwGeom.degrees),afwGeom.Angle(0.,afwGeom.degrees))
        crpix   =   afwGeom.Point2D(0.0, 0.0)
        dataWcs =   afwGeom.makeSkyWcs(crpix,crval,CD)
        exposure.setWcs(dataWcs)

        #prepare the frc
        dataCalib = afwImg.makePhotoCalibFromCalibZeroPoint(63095734448.0194)
        exposure.setPhotoCalib(dataCalib)
        self.log.info('writing exposure')
        exposure.writeFits(outFname)
        del exposure
        t = Table(rows=data_rows, names=('ipos','COSMOS_INDEX','IDENT','RA','DEC','COSMOS_MAG','a_i','parent_seeing'))
        catPreName  =   os.path.join(self.config.catPre,'catalog-%05d.fits' %(ihpix))
        t.write(catPreName,overwrite=True)
        del t
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

class cgcSimBatchConfig(pexConfig.Config):
    perGroup =   pexConfig.Field(dtype=int, default=100, doc = 'data per field')
    cgcSim = pexConfig.ConfigurableField(
        target = cgcSimTask,
        doc = "cgcSim task to run on multiple cores"
    )
class cgcSimRunner(TaskRunner):
    @staticmethod
    def getTargetList(parsedCmd, **kwargs):
        minGroup    =  parsedCmd.minGroup
        maxGroup    =  parsedCmd.maxGroup
        return [(ref, kwargs) for ref in range(minGroup,maxGroup)]
def unpickle(factory, args, kwargs):
    """Unpickle something by calling a factory"""
    return factory(*args, **kwargs)

class cgcSimBatchTask(BatchPoolTask):
    ConfigClass = cgcSimBatchConfig
    RunnerClass = cgcSimRunner
    _DefaultName = "cgcSimBatch"
    def __reduce__(self):
        """Pickler"""
        return unpickle, (self.__class__, [], dict(config=self.config, name=self._name,
                parentTask=self._parentTask, log=self.log))
    def __init__(self,**kwargs):
        BatchPoolTask.__init__(self, **kwargs)
        self.makeSubtask("cgcSim")
        return

    @abortOnError
    def runDataRef(self,Id):
        self.log.info('beginning group %d' %(Id))
        perGroup=   self.config.perGroup
        fMin    =   perGroup*Id
        fMax    =   perGroup*(Id+1)
        #Prepare the pool
        pool    =   Pool("cgcSim")
        pool.cacheClear()
        fieldList=  range(fMin,fMax)
        pool.map(self.process,fieldList)
        self.log.info('finish group %d'%(Id) )
        return

    def process(self,cache,isim):
        """
        try:
        except:
            self.log.info('fail to finish field %05d'%(isim) )
        """
        self.cgcSim.runDataRef(isim)
        self.log.info('finish field %05d' %(isim))
        return

    @classmethod
    def _makeArgumentParser(cls, *args, **kwargs):
        kwargs.pop("doBatch", False)
        parser = pipeBase.ArgumentParser(name=cls._DefaultName)
        parser.add_argument('--minGroup', type= int,
                        default=0,
                        help='minimum group number')
        parser.add_argument('--maxGroup', type= int,
                        default=10,
                        help='maximum group number')
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
