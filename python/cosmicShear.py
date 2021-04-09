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
import gc
import fitsio
import catutil
import treecorr
import numpy as np

# lsst Tasks
import lsst.daf.base as dafBase
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase

from lsst.pipe.base import TaskRunner
from lsst.ctrl.pool.parallel import BatchPoolTask
from lsst.ctrl.pool.pool import Pool, abortOnError

import lsst.obs.subaru.filterFraction

class cosmicShearConfig(pexConfig.Config):
    def setDefaults(self):
        pexConfig.Config.setDefaults(self)
class cosmicShearRunner(TaskRunner):
    @staticmethod
    def getTargetList(parsedCmd, **kwargs):
        minGroup    =  parsedCmd.minGroup
        maxGroup    =  parsedCmd.maxGroup
        return [(ref, kwargs) for ref in range(minGroup,maxGroup)]

def unpickle(factory, args, kwargs):
    """Unpickle something by calling a factory"""
    return factory(*args, **kwargs)
class cosmicShearTask(BatchPoolTask):
    ConfigClass = cosmicShearConfig
    RunnerClass = cosmicShearRunner
    _DefaultName = "cosmicShear"
    def __reduce__(self):
        """Pickler"""
        return unpickle, (self.__class__, [], dict(config=self.config, name=self._name,
                parentTask=self._parentTask, log=self.log))
    def __init__(self,**kwargs):
        BatchPoolTask.__init__(self, **kwargs)

    @abortOnError
    def runDataRef(self,isim):
        outDir  =   'anaCat_newS19Mask_fdeltacut/cosmicShear/covariance/sim-%03d' %isim
        if not os.path.isdir(outDir):
            os.mkdir(outDir)
        self.log.info('beginning isim: %03d' %(isim))
        # 108 sims x 13 rotations
        rotList =  range(7)
        # Prepare the pool
        pool    =   Pool("cosmicShear")
        pool.cacheClear()
        pool.storeSet(isim=isim)
        pool.storeSet(outDir=outDir)
        pool.map(self.process,rotList)
        return

    def prepareData(self,isim,irot,outDir):
        nameList=   ['ra_mock','dec_mock','noise1_int','noise2_int','noise1_mea','noise2_mea',\
                     'shear1_sim','shear2_sim','e1_mock','e2_mock','weight']
        wDir    =   'anaCat_newS19Mask_fdeltacut/catalog_compressed/'
        fname   =   'mock_nres13_r%03d_rotmat%d_shear_catalog.fits' %(isim,irot)
        dataG   =   []
        for fieldname in catutil.field_names:
            # selection
            mask    =   fitsio.read(os.path.join(wDir,'%s_zmask.fits' %fieldname))['mask']
            gc.collect()
            catFname=   'anaCat_newS19Mask_fdeltacut/catalog_mock/shape/%s/%s' %(fieldname,fname)
            data    =   fitsio.read(os.path.join(catFname))[nameList]
            mask    =   mask&(data['e1_mock']**2.+data['e2_mock']**2.<10.)
            data    =   data[mask]
            dataG.append(data)
            del data,mask,catFname
            gc.collect()
        dataG=  np.hstack(dataG)
        return dataG

    def process(self,cache,irot):
        isim    =   cache.isim
        outDir  =   cache.outDir
        dataG   =   self.prepareData(isim,irot,outDir)
        if os.path.isfile(os.path.join(outDir,'rotmat%d_corG.csv' %irot)):
            return

        # correlations
        catE1   =   treecorr.Catalog(g1=dataG['e1_mock'],g2=-dataG['e2_mock'],w=dataG['weight'],\
                ra=dataG['ra_mock'],dec=dataG['dec_mock'],ra_units='deg',dec_units='deg')
        corE1   =   treecorr.GGCorrelation(nbins=30,min_sep=0.281838293,max_sep=354.8133892,sep_units='arcmin')
        corE1.process(catE1)
        corE1.write(os.path.join(outDir,'rotmat%d_corE1.csv' %irot))
        del corE1,catE1
        gc.collect()

        eI1     =   dataG['noise1_int']+dataG['noise1_mea']
        eI2     =   dataG['noise2_int']+dataG['noise2_mea']
        catE0   =   treecorr.Catalog(g1=eI1,g2=-eI2,w=dataG['weight'],\
                ra=dataG['ra_mock'],dec=dataG['dec_mock'],ra_units='deg',dec_units='deg')
        del eI1,eI2
        gc.collect()

        corE0   =   treecorr.GGCorrelation(nbins=30,min_sep=0.281838293,max_sep=354.8133892,sep_units='arcmin')
        corE0.process(catE0)
        corE0.write(os.path.join(outDir,'rotmat%d_corE0.csv' %irot))
        del corE0
        gc.collect()

        catG=   treecorr.Catalog(g1=dataG['shear1_sim'],g2=-dataG['shear2_sim'],w=dataG['weight'],\
                ra=dataG['ra_mock'],dec=dataG['dec_mock'],ra_units='deg',dec_units='deg')
        del dataG
        gc.collect()

        corGE0  =   treecorr.GGCorrelation(nbins=30,min_sep=0.281838293,max_sep=354.8133892,sep_units='arcmin')
        corGE0.process(catG,catE0)
        corGE0.write(os.path.join(outDir,'rotmat%d_corGE0.csv' %irot))
        del corGE0,catE0
        gc.collect()

        corG    =   treecorr.GGCorrelation(nbins=30,min_sep=0.281838293,max_sep=354.8133892,sep_units='arcmin')
        corG.process(catG)
        corG.write(os.path.join(outDir,'rotmat%d_corG.csv' %irot))
        del corG,catG
        gc.collect()
        self.log.info('finish irot: %02d' %(irot))
        return

    @classmethod
    def _makeArgumentParser(cls, *args, **kwargs):
        kwargs.pop("doBatch", False)
        parser = pipeBase.ArgumentParser(name=cls._DefaultName)
        parser.add_argument('--minGroup', type= int,
                        default=0,
                        help='minimum group number')
        parser.add_argument('--maxGroup', type= int,
                        default=9,
                        help='maximum group number')
        return parser
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
