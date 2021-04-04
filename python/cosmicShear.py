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
        rotList =  range(13)
        # Prepare the pool
        pool    =   Pool("cosmicShear")
        pool.cacheClear()
        pool.storeSet(isim=isim)
        pool.storeSet(outDir=outDir)
        pool.map(self.process,rotList)
        return

    def process(self,cache,irot):
        isim    =   cache.isim
        outDir  =   cache.outDir
        fname   =   'mock_nres13_r%03d_rotmat%d_shear_catalog.fits' %(isim,irot)
        wDir    =   'anaCat_newS19Mask_fdeltacut/catalog_compressed/'
        fname   =   'mock_nres13_r003_rotmat12_shear_catalog.fits'
        dataG   =   []
        wG      =   []
        for fieldname in catutil.field_names:
            wFname=os.path.join(wDir,'%s_weights.fits' %fieldname)
            weights=   fitsio.read(wFname)
            catFname=   'anaCat_newS19Mask_fdeltacut/catalog_mock/shape/%s/%s' %(fieldname,fname)
            data   =   fitsio.read(os.path.join(catFname))
            assert np.all(data['object_id']==weights['object_id'])
            # selection
            mask   =   data['z_source_mock']>0.3
            mask   =   mask&(data['z_source_mock']<1.5)
            mask   =   mask&(data['e1_mock']**2.+data['e2_mock']**2.<10.)
            data   =   data[mask]
            weights=   weights[mask]['i_hsmshaperegauss_derived_weight']
            dataG.append(data)
            wG.append(weights)
            del data,mask
            gc.collect()
        dataG=  np.hstack(dataG)
        wG  =   np.hstack(wG)

        ra      =   dataG['ra_mock']
        dec     =   dataG['dec_mock']
        eI1     =   dataG['noise1_int']+dataG['noise1_mea']
        eI2     =   dataG['noise2_int']+dataG['noise2_mea']
        # define catalog (e,e_noi,g)
        catE0   =   treecorr.Catalog(g1=eI1,g2=eI2,ra=ra,dec=dec,w=wG,ra_units='deg',dec_units='deg')
        catE1   =   treecorr.Catalog(g1=dataG['e1_mock'],g2=dataG['e2_mock'],w=wG,ra=ra,dec=dec,ra_units='deg',dec_units='deg')
        catG    =   treecorr.Catalog(g1=dataG['shear1_sim'],g2=dataG['shear2_sim'],w=wG,ra=ra,dec=dec,ra_units='deg',dec_units='deg')
        # correlations
        corE0   =   treecorr.GGCorrelation(nbins=30,min_sep=0.281838293,max_sep=354.8133892,sep_units='arcmin')
        corE0.process(catE0)
        corE0.write(os.path.join(outDir,'rotmat%d_corE0.csv' %irot))

        corG    =   treecorr.GGCorrelation(nbins=30,min_sep=0.281838293,max_sep=354.8133892,sep_units='arcmin')
        corG.process(catG)
        corG.write(os.path.join(outDir,'rotmat%d_corG.csv' %irot))

        corGE0  =   treecorr.GGCorrelation(nbins=30,min_sep=0.281838293,max_sep=354.8133892,sep_units='arcmin')
        corGE0.process(catG,catE0)
        corGE0.write(os.path.join(outDir,'rotmat%d_corGE0.csv' %irot))

        corE1   =   treecorr.GGCorrelation(nbins=30,min_sep=0.281838293,max_sep=354.8133892,sep_units='arcmin')
        corE1.process(catE1)
        corE1.write(os.path.join(outDir,'rotmat%d_corE1.csv' %irot))

        del dataG,catG,catE0,catE1
        del corE0,corG,corGE0,corE1
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
                        default=108,
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
