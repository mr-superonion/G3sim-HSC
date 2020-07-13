#!/usr/bin/env python
import os
import gc
import sys

import numpy as np
import pandas as pd
import sim_analysis_utilities as utilities

import lsst.daf.base as dafBase
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase

from lsst.pipe.base import TaskRunner
from lsst.ctrl.pool.parallel import BatchPoolTask
from lsst.ctrl.pool.pool import Pool, abortOnError

class reGaus_sanityTestConfig(pexConfig.Config):
    anaDir      =   pexConfig.Field(dtype=str, default='anaCat_newS19Mask', doc = 'simulation catalogs')
    pltDir      =   pexConfig.Field(dtype=str, default='plot/sanityTest/', doc = 'simulation catalogs')
    obsname     =   pexConfig.Field(dtype=str, default='mag', doc = 'name of the selection function')
    def setDefaults(self):
        pexConfig.Config.setDefaults(self)
    def validate(self):
        pexConfig.Config.validate(self)
        # Checks that the paths provided are valid
        if not os.path.isdir(self.anaDir):
            raise RuntimeError("Simulation catalog directory %s does not exist "
                         "or is not a directory" % (self.anaDir))

class reGaus_sanityTestRunner(TaskRunner):
    @staticmethod
    def getTargetList(parsedCmd, **kwargs):
        return [(ref, kwargs) for ref in range(1)]
def unpickle(factory, args, kwargs):
    """Unpickle something by calling a factory"""
    return factory(*args, **kwargs)

class reGaus_sanityTestTask(BatchPoolTask):
    ConfigClass     =   reGaus_sanityTestConfig
    RunnerClass     =   reGaus_sanityTestRunner
    _DefaultName    =   "reGaus_sanityTest"
    def __init__(self,**kwargs):
        BatchPoolTask.__init__(self, **kwargs)
    def __reduce__(self):
        """Pickler"""
        return unpickle, (self.__class__, [], dict(config=self.config, name=self._name,
                parentTask=self._parentTask, log=self.log))


    def determineBinMask(self,dataU,obsname):
        slide_test_bins =   np.zeros(len(dataU),dtype=np.uint64)
        assert np.all(dataU[0::2]['ipos']==(dataU['ipos'][1::2]-1))
        assert np.all(dataU[0::2]['weight']==(dataU['weight'][1::2]))
        if obsname=="mag":
            min_obs=21.5
            max_obs=24.5
            nbin_1d=3
            obs_bounds  =   np.linspace(min_obs,max_obs,nbin_1d+1)
        elif obsname=="magA10":
            min_obs=22.5
            max_obs=25.5
            nbin_1d=3
            obs_bounds  =   np.linspace(min_obs,max_obs,nbin_1d+1)
        elif obsname=="res":
            min_obs=0.3
            max_obs=1.0
            nbin_1d=3
            obs_bounds  =   np.linspace(min_obs,max_obs,nbin_1d+1)
        elif obsname=="snr":
            min_obs=10.
            max_obs=100.
            nbin_1d=3
            obs_bounds  =   np.linspace(min_obs,max_obs,nbin_1d+1)
        elif obsname=="COSMOS_Z":
            obs_bounds  =   np.array([0.0,0.5,1.0,1.5,2.5,4.])
            nbin_1d     =   len(obs_bounds)-1
        elif obsname=="COSMOS15_Z":
            obs_bounds  =   np.array([0.0,0.5,1.0,1.5,2.5,4.])
            nbin_1d     =   len(obs_bounds)-1
        elif obsname=="mizuki_Z":
            obs_bounds  =   np.array([0.0,0.5,1.0,1.5,2.5,4.])
            nbin_1d     =   len(obs_bounds)-1
        elif obsname=="COSMOS_Mag":
            min_obs=21.5
            max_obs=24.5
            nbin_1d=3
            obs_bounds  =   np.linspace(min_obs,max_obs,nbin_1d+1)
        elif obsname=="fwhm":
            min_obs=0.3
            max_obs=0.9
            nbin_1d=3
            obs_bounds  =   np.linspace(min_obs,max_obs,nbin_1d+1)
        obs_lowers  =   obs_bounds[:-1]
        obs_uppers  =   obs_bounds[1:]

        for iu in range(nbin_1d):
            tmp_msk     =   (dataU[obsname]>obs_lowers[iu])
            tmp_msk     &=  (dataU[obsname]<obs_uppers[iu])
            tmp_msk[1::2]=  tmp_msk[0::2]
            slide_test_bins[tmp_msk]=(slide_test_bins[tmp_msk]|np.uint(2**iu))
        return slide_test_bins,nbin_1d

    @abortOnError
    def runDataRef(self,Id):
        obsname =   self.config.obsname
        pool    =   Pool("reGaus_sanityTest")
        pool.cacheClear()
        pool.storeSet(anaDir=self.config.anaDir,
                obsname=obsname)

        indexList=  range(2500)
        outcomes=   pool.map(self.process,indexList)
        outcomes=   [x for x in outcomes if x is not None]
        outFname=   os.path.join(self.config.pltDir,'shear_sanitytest_%s_90rot-2.npy' %obsname)
        np.save(outFname,outcomes)
        return

    def process(self,cache,isim):
        # Get catalog
        inFname =   os.path.join(cache.anaDir,'allSim.h5')
        try:
            catalog =   pd.read_hdf(inFname).xs(isim,level='isim').to_records(index=True)
        except:
            self.log.info('Do not have sim-%05d' %isim)
            return None
        if len(catalog)<2:
            return None
        fN      =   catalog[0]['fieldname']
        self.log.info('Beginning for field: %s, sim: %s' %(fN,isim))

        self.log.info('Applying WL mask')
        obsname =   cache.obsname
        wlmaskS =   catalog['wlmask_90rot']
        catalog =   catalog[wlmaskS]

        self.log.info('Determining bin number')
        binMsk,nbin_1d  =   self.determineBinMask(catalog,obsname)

        # Do a round trip
        model_m =   utilities.get_m_model(catalog,weight_bias=True)
        model_c1,model_c2 = utilities.get_c_model(catalog,weight_bias=True)

        # The shear estimation routine is built to work on simulation subfields
        # We need to iterate on all independent subfields
        g1_bins = np.zeros(nbin_1d)
        g2_bins = np.zeros(nbin_1d)
        sigma_g_bins= np.zeros(nbin_1d)

        # Get true shear
        g1_true, g2_true= utilities.get_true_shear(catalog[0])
        # Define the PSF shape
        psf_e1, psf_e2  =   utilities.get_psf_ellip(catalog[0])
        reweight    =   catalog['weight'][0]
        for j in range(nbin_1d):
            # 90_degrees_rotation has already been forced
            mm  =   (binMsk&(2**j))!=0
            if mm.astype(int).sum()>4 and catalog['weight'][0]>0:
                roundtrip_data = (model_m[mm], model_c1[mm], model_c2[mm])
                tmp_g1, tmp_g2, tmp_resp = utilities.estimate_subfield_shear(
                    catalog[mm], use_model=True,
                    simple_resp=True, report_resp=True,
                    reweight=True,
                    roundtrip_data=roundtrip_data,
                    force_weight=False)
                g1_bins[j] = tmp_g1
                g2_bins[j] = tmp_g2
                sigma_g_bins[j] = 0.25/np.sqrt(mm.astype(int).sum())
            else:
                g1_bins[j] = 0.0
                g2_bins[j] = 0.0
                sigma_g_bins[j] = 1.e12

        oDat    =   (g1_true,g2_true,g1_bins,g2_bins,sigma_g_bins,psf_e1,psf_e2,reweight)
        return oDat

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
