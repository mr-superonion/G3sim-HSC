#!/usr/bin/env python
import os
import sys

import numpy as np
import pandas as pd
import scipy.stats as stats
import sim_analysis_utilities as utilities

import lsst.daf.base as dafBase
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase

from lsst.pipe.base import TaskRunner
from lsst.ctrl.pool.parallel import BatchPoolTask
from lsst.ctrl.pool.pool import Pool, abortOnError

class reGaus_selBiasConfig(pexConfig.Config):
    anaDir      =   pexConfig.Field(dtype=str, default='anaCat', doc = 'simulation catalogs')
    pltDir      =   pexConfig.Field(dtype=str, default='plot/reGausCalib', doc = 'simulation catalogs')
    nbin_1d     =   pexConfig.Field(dtype=int, default=8, doc = 'number of bins for 1-D parameter space')
    use_model   =   pexConfig.Field(dtype=bool, default=True, doc = 'whether to use model to estimate weight')
    reweight    =  pexConfig.Field(dtype=bool, default=True, doc = 'reweigting?')
    force_weight=   pexConfig.Field(dtype=bool, default=False, doc = 'for weight bias')
    test_weight_bias  =   pexConfig.Field(dtype=bool, default=True, doc = 'test weight bias?')

    def setDefaults(self):
        pexConfig.Config.setDefaults(self)
    def validate(self):
        pexConfig.Config.validate(self)
        if self.test_weight_bias and self.force_weight:
            raise RuntimeError("Cannot test weight bias when forcing weight")
        # Checks that the paths provided are valid
        if not os.path.isdir(self.anaDir):
            raise RuntimeError("Simulation catalog directory %s does not exist "
                         "or is not a directory" % (self.anaDir))

class reGaus_selBiasRunner(TaskRunner):
    @staticmethod
    def getTargetList(parsedCmd, **kwargs):
        return [(ref, kwargs) for ref in range(1)]
def unpickle(factory, args, kwargs):
    """Unpickle something by calling a factory"""
    return factory(*args, **kwargs)

class reGaus_selBiasTask(BatchPoolTask):
    ConfigClass     =   reGaus_selBiasConfig
    RunnerClass     =   reGaus_selBiasRunner
    _DefaultName    =   "reGaus_selBias"

    def __init__(self,**kwargs):
        BatchPoolTask.__init__(self, **kwargs)

        # Sliding bins
        self.snr_bins   =   np.zeros(self.nbin_1d)
        self.res_bins   =   np.zeros(self.nbin_1d)

        self.log.info(' Beginning shearCalib analysis for the following scenario:')
        self.log.info('  %d bins in both SNR and resolution'%self.config.nbin_1d)

    def __reduce__(self):
        """Pickler"""
        return unpickle, (self.__class__, [], dict(config=self.config, name=self._name,
                parentTask=self._parentTask, log=self.log))

    @abortOnError
    def run(self,Id):
        self.log.info('Estimating bias for each bin')
        inFname =   os.path.join(self.config.anaDir,'allSim.h5')
        dataF   =   pd.read_hdf(inFname)
        dataF   =   dataF[dataF['wlmask_90rot']]
        pool    =   Pool("reGaus_selBias")
        pool.cacheClear()
        pool.storeSet(dataF=dataF)
        indexList=  dataF.index.levels[1]
        outcomes=   pool.map(self.process,indexList)
        outcomes=   [x for x in outcomes if x is not None]
        outFname=   os.path.join(self.config.pltDir,'shear_selbias_res.npy')
        np.save(outFname,outcomes)
        return

    def process(self,cache,index):
        # Get catalog
        catalog =   cache.dataF.xs(index,level='isim').to_records(index=True)
        if len(catalog)<2:
            return None
        fN=catalog[0]['fieldname']
        isim=index
        self.log.info('Beginning for field: %s, sim: %s' %(fN,isim))
        shape_weights = utilities.get_weight_model(catalog)*catalog['weight']

        if self.config.force_weight:
            shape_weights[1::2] = shape_weights[0::2]

        # Do a round trip
        model_m = utilities.get_m_model(catalog,weight_bias=self.config.test_weight_bias)
        model_c1,model_c2 = utilities.get_c_model(catalog,weight_bias=self.config.test_weight_bias)

        # The shear estimation routine is built to work on simulation subfields
        # We need to iterate on all independent subfields
        g1_bins = np.zeros(self.nbin_1d)
        g2_bins = np.zeros(self.nbin_1d)
        sigma_g_bins= np.zeros(self.nbin_1d)

        # Get true shear
        g1_true, g2_true= utilities.get_true_shear(catalog[0])
        # Define the PSF shape
        psf_e1, psf_e2  =   utilities.get_psf_ellip(catalog[0])
        reweight     =   catalog['weight'][0]

        # Get the estimate in bins in SNR/res as well.
        for j in range(self.config.nbin_1d):
            for k in range(self.config.nbin_1d):
                # 90_degrees_rotation has already been forced
                sbin=j*self.config.nbin_1d+k
                mm  =   (catalog['calib_slide_ibins']==sbin)
                if mm.astype(int).sum()>4 and np.sum(shape_weights[mm])>0:
                    roundtrip_data = (model_m[mm], model_c1[mm], model_c2[mm])
                    tmp_g1, tmp_g2, tmp_resp = utilities.estimate_subfield_shear(
                        catalog[mm], use_model=self.config.use_model,
                        simple_resp=True, report_resp=True,
                        reweight=self.config.reweight,
                        roundtrip_data=roundtrip_data,
                        force_weight=self.config.force_weight)
                    g1_bins[j,k] = tmp_g1
                    g2_bins[j,k] = tmp_g2
                    sigma_g_bins[j,k] = 0.25/np.sqrt(mm.astype(int).sum())
                else:
                    g1_bins[j,k] = 0.0
                    g2_bins[j,k] = 0.0
                    sigma_g_bins[j,k] = 1.e12

        oDat    =   (g1_true,g2_true,g1_bins,g2_bins,sigma_g_bins,
                    psf_e1,psf_e2,reweight)
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
