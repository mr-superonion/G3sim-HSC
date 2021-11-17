#!/usr/bin/env python
import os
import glob

import numpy as np
import pandas as pd
import fitsio
import sim_analysis_utilities as utilities

import lsst.daf.base as dafBase
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase

from lsst.pipe.base import TaskRunner
from lsst.ctrl.pool.parallel import BatchPoolTask
from lsst.ctrl.pool.pool import Pool, abortOnError

class reGaus_selBiasConfig(pexConfig.Config):
    anaDir      =   pexConfig.Field(dtype=str, default='anaCat_newS19Mask_fdeltacut', doc = 'simulation catalogs')
    pltDir      =   pexConfig.Field(dtype=str, default='anaCat_newS19Mask_fdeltacut/plot/selBias/',
                        doc = 'simulation catalogs')
    obsname     =   pexConfig.Field(dtype=str, default='mag', doc = 'name of the selection function')
    subname     =   pexConfig.Field(dtype=str, default='mizuki_Z', doc = 'name of the selection function')
    nbin_1d     =   pexConfig.Field(dtype=int, default=5, doc = 'number of bins for 1-D parameter space')
    force_90_pair     =   pexConfig.Field(dtype=bool, default=False, doc = 'wheter to force 90 degree galaxies')
    def setDefaults(self):
        pexConfig.Config.setDefaults(self)
    def validate(self):
        pexConfig.Config.validate(self)
        # Checks that the paths provided are valid
        if not os.path.isdir(self.anaDir):
            raise RuntimeError("Simulation catalog directory %s does not exist "
                 "or is not a directory" % (self.anaDir))
        if not os.path.isdir(self.pltDir):
            os.system('mkdir -p %s' %self.pltDir)

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
    def __reduce__(self):
        """Pickler"""
        return unpickle, (self.__class__, [], dict(config=self.config, name=self._name,
                parentTask=self._parentTask, log=self.log))

    def determineMask(self,data,obsname,force_90_pair):
        # Read catalog
        # Apply WL cut for simulations
        max_mag =   24.5
        max_magA10= 25.5
        min_res =   0.3
        min_snr =   10
        max_snr =   None
        max_logb=   -0.38
        max_sigma_e=0.4
        max_e   =   2.

        if not force_90_pair:
            if obsname=="mag":
                max_mag =   None
            elif obsname=="magA10":
                max_magA10= None
            elif obsname=="res":
                min_res =   None
            elif obsname=="snr":
                min_snr =   None

        wlmaskS    =   utilities.galaxy_selector(
                        catalog     =   data,
                        min_snr     =   min_snr, #snr
                        max_snr     =   max_snr,
                        max_magA10  =   max_magA10, #magA10
                        min_res     =   min_res, #res
                        max_mag     =   max_mag, #mag
                        max_logb    =   max_logb,#logb
                        force_90_pair=  True,
                        max_e       =   max_e,     #absE
                        max_sigma_e =   max_sigma_e #sigma_e
                        )
        return wlmaskS

    def determineBinMask(self,dataU,subname,nbin_1d):
        if subname=="mag":
            obs_lowers  =   np.linspace(20,22.5,nbin_1d)
            obs_uppers  =   100.*np.ones(nbin_1d)
        elif subname=="magA10":
            obs_lowers  =   np.linspace(23,24.5,nbin_1d)
            obs_uppers     =   100.*np.ones(nbin_1d)
        elif subname=="res":
            obs_lowers  =   -1.*np.ones(nbin_1d)
            obs_uppers  =   np.linspace(0.5,1.,nbin_1d)
        elif subname=="snr":
            obs_lowers  =   -10.*np.ones(nbin_1d)
            obs_uppers  =   np.linspace(80.,200.,nbin_1d)
        elif subname=="mizuki_Z":
            obs_lowers  =   np.array([0.0,0.5,1.0,1.5,2.5])
            obs_uppers  =   np.array([0.5,1.0,1.5,2.5,4.0])
        return self.getBinNum(dataU,obs_lowers,obs_uppers,subname)

    def getBinNum(self,dataU,cut_lower,cut_upper,subname):
        """
        Parameters
        dataU: data after basic wl cut without the cut to test
        cut: the cut at which the selection bias is estimated
        cuts:the cuts at the other side (will force_90_pair)
        """
        slide_test_bins =   np.zeros(len(dataU),dtype=np.uint64)
        iu  =   0
        for lc,uc in zip(cut_lower,cut_upper):
            tmp_msk     =   (dataU[subname]>=lc)
            tmp_msk    &=   (dataU[subname]<uc)
            # force orthogonal to remove selection bias
            tmp_msk[1::2]=  tmp_msk[0::2]
            slide_test_bins[tmp_msk]=(slide_test_bins[tmp_msk]|np.uint(2**iu))
            iu  +=  1
        return slide_test_bins

    @abortOnError
    def runDataRef(self,Id):
        obsname =   self.config.obsname
        subname =   self.config.subname
        pool    =   Pool("reGaus_selBias")
        pool.cacheClear()
        pool.storeSet(anaDir=self.config.anaDir,
                obsname=obsname,subname=subname,
                nbin_1d=self.config.nbin_1d,
                force_90_pair=self.config.force_90_pair)

        gfn     =   os.path.join(self.config.anaDir,'*/smallcats/*.fits')
        fnList  =   glob.glob(gfn)
        outcomes=   pool.map(self.process,fnList)
        outcomes=   [x for x in outcomes if x is not None]
        if not self.config.force_90_pair:
            outFname=   os.path.join(self.config.pltDir,'shear_selbias_%s_%s.npy'
                    %(obsname,subname))
        else:
            outFname=   os.path.join(self.config.pltDir,'shear_selbias_%s_%s_90rot.npy'
                    %(obsname,subname))
        np.save(outFname,outcomes)
        return

    def process(self,cache,fname):
        # Get catalog
        catalog =   fitsio.read(fname)
        if len(catalog)<2:
            return None
        self.log.info('Processing file: %s' %(fname))

        self.log.info('Applying WL mask')
        obsname =   cache.obsname
        subname =   cache.subname
        nbin_1d =   cache.nbin_1d
        force_90_pair =   cache.force_90_pair
        wlmaskS =   self.determineMask(catalog,obsname,force_90_pair)
        wlmaskS =   wlmaskS&(~np.isnan(utilities.get_photo_z(catalog,'mizuki')))
        catalog =   catalog[wlmaskS]
        assert np.all(catalog['ipos'][0::2]==catalog['ipos'][1::2]-1)
        self.log.info('Determining bin number')
        binMsk  =   self.determineBinMask(catalog,subname,nbin_1d)
        if not force_90_pair:
            if obsname=="mag":
                maskCut=(catalog["mag"]-catalog['a_i'])<=24.5
            elif obsname=="magA10":
                maskCut=catalog["magA10"]<=25.5
            elif obsname=="res":
                maskCut=catalog["res"]>=0.3
            elif obsname=="snr":
                maskCut=catalog["snr"]>=10.
            catalog=catalog[maskCut]
            binMsk=binMsk[maskCut]

        # Do a round trip
        model_m =   utilities.get_m_model(catalog,weight_bias=True,photo_z_dep=True)
        model_c1,model_c2 = utilities.get_c_model(catalog,weight_bias=True,photo_z_dep=True)

        # The shear estimation routine is built to work on simulation subfields
        # We need to iterate on all independent subfields
        g1_bins = np.zeros(self.config.nbin_1d)
        g2_bins = np.zeros(self.config.nbin_1d)
        sigma_g_bins= np.zeros(self.config.nbin_1d)

        # Get true shear
        g1_true, g2_true= utilities.get_true_shear(catalog[0])
        # Define the PSF shape
        psf_e1, psf_e2  =   utilities.get_psf_ellip(catalog[0])
        reweight    =   catalog['weight'][0]
        for j in range(self.config.nbin_1d):
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
