#!/usr/bin/env python
import os
import sys
import matplotlib
matplotlib.use('agg')
matplotlib.rcParams['ps.useafm'] = True
matplotlib.rcParams['pdf.use14corefonts'] = True
matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt

import scipy.stats as stats
import sim_analysis_utilities as utilities
import numpy as np
from astropy.table import Table, vstack

import lsst.daf.base as dafBase
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase

from lsst.pipe.base import TaskRunner
from lsst.ctrl.pool.parallel import BatchPoolTask
from lsst.ctrl.pool.pool import Pool, abortOnError


class reGaus_shearMcConfig(pexConfig.Config):
    anaDir  =   pexConfig.Field(dtype=str, default='anaCat', doc = 'simulation catalogs')
    obsDir  =   pexConfig.Field(dtype=str, default='/work/xiangchong.li/work/S16AFPFS/final/noBg/field/', doc = 'observed catalogs')
    pltDir  =   pexConfig.Field(dtype=str, default='plot/shearMc/', doc = 'plot directory')
    nbin_1d =   pexConfig.Field(dtype=int, default=8, doc = 'number of bins for 1-D parameter space')
    use_model=  pexConfig.Field(dtype=bool, default=True, doc = 'whether to use model to estimate weight')
    
    max_snr =   pexConfig.Field(dtype=float, default=80., doc = 'Maximum value of signal to noise ratio')
    max_e   =   pexConfig.Field(dtype=float, default=2., doc = 'Maximum value of distortion |e|')
    max_sigma_e =   pexConfig.Field(dtype=float, default=0.4, doc = 'Maximum value of sigmaE')
    lower_psf_size  =   pexConfig.Field(dtype=float, default=0., doc = 'Minimum of PSF size')
    upper_psf_size  =   pexConfig.Field(dtype=float, default=100., doc = 'Maximum of PSF size')
    
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

class reGaus_shearMcRunner(TaskRunner):
    @staticmethod
    def getTargetList(parsedCmd, **kwargs):
        return [(ref, kwargs) for ref in range(1)] 

def unpickle(factory, args, kwargs):
    """Unpickle something by calling a factory"""
    return factory(*args, **kwargs)

class reGaus_shearMcTask(BatchPoolTask):
    ConfigClass     =   reGaus_shearMcConfig
    RunnerClass     =   reGaus_shearMcRunner
    _DefaultName    =   "reGaus_shearMc"
    
    def __init__(self,**kwargs):
        BatchPoolTask.__init__(self, **kwargs)
        self.max_snr = self.config.max_snr
        self.max_e = self.config.max_e
        self.max_sigma_e = self.config.max_sigma_e
        self.lower_psf_size = self.config.lower_psf_size
        self.upper_psf_size = self.config.upper_psf_size

        # Initialize array to store multiplicative bias and associated quantities.
        self.m1_bins = np.zeros((self.config.nbin_1d, self.config.nbin_1d))
        self.a1_bins = np.zeros((self.config.nbin_1d, self.config.nbin_1d))
        self.m2_bins = np.zeros((self.config.nbin_1d, self.config.nbin_1d))
        self.a2_bins = np.zeros((self.config.nbin_1d, self.config.nbin_1d))
        self.sigm1_bins = np.zeros((self.config.nbin_1d, self.config.nbin_1d))
        self.siga1_bins = np.zeros((self.config.nbin_1d, self.config.nbin_1d))
        
        
        
        
        self.log.info('Beginning shearMc analysis for the following scenario:')
        self.log.info('   %d bins in both SNR and resolution'%self.config.nbin_1d)
        self.log.info('   Percentile ranges for PSF size: %.2f-%.2f'%\
            (self.lower_psf_size, self.upper_psf_size))
        self.log.info('')
        
    def __reduce__(self):
        """Pickler"""
        return unpickle, (self.__class__, [], dict(config=self.config, name=self._name,
                parentTask=self._parentTask, log=self.log))
    
    @abortOnError
    def run(self,Id):
        #Prepare the pool
        rows    =   Table.read('noiVar_rPsf_s18a.csv')
        pool    =   Pool("reGaus_shearMc")
        pool.cacheClear()
        outcomes=   pool.map(self.process,utilities.field_keys)
        outcomes    =   [x for x in outcomes if x is not None]
        outtable=   vstack(outcomes)
        self.finalise(outtable)
        return
        
        
    def process(self,cache,fieldKey,reweight_s16a=True,use_wl_flag=True):
        fieldname=  utilities.field_keys[fieldKey]
        if fieldname==  'AEGIS':
            return None
        self.log.info('Processing field %s' %(fieldname))
        fname   =   os.path.join(self.config.anaDir,fieldname,'%s.fits'%(fieldname))
        data =   Table.read(fname)
        """
        Get estimated and true shape information, lensing shear, and other information.
        """
        # Define a mask with some set of cuts.  Two notes:
        # 1) cancel_shape_noise MUST be true for this test (so we don't use
        #    self.cancel_shape_noise).  It's not an option that the user
        #    can change, because the 90 degree rotated pairs are required.
        # 2) bright_star_mask is set to False (instead of
        #    self.bright_star_mask) because we are running on sims.
        mask    =   utilities.galaxy_selector(catalog=data,
                                         use_wl_flag=use_wl_flag,
                                         bright_star_mask=False,
                                         max_snr=self.max_snr,
                                         cancel_shape_noise=True,
                                         max_e=self.max_e,
                                         max_sigma_e=self.max_sigma_e)
                                         
        return data[mask]
        
    def finalise(self,catalog):
        """
        Compute shear biases and make plots
        """
        # Concatenate all the catalogs across fields.
        catalog = astropy.table.vstack(self.catalog_list)

        # Get mask for imaging conditions.  Note this will always preserve shape noise cancellation
        # due to how the sims are set up, though we can do an assert just to check.
        psf_size = utilities.get_psf_size(catalog)
        cond_mask = np.ones_like(catalog).astype(bool)
        if self.lower_psf_size > 0.:
            cut_val = np.percentile(psf_size, self.lower_psf_size)
            cond_mask &= (psf_size >= cut_val)
        if self.upper_psf_size < 100.:
            cut_val = np.percentile(psf_size, self.upper_psf_size)
            cond_mask &= (psf_size < cut_val)
        catalog = catalog[cond_mask]
        
        # Set up sliding window definition and make a mask for each SNR, resolution bin.
        snr = utilities.get_snr(catalog)
        res = utilities.get_res(catalog)
        # We need to get the weights here, to use when averaging over galaxies to get mean effective
        # SNR and resolution.
        # Note, it might be worth making a utility for this combined shape + S16A weight, since it's
        # done here and in estimate_subfield_shear.
        if self.config.use_model: 
            shape_weights = utilities.get_weight_model(catalog)
        else:
            shape_weights = catalog['ishape_hsm_regauss_derived_shape_weight']
        
        
        
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
