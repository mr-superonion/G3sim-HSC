#!/usr/bin/env python
import os
import sys
import matplotlib
matplotlib.use('agg')
matplotlib.rcParams['ps.useafm'] = True
matplotlib.rcParams['pdf.use14corefonts'] = True
matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt

import scipy.optimize as opt
import weightedstats as ws
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


class reGaus_shearTestConfig(pexConfig.Config):
    anaDir  =   pexConfig.Field(dtype=str, default='anaCat', doc = 'simulation catalogs')
    obsDir  =   pexConfig.Field(dtype=str, default='/work/xiangchong.li/work/S16AFPFS/final/noBg/field/', doc = 'observed catalogs')
    pltDir  =   pexConfig.Field(dtype=str, default='plot/shearTest/', doc = 'plot directory')
    nbin_1d =   pexConfig.Field(dtype=int, default=8, doc = 'number of bins for 1-D parameter space')
    use_model=  pexConfig.Field(dtype=bool, default=True, doc = 'whether to use model to estimate weight')
    cancel_shape_noise=  pexConfig.Field(dtype=bool, default=True, doc = 'whether to use orthogonal galaxies to cancel shape noise')
    separate_components= pexConfig.Field(dtype=bool, default=False, doc = 'whether to seperate g1 g2 components')
    bin_plot= pexConfig.Field(dtype=bool, default=False, doc = 'whether to seperate g1 g2 components')
    
    max_snr =   pexConfig.Field(dtype=float, default=100000., doc = 'Maximum value of signal to noise ratio')
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

class reGaus_shearTestRunner(TaskRunner):
    @staticmethod
    def getTargetList(parsedCmd, **kwargs):
        return [(ref, kwargs) for ref in range(1)] 

def unpickle(factory, args, kwargs):
    """Unpickle something by calling a factory"""
    return factory(*args, **kwargs)

class reGaus_shearTestTask(BatchPoolTask):
    ConfigClass     =   reGaus_shearTestConfig
    RunnerClass     =   reGaus_shearTestRunner
    _DefaultName    =   "reGaus_shearTest"
    
    def __init__(self,**kwargs):
        BatchPoolTask.__init__(self, **kwargs)
        self.max_snr = self.config.max_snr
        self.max_e = self.config.max_e
        self.max_sigma_e = self.config.max_sigma_e
        self.lower_psf_size = self.config.lower_psf_size
        self.upper_psf_size = self.config.upper_psf_size

        # Initialize array to store multiplicative bias and associated quantities.
        self.snr_bins = np.zeros((self.config.nbin_1d, self.config.nbin_1d))
        self.res_bins = np.zeros((self.config.nbin_1d, self.config.nbin_1d))
        self.m1_bins = np.zeros((self.config.nbin_1d, self.config.nbin_1d))
        self.a1_bins = np.zeros((self.config.nbin_1d, self.config.nbin_1d))
        self.m2_bins = np.zeros((self.config.nbin_1d, self.config.nbin_1d))
        self.a2_bins = np.zeros((self.config.nbin_1d, self.config.nbin_1d))
        self.sigm1_bins = np.zeros((self.config.nbin_1d, self.config.nbin_1d))
        self.siga1_bins = np.zeros((self.config.nbin_1d, self.config.nbin_1d))
        
        #Estimated g1,g2 for the whole isim
        self.g1 = []
        self.g2 = []
        #True g1,g2
        self.g1_true = []
        self.g2_true = []
        self.psf_e1 = []
        self.psf_e2 = []
        self.sigma_g = []
        self.weight_s16a = []
        
        self.log.info('Beginning shearTest analysis for the following scenario:')
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
        pool    =   Pool("reGaus_shearTest")
        pool.cacheClear()
        outcomes=   pool.map(self.process,utilities.field_keys)
        outcomes=   [x for x in outcomes if x is not None]
        outtable=   vstack(outcomes)
        self.finalise(outtable)
        return
        
    def process(self,cache,fieldKey):
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
                                         use_wl_flag=True,
                                         max_snr=self.max_snr,
                                         cancel_shape_noise=self.config.cancel_shape_noise,
                                         max_e=self.max_e,
                                         max_sigma_e=self.max_sigma_e)
                                         
        return data[mask]
        
    def finalise(self,catalog,reweight_s16a=True):
        """
        Compute shear biases and make plots
        """
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
        assert np.all(cond_mask[0::2]==cond_mask[1::2]),\
                "Mask for imaging conditions does not obey shape noise cancellation sanity check"
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
        if reweight_s16a:
            shape_weights *= catalog['reweight_s16a']
        if self.config.cancel_shape_noise:
            shape_weights[1::2] = shape_weights[0::2]
            assert np.all(catalog['reweight_s16a'][0::2]==catalog['reweight_s16a'][1::2]), \
                'Galaxies in pairs do not have the same reweight_s16a factors'
        
        # Sort out binning. 
        snr_lower, snr_upper, res_lower, res_upper = \
            utilities.sliding_window_def(snr, res,
            n_bin_1d=self.config.nbin_1d)
        # Store Mask for each (snr,res) bin
        bin_mask = {}
        for j in range(self.config.nbin_1d):
            for k in range(self.config.nbin_1d):
                # This line below doesn't preserve shape noise cancellation, since each galaxy in
                # the pair has the cut imposed separately.
                tmp_mask = (snr >= snr_lower[j]) & (snr < snr_upper[j]) & \
                    (res >= res_lower[j,k]) & (res < res_upper[j,k])
                # This line will re-impose the shape noise cancellation (by only using the first in
                # the pair for the cut).
                if self.config.cancel_shape_noise:
                    tmp_mask[1::2] = tmp_mask[0::2]
                mean_snr = np.sum(shape_weights[tmp_mask]*snr[tmp_mask]) / \
                    np.sum(shape_weights[tmp_mask])
                mean_res = np.sum(shape_weights[tmp_mask]*res[tmp_mask]) / \
                    np.sum(shape_weights[tmp_mask])
                self.snr_bins[j,k] = mean_snr
                self.res_bins[j,k] = mean_res
                bin_mask['%d_%d'%(j,k)] = tmp_mask.copy()

        # The shear estimation routine is built to work on simulation subfields
        # We need to iterate on all independent subfields
        subfield_id = catalog['isim']
        sids = np.unique(subfield_id)
        self.g1_bins = np.zeros((self.config.nbin_1d, self.config.nbin_1d, len(sids)))
        self.g2_bins = np.zeros((self.config.nbin_1d, self.config.nbin_1d, len(sids)))
        self.sigma_g_bins = np.zeros((self.config.nbin_1d, self.config.nbin_1d, len(sids)))

        resp_list = []
        #Estimate shear for each isim
        for tmp_ind, i in enumerate(sids):
            m = (subfield_id == i)
            # Get true shear
            g1_true, g2_true = utilities.get_true_shear(catalog[m][0])
            # Define the PSF shape and rotation angle to apply to get into that frame
            psf_e1, psf_e2 = utilities.get_psf_ellip(catalog[m][0])
            self.g1_true.append(g1_true)
            self.g2_true.append(g2_true)
            self.psf_e1.append(psf_e1)
            self.psf_e2.append(psf_e2)
            self.weight_s16a.append(catalog['reweight_s16a'][m][0])
            
            
            #For overall estimation
            if np.sum(shape_weights[m])>0:
                # Estimate the shear
                if i==1: verbose=True
                else: verbose=False
                g1, g2, resp = utilities.estimate_subfield_shear(
                    catalog[m], use_model=self.config.use_model,
                    cancel_shape_noise=self.config.cancel_shape_noise,
                    simple_resp=True, report_resp=True,
                    reweight_s16a=reweight_s16a,
                    roundtrip_data=None,
                    test_weight_bias=False,
                    verbose=verbose)
                # Get an approximate shape noise per component:
                # sigma_g / sqrt(N) ~ 0.25/sqrt(N)
                # Note: we have cancelled shape noise here!  So there is not really
                # this much noise.  But measurement error goes like 1/sqrt(N) as well,
                # and capturing that trend is really the important part.
                sigma_g = 0.25/np.sqrt(shape_weights[m].sum())
                self.g1.append(g1)
                self.g2.append(g2)
                resp_list.append(resp)
                self.sigma_g.append(sigma_g)
            else:
                self.g1.append(0)
                self.g2.append(0)
                resp_list.append(0)
                self.sigma_g.append(1e12)
            
            
            # Get the estimate in bins in SNR/res as well.
            for j in range(self.config.nbin_1d):
                for k in range(self.config.nbin_1d):
                    mm = bin_mask['%d_%d'%(j,k)] & (subfield_id == i)
                    # Some weird bins might have no objects in particular subfields.
                    # Don't use estimate_subfield_shear on those!
                    if mm.astype(int).sum() > 2:
                        # Another criterion: make sure the objects have nonzero weight
                        if np.sum(shape_weights[mm])>0:
                            tmp_g1, tmp_g2, tmp_resp = utilities.estimate_subfield_shear(
                                catalog[mm], use_model=self.config.use_model,
                                cancel_shape_noise=self.config.cancel_shape_noise,
                                simple_resp=True, report_resp=True,
                                reweight_s16a=reweight_s16a,
                                roundtrip_data=None,
                                test_weight_bias=False)
                            self.g1_bins[j,k,tmp_ind] = tmp_g1
                            self.g2_bins[j,k,tmp_ind] = tmp_g2
                            self.sigma_g_bins[j,k,tmp_ind] = 0.25/np.sqrt(mm.astype(int).sum())
                        else:
                            self.g1_bins[j,k,tmp_ind] = 0.0
                            self.g2_bins[j,k,tmp_ind] = 0.0
                            self.sigma_g_bins[j,k,tmp_ind] = 1.e12
                    else:
                        self.g1_bins[j,k,tmp_ind] = 0.0
                        self.g2_bins[j,k,tmp_ind] = 0.0
                        self.sigma_g_bins[j,k,tmp_ind] = 1.e12

        self.sigma_g = np.array(self.sigma_g)
        resp_list = np.array(resp_list)
        self.g1 = np.array(self.g1)
        self.g2 = np.array(self.g2)
        self.weight_s16a = np.array(self.weight_s16a)
        self.g1_true = np.array(self.g1_true)
        self.g2_true = np.array(self.g2_true)
        self.psf_e1 = np.array(self.psf_e1)
        self.psf_e2 = np.array(self.psf_e2)
        overall_use = self.sigma_g < 1e10
        self.log.info('For overall analysis, using %d subfields with nonzero total weight'%\
            overall_use.astype(int).sum())
        sigma_vals = self.sigma_g[overall_use]
        if reweight_s16a:
            # We are going to weight by 1/sigma_vals^2.
            # But self.weight_s16a should be a multiplier for the weight.
            # So we take sigma_vals, and divide by sqrt(self.weight_s16a),
            # controlling for the cases where it leads to division by zero.
            division_val = np.sqrt(self.weight_s16a[overall_use])
            division_val[division_val<1.e-6] = 1.e-6
            sigma_vals /= division_val
        
        self.log.info('Overall bias analysis:')
        if self.config.separate_components:
            m1, a1, sig_m1, sig_a1, m2, a2, sig_m2, sig_a2, fig1, fig2, fig3, fig4 = \
                utilities.get_m_c_a(self.g1[overall_use], self.g2[overall_use],
                                    self.g1_true[overall_use], self.g2_true[overall_use],
                                    self.psf_e1[overall_use], self.psf_e2[overall_use],
                                    sigma_vals, separate_components=self.config.separate_components)
            self.log.info('m1 = %.3f +/- %.3f' % (m1, sig_m1))
            self.log.info('a1 = %.3f +/- %.3f' % (a1, sig_a1))
            self.log.info('m2 = %.3f +/- %.3f' % (m2, sig_m2))
            self.log.info('a2 = %.3f +/- %.3f' % (a2, sig_a2))
            fig1.savefig(self.config.pltDir+'shear_bias_1.eps')
            fig2.savefig(self.config.pltDir+'shear_bias_2.eps')
            fig3.savefig(self.config.pltDir+'shear_additive_residual_1.eps')
            fig4.savefig(self.config.pltDir+'shear_additive_residual_2.eps')
        else:
            m1, a1, sig_m1, sig_a1, fig1, fig2 = \
                utilities.get_m_c_a(self.g1[overall_use], self.g2[overall_use],
                                    self.g1_true[overall_use], self.g2_true[overall_use],
                                    self.psf_e1[overall_use], self.psf_e2[overall_use],
                                    sigma_vals, separate_components=self.config.separate_components)

            self.log.info('m = %.3f +/- %.3f' % (m1, sig_m1))
            self.log.info('a = %.3f +/- %.3f' % (a1, sig_a1))
            fig1.savefig(self.config.pltDir+'shear_bias.eps')
            fig2.savefig(self.config.pltDir+'shear_additive_residual.eps')
            

        # Histogram of responsivity and weighted mean etc.
        fig = plt.figure()
        ax = fig.add_subplot(111)
        h, ed = np.histogram(np.array(resp_list[overall_use]),
                             bins=int(len(resp_list[overall_use])/20.),
                             weights = 1./sigma_vals**2, density=True,
                             range = np.percentile(resp_list[overall_use], [1.,99.]))
        x = 0.5*(ed[:-1] + ed[1:])
        median_val = ws.weighted_median(resp_list[overall_use], weights=1./sigma_vals**2)
        ax.plot(x,h)
        ax.plot([median_val,median_val],[np.min(h),np.max(h)], 'r--')
        plt.xlabel('responsivity')
        ax.set_title('Median responsivity: %.3f'%(median_val))
        plt.savefig(self.config.pltDir+'shear_responsivity.eps')
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
