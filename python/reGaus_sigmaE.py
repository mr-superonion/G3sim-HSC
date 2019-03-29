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


"""
This script is for tests of the catalog sigma_e values.
"""

import math
import numpy as np
import scipy.stats as stats
import sim_analysis_utilities as utilities
import matplotlib.pyplot as plt
import matplotlib


import lsst.daf.base as dafBase
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase

from lsst.pipe.base import TaskRunner
from lsst.ctrl.pool.parallel import BatchPoolTask
from lsst.ctrl.pool.pool import Pool, abortOnError


matplotlib.rcParams['ps.useafm'] = True
matplotlib.rcParams['pdf.use14corefonts'] = True
matplotlib.rcParams['text.usetex'] = True



def do_ks_test(x,y,sample_names="two samples"):
    """Routine to do a 2-sample KS test and report on consistency of distributions"""

    (ks_d, ks_p) = stats.ks_2samp(x,y)
    if ks_p >= 0.05:
        self.log.info('Distributions for %s are consistent!'%sample_names)
    else:
        self.log.info('Distributions for %s differ!'%sample_names)
    self.log.info('(d, p)=(%f, %f)'%(ks_d, ks_p))

class reGaus_sigmaE_Config(pexConfig.Config):
    nbin_1d =   pexConfig.Field(dtype=int, default=20, doc = 'number of bins for 1-D parameter space')
    max_snr =   pexConfig.Field(dtype=float, default=80., doc = 'Maximum value of signal to noise ratio')
    max_e   =   pexConfig.Field(dtype=float, default=2., doc = 'Maximum value of distortion |e|')
    max_sigma_e =   pexConfig.Field(dtype=float, default=0.4, doc = 'Maximum value of sigmaE')
    lower_psf_size  =   pexConfig.Field(dtype=float, default=0., doc = 'Minimum of PSF size')
    upper_psf_size  =   pexConfig.Field(dtype=float, default=100., doc = 'Maximum of PSF size')
    def setDefaults(self):
        pexConfig.Config.setDefaults(self)

    def validate(self):
        pexConfig.Config.validate(self)


class reGaus_sigmaE_Task(pipeBase.CmdLineTask):
    
    """
    Infer shape measurement error distribution from 90 degree rotated pairs; compare with the one
    that is implied from the catalog; and build up a model for the effective sigma_e as a function
    of galaxy parameters.
    """
    _DefaultName = "reGaus_sigmaE"
    ConfigClass = reGaus_sigmaE_Config
    
    def __init__(self,**kwargs):
        pipeBase.CmdLineTask.__init__(self, **kwargs)
    
        """
        Infer shape measurement error distribution (at least initially over all subfields) and
        compare with catalog information.
        """

        self.max_snr = self.config.max_snr
        self.max_e = self.config.max_e
        self.max_sigma_e = self.config.max_sigma_e
        self.lower_psf_size = self.config.lower_psf_size
        self.upper_psf_size = self.config.upper_psf_size

        # Initialise lists to store data
        self.e1 = []
        self.e2 = []
        self.sigma_e = []
        self.g1_true = []
        self.g2_true = []
        self.resolution = []
        self.snr = []
        self.weight = []
        # Initialize array to store effective sigma_e values and associated quantities.
        self.cat_est_sigma_e = np.zeros((self.config.nbin_1d, self.config.nbin_1d))
        self.est_sigma_e = np.zeros((self.config.nbin_1d, self.config.nbin_1d))
        self.higherorder_var_e = np.zeros((self.config.nbin_1d, self.config.nbin_1d))
        self.mean_res = np.zeros((self.config.nbin_1d, self.config.nbin_1d))
        self.mean_snr = np.zeros((self.config.nbin_1d, self.config.nbin_1d))
        self.true_rms = np.zeros((self.config.nbin_1d, self.config.nbin_1d))
        self.est_rms = np.zeros((self.config.nbin_1d, self.config.nbin_1d))
        self.cat_est_rms = np.zeros((self.config.nbin_1d, self.config.nbin_1d))
        self.corr_coeff = np.zeros((self.config.nbin_1d, self.config.nbin_1d))
        self.psf_size = []

        self.log.info('Beginning sigma_e analysis for the following scenario:')
        self.log.info('   %d bins in both SNR and resolution'%self.config.nbin_1d)
        self.log.info('   Percentile ranges for PSF size: %.2f-%.2f'%\
            (self.lower_psf_size, self.upper_psf_size))
        self.log.info('')

    @pipeBase.timeMethod
    def run(self,ifield,fieldname,reweight_s16a=False):
        """
        Get estimated and true shape information, lensing shear, and other information.
        """
        # Define a mask with some set of cuts.  Two notes:
        # 1) cancel_shape_noise MUST be true for this test (so we don't use
        #    self.cancel_shape_noise).  It's not an option that the user
        #    can change, because the 90 degree rotated pairs are required.
        # 2) bright_star_mask is set to False (instead of
        #    self.bright_star_mask) because we are running on sims.
        data    =   
        mask    =   utilities.galaxy_selector(catalog=data,
                                         use_wl_flag=self.use_wl_flag,
                                         bright_star_mask=False,
                                         max_snr=self.max_snr,
                                         cancel_shape_noise=True,
                                         max_e=self.max_e,
                                         max_sigma_e=self.max_sigma_e)

        # Get true shear for the subfield.
        g1_true, g2_true = utilities.get_true_shear(data[mask])
        # Store NumPy arrays for observed and true e1, e2 for each subfield, and sigma_e
        self.e1.append(data[mask]['ishape_hsm_regauss_e1'])
        self.e2.append(data[mask]['ishape_hsm_regauss_e2'])
        self.sigma_e.append(utilities.get_sigma_e(data[mask])) # applies noise rescaling
        self.resolution.append(data[mask]['ishape_hsm_regauss_resolution'])
        self.snr.append(utilities.get_snr(data[mask]))
        if reweight_s16a:
            self.weight.append(data[mask]['reweight_s16a'])
        else:
            self.weight.append(np.ones_like(true_e1).astype(float))
        self.psf_size.append(utilities.get_psf_size(data[mask]))
        # Store subfield true shears.
        self.g1_true.append(g1_true)
        self.g2_true.append(g2_true)
        return

    def finalise(self,pltDir,pltType):
        """
        Make distributions of the quantities that we need to calculate
        in order to estimate the distribution of shape measurement
        errors.  Everything here is unweighted; the only exception is
        weighting to account for the S16A conditions.
        """
        # First: assemble the arrays for individual subfields into one big master catalog.
        n_subfields = len(self.e1)
        master_e1 = np.concatenate(self.e1)
        master_e2 = np.concatenate(self.e2)
        master_sigma_e = np.concatenate(self.sigma_e)
        master_resolution = np.concatenate(self.resolution)
        master_snr = np.concatenate(self.snr)
        master_weight = np.concatenate(self.weight)
        tmp_e = np.sqrt(master_e1**2 + master_e2**2)
        master_g1_true = np.concatenate(self.g1_true)
        master_g2_true = np.concatenate(self.g2_true)
        master_psf_size = np.concatenate(self.psf_size)


        # Figure out any other mask to apply for observing conditions.
        #     cond_mask is a mask for the *pairs*
        #     global_cond_mask is a mask for the *individual objects*
        cond_mask = np.ones_like(master_g1_true[0::2]).astype(bool)
        global_cond_mask = np.ones_like(master_g1_true).astype(bool)
        if self.lower_psf_size > 0.:
            cut_val = np.percentile(master_psf_size, self.lower_psf_size)
            cond_mask &= master_psf_size[0::2] >= cut_val
            global_cond_mask &= master_psf_size >= cut_val
        if self.upper_psf_size < 100.:
            cut_val = np.percentile(master_psf_size, self.upper_psf_size)
            cond_mask &= master_psf_size[0::2] < cut_val
            global_cond_mask &= master_psf_size < cut_val

        # Calculate LHS of equation A8 in writeup, for each component
        # separately.
        sum_e_erot_1 = master_e1[0::2] + master_e1[1::2]
        sum_e_erot_2 = master_e2[0::2] + master_e2[1::2]

        # Then start calculating the correction terms with signs given in equation 15.
        # We start with the shear itself.
        correction_1 = -4.*master_g1_true[0::2]
        correction_2 = -4.*master_g2_true[0::2]
        # Since this doesn't depend on the intrinsic shape, let's ALWAYS apply this correction.
        # Fortunately, this is the lowest order correction, so it's good we can always apply it.
        # Eq.(A9)
        sum_e_erot_1 += correction_1
        sum_e_erot_2 += correction_2
        sum_e_erot_1_corr = sum_e_erot_1.copy()
        sum_e_erot_2_corr = sum_e_erot_2.copy()
        # So technically the "sum_e_erot_1" vs. "sum_e_erot_1_corr" is a bit confusing.  BOTH
        # include the lowest order correction due to lensing.  The "corr" refers to whether it has
        # the higher order correction (which is not so important anyway).


        # So the sum_e_erot_#_corr should be distributed like the convolution of p_err(e) with
        # itself.  We can build up a theoretical expectation of what this looks like assuming the
        # Gaussian approximation.
        expected_sigma = np.concatenate((master_sigma_e[0::2][cond_mask],
                                         master_sigma_e[1::2][cond_mask]))

        # Do K-S test for consistency of distributions for the two components.
        do_ks_test(sum_e_erot_1_corr[cond_mask], sum_e_erot_2_corr[cond_mask],
                   sample_names="two shape components")

        # Do calculation overall.
        e_uncorr_input = np.concatenate((
                sum_e_erot_1[cond_mask],
                sum_e_erot_2[cond_mask]))
        e_input = np.concatenate((
                sum_e_erot_1_corr[cond_mask],
                sum_e_erot_2_corr[cond_mask]))
        weight_input = np.concatenate((
                master_weight[0::2][cond_mask],
                master_weight[0::2][cond_mask]))# only use the first of the first galaxy in a pair
        # Get RMS of underlying distribution, i.e., not the
        # auto-convolution but the original one, hence the 0.5 inside
        # the sqrt.
        rms = np.sqrt(0.5*np.sum(weight_input*e_input**2)/np.sum(weight_input))
        self.log.info('Overall N, <w>, RMS: ',len(weight_input), np.mean(weight_input), rms)

        # Get implied RMS from catalog.
        input_rms = np.sqrt(np.sum(weight_input*expected_sigma**2)/np.sum(weight_input))
        self.log.info('Implied from catalog: ',input_rms)

        # Now set up sliding window definition
        # Only use the first in pair
        snr_lower, snr_upper, res_lower, res_upper = \
            utilities.sliding_window_def(master_snr[0::2][cond_mask],
                                         master_resolution[0::2][cond_mask],
                                         n_bin_1d=self.config.nbin_1d)

        # Now loop over SNR and resolution bins, using the first
        # galaxy in a pair to stick them in a bin.
        for i in range(len(snr_lower)):
            for j in range(len(snr_lower)):
                self.log.info('Beginning bin: SNR in [%.3f, %.3f], resolution in [%.3f, %.3f]'%\
                    (snr_lower[i], snr_upper[i], res_lower[i,j], res_upper[i,j]))
                bin_mask = \
                    (master_snr[0::2] >= snr_lower[i]) & \
                    (master_snr[0::2] < snr_upper[i]) & \
                    (master_resolution[0::2] >= res_lower[i,j]) & \
                    (master_resolution[0::2] < res_upper[i,j]) & \
                    cond_mask
                global_bin_mask = \
                    (master_snr >= snr_lower[i]) & \
                    (master_snr < snr_upper[i]) & \
                    (master_resolution >= res_lower[i,j]) & \
                    (master_resolution < res_upper[i,j]) & \
                    global_cond_mask
                self.log.info('Using %d pairs'%bin_mask.astype(int).sum())
                self.mean_snr[i,j] = \
                    np.sum(master_weight[0::2][bin_mask]*master_snr[0::2][bin_mask])/ \
                    np.sum(master_weight[0::2][bin_mask])
                self.mean_res[i,j] = \
                    np.sum(master_weight[0::2][bin_mask]*master_resolution[0::2][bin_mask])/ \
                    np.sum(master_weight[0::2][bin_mask])
                # Sanity check weights.
                assert np.all(master_weight[0::2][bin_mask]==master_weight[1::2][bin_mask])

                # Check for consistency between components.
                do_ks_test(sum_e_erot_1_corr[bin_mask], sum_e_erot_2_corr[bin_mask],
                           sample_names="two shape components")
                
                e_uncorr_input = np.concatenate((
                        sum_e_erot_1[bin_mask],
                        sum_e_erot_2[bin_mask]))
                e_input = np.concatenate((
                        sum_e_erot_1_corr[bin_mask],
                        sum_e_erot_2_corr[bin_mask]))
                weight_input = np.concatenate((
                        master_weight[0::2][bin_mask],
                        master_weight[0::2][bin_mask]))
                # Get RMS of underlying distribution, i.e., not the auto-convolution but the original one
                bin_rms = np.sqrt(0.5*np.sum(weight_input*e_input**2)/np.sum(weight_input))
                self.log.info('N, <w>, RMS in bin: ',len(weight_input), np.mean(weight_input), bin_rms)

                # Get implied RMS from catalog.
                expected_sigma = np.concatenate((master_sigma_e[0::2][bin_mask],
                                                 master_sigma_e[1::2][bin_mask]))
                input_bin_rms = np.sqrt(
                    np.sum(weight_input*expected_sigma**2)/np.sum(weight_input))
                self.log.info('RMS in bin implied from catalog: ',input_bin_rms)

                self.est_sigma_e[i,j] = bin_rms
                self.cat_est_sigma_e[i,j] = input_bin_rms

                # Compute RMS e values per component, including:
                # 1) Observed
                # 2) Estimated intrinsic, using 
                #   (a) the catalog sigma_e values, and
                #   (b) the one from these distributions.
                #
                # (1) Observed
                rms_observed = np.sqrt(0.5*np.sum(
                        master_weight[global_bin_mask]* \
                            (master_e1[global_bin_mask]**2+master_e2[global_bin_mask]**2)) / \
                                           np.sum(master_weight[global_bin_mask]))
                # (2a) Estimated intrinsic using catalog sigma_e values.
                est_true_cat = np.sqrt(
                    0.5*np.sum(
                        master_weight[global_bin_mask]*\
                            (master_e1[global_bin_mask]**2+ \
                                    master_e2[global_bin_mask]**2- \
                                    2.*master_sigma_e[global_bin_mask]**2)) / \
                        np.sum(master_weight[global_bin_mask]))
                # (2b) Estimate intrinsic using the estimated sigma_e values
                est_true_dist = np.sqrt(
                    0.5*np.sum(
                        master_weight[global_bin_mask]*\
                            (master_e1[global_bin_mask]**2+ \
                                 master_e2[global_bin_mask]**2- \
                                 2.*bin_rms**2)) / \
                        np.sum(master_weight[global_bin_mask]))

                self.cat_est_rms[i,j] = est_true_cat
                self.est_rms[i,j] = est_true_dist
                self.log.info('RMS values and responsivities per component in this bin:')
                self.log.info('  observed, including shape measurement error: %.4f, %.4f'%\
                    (rms_observed, 2.*(1.-rms_observed**2)))
                self.log.info('  estimated using distribution to remove measurement error: %.4f, %.4f'%\
                    (est_true_dist, 2.*(1.-est_true_dist**2)))

        self.log.info('')
        self.log.info('')
        x_vals = np.vstack([np.log(self.mean_snr/20.).flatten(),
                            np.log(self.mean_res/0.5).flatten(),
                            np.ones(len(self.mean_res.flatten()))]).transpose()
        est_log_sig_e = np.log(self.est_sigma_e.flatten())
        result  = np.linalg.lstsq(x_vals, est_log_sig_e)
        #Eq.(18)
        model_sigma_e = np.exp(result[0][2])*((self.mean_snr/20.)**result[0][0])*\
            ((self.mean_res/0.5)**result[0][1])
        self.log.info('model sigma_e = %.3f (SNR/20)^%.3f (res/0.5)^%.3f'%\
            (np.exp(result[0][2]), result[0][0], result[0][1]))

        # Plot estimated sigma_e as a function of SNR, resolution.
        cm = plt.cm.get_cmap('RdYlBu')
        fig = plt.figure()
        ax = fig.add_subplot(111)
        vmin = min(np.min(self.est_sigma_e), np.min(model_sigma_e), np.min(self.cat_est_sigma_e))
        vmax = max(np.max(self.est_sigma_e), np.max(model_sigma_e), np.max(self.cat_est_sigma_e))
        sc = ax.scatter(self.mean_snr, self.mean_res, c=self.est_sigma_e, cmap=cm,
                        vmin=vmin, vmax=vmax)
        plt.colorbar(sc)
        ax.set_title('Estimated RMS measurement error')
        plt.xscale('log')
        plt.xlabel('SNR')
        plt.ylabel('resolution')
        plt.savefig(os.path.join(pltDir,'sigmae_estimated.%s'%(pltType)))


        # Plot catalog sigma_e as a function of SNR, resolution.
        fig = plt.figure()
        ax = fig.add_subplot(111)
        sc = ax.scatter(self.mean_snr, self.mean_res, c=self.cat_est_sigma_e, cmap=cm,
                        vmin=vmin, vmax=vmax)
        plt.colorbar(sc)
        ax.set_title('Catalog RMS measurement error')
        plt.xscale('log')
        plt.xlabel('SNR')
        plt.ylabel('resolution')
        plt.savefig(os.path.join(pltDir,'sigmae_cat.%s'%(pltType)) )

        # Plot ratio of catalog to estimated sigma_e.
        fig = plt.figure()
        ax = fig.add_subplot(111)
        sc = ax.scatter(self.mean_snr, self.mean_res, c=self.cat_est_sigma_e/self.est_sigma_e, cmap=cm)
        plt.colorbar(sc)
        ax.set_title('Catalog/estimated RMS measurement error')
        plt.xscale('log')
        plt.xlabel('SNR')
        plt.ylabel('resolution')
        plt.savefig(os.path.join(pltDir,'sigmae_cat_estimated_ratio.%s' %(pltType)))


        # Model sigma_e as a function of SNR, resolution (including
        # extrapolation of the model).
        fig = plt.figure()
        ax = fig.add_subplot(111)
        sc = ax.scatter(self.mean_snr, self.mean_res, c=model_sigma_e, cmap=cm,
                        vmin=vmin, vmax=vmax)
        plt.colorbar(sc)
        ax.set_title('Model measurement RMS')
        plt.xscale('log')
        plt.xlabel('SNR')
        plt.ylabel('resolution')
        plt.savefig(os.path.join(pltDir,'sigmae_model.%s'%(pltType)))
        # model / estimated sigma_e as a function of SNR, resolution.  Should be close to 1.
        fig = plt.figure()
        ax = fig.add_subplot(111)
        sc = ax.scatter(self.mean_snr, self.mean_res, c=self.est_sigma_e/model_sigma_e, cmap=cm)
        plt.colorbar(sc)
        ax.set_title('Estimated/model RMS')
        plt.xscale('log')
        plt.xlabel('SNR')
        plt.ylabel('resolution')
        plt.savefig(os.path.join(pltDir,'sigmae_ratio.%s' %(pltType)))

        # save data table for model/actual
        out_data = (self.mean_snr.flatten(), self.mean_res.flatten(),
                    (self.est_sigma_e/model_sigma_e).flatten())
        out_data = np.column_stack(out_data)
        np.savetxt(os.path.join(pltDir,'sigmae_ratio.dat'), out_data)
        return
        
    @classmethod
    def _makeArgumentParser(cls):
        """Create an argument parser
        """
        parser = pipeBase.ArgumentParser(name=cls._DefaultName)
        return parser

    def writeConfig(self, butler, clobber=False, doBackup=False):
        pass

    def writeSchemas(self, butler, clobber=False, doBackup=False):
        pass

    def writeMetadata(self, ifield):
        pass

    def writeEupsVersions(self, butler, clobber=False, doBackup=False):
        pass
