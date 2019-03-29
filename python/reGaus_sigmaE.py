"""
This script is for tests of the catalog sigma_e values.
"""

import galsim
import sim_analysis_utilities as utilities
from base_analysis import base_analysis
import numpy as np
import scipy.stats as stats
import math
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['ps.useafm'] = True
matplotlib.rcParams['pdf.use14corefonts'] = True
matplotlib.rcParams['text.usetex'] = True

default_n_bin_1d = 20

def do_ks_test(x,y,sample_names="two samples"):
    """Routine to do a 2-sample KS test and report on consistency of distributions"""

    (ks_d, ks_p) = stats.ks_2samp(x,y)
    if ks_p >= 0.05:
        print('Distributions for %s are consistent!'%sample_names)
    else:
        print('Distributions for %s differ!'%sample_names)
    print('(d, p)=(%f, %f)'%(ks_d, ks_p))

class sigmae_analysis(base_analysis):
    """
    Infer shape measurement error distribution from 90 degree rotated pairs; compare with the one
    that is implied from the catalog; and build up a model for the effective sigma_e as a function
    of galaxy parameters.
    """
    def __init__(self, max_snr=80., max_e=2., max_sigma_e=0.4,
                 param_cat_dir="/home/rmandelb/git/faint-data",
                 lower_psf_size=0., upper_psf_size=100.,
                 **kwargs):
        """
        Infer shape measurement error distribution (at least initially over all subfields) and
        compare with catalog information.

        Parameters
        ----------

        max_snr: float (optional)
            Maximum cmodel SNR to use when dividing up the sample and inferring sigma_e

        max_e: float (optional)
            Maximum |e| to use for the sample.

        max_sigma_e: float (optional)
            Maximum catalog sigma_e to use for the sample.

        param_cat_dir: string (optional)
            Directory containing the GalSim parametric catalog for the I<25.2 sample

        lower_psf_size: float (optional)
            Lower percentile to use for PSF size (i.e., use 50 to get everything above the median).

        upper_psf_size: float (optional)
            Upper percentile to use for PSF size (i.e., use 50 to get everything below the median).

        """
        super(sigmae_analysis, self).__init__(**kwargs)

        self.max_snr = max_snr
        self.max_e = max_e
        self.max_sigma_e = max_sigma_e
        self.param_cat_dir = param_cat_dir
        self.lower_psf_size = lower_psf_size
        self.upper_psf_size = upper_psf_size

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
        self.cat_est_sigma_e = np.zeros((default_n_bin_1d, default_n_bin_1d))
        self.est_sigma_e = np.zeros((default_n_bin_1d, default_n_bin_1d))
        self.higherorder_var_e = np.zeros((default_n_bin_1d, default_n_bin_1d))
        self.mean_res = np.zeros((default_n_bin_1d, default_n_bin_1d))
        self.mean_snr = np.zeros((default_n_bin_1d, default_n_bin_1d))
        self.true_rms = np.zeros((default_n_bin_1d, default_n_bin_1d))
        self.est_rms = np.zeros((default_n_bin_1d, default_n_bin_1d))
        self.cat_est_rms = np.zeros((default_n_bin_1d, default_n_bin_1d))
        self.corr_coeff = np.zeros((default_n_bin_1d, default_n_bin_1d))
        self.psf_size = []

        print('Beginning sigma_e analysis for the following scenario:')
        print('   %d bins in both SNR and resolution'%default_n_bin_1d)
        print('   Percentile ranges for PSF size: %.2f-%.2f'%\
            (self.lower_psf_size, self.upper_psf_size))
        print('')

    def select_data_sims(self):
        """
        Decide whether the analysis is going to use sims, data, or both.
        """
        self.use_data = False
        self.use_sims = True

    def process_data(self, data, fieldname):
        """
        This analysis does not apply to real data, do nothing
        """
        pass

    def process_sims(self, data, fieldname):
        """
        Get estimated and true shape information, lensing shear, and other information.
        """
        # Define a mask with some set of cuts.  Two notes:
        # 1) cancel_shape_noise MUST be true for this test (so we don't use
        #    self.cancel_shape_noise).  It's not an option that the user
        #    can change, because the 90 degree rotated pairs are required.
        # 2) bright_star_mask is set to False (instead of
        #    self.bright_star_mask) because we are running on sims.
        mask = utilities.galaxy_selector(catalog=data,
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
        if self.reweight_s16a:
            self.weight.append(data[mask]['reweight_s16a'])
        else:
            self.weight.append(np.ones_like(true_e1).astype(float))
        self.psf_size.append(utilities.get_psf_size(data[mask]))
        # Store subfield true shears.
        self.g1_true.append(g1_true)
        self.g2_true.append(g2_true)

    def finalise(self):
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
        print('Overall N, <w>, RMS: ',len(weight_input), np.mean(weight_input), rms)

        # Get implied RMS from catalog.
        input_rms = np.sqrt(np.sum(weight_input*expected_sigma**2)/np.sum(weight_input))
        print('Implied from catalog: ',input_rms)

        # Now set up sliding window definition
        # Only use the first in pair
        snr_lower, snr_upper, res_lower, res_upper = \
            utilities.sliding_window_def(master_snr[0::2][cond_mask],
                                         master_resolution[0::2][cond_mask],
                                         n_bin_1d=default_n_bin_1d)

        # Now loop over SNR and resolution bins, using the first
        # galaxy in a pair to stick them in a bin.
        for i in range(len(snr_lower)):
            for j in range(len(snr_lower)):
                print('Beginning bin: SNR in [%.3f, %.3f], resolution in [%.3f, %.3f]'%\
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
                print('Using %d pairs'%bin_mask.astype(int).sum())
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
                print('N, <w>, RMS in bin: ',len(weight_input), np.mean(weight_input), bin_rms)

                # Get implied RMS from catalog.
                expected_sigma = np.concatenate((master_sigma_e[0::2][bin_mask],
                                                 master_sigma_e[1::2][bin_mask]))
                input_bin_rms = np.sqrt(
                    np.sum(weight_input*expected_sigma**2)/np.sum(weight_input))
                print('RMS in bin implied from catalog: ',input_bin_rms)

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
                print('RMS values and responsivities per component in this bin:')
                print('  observed, including shape measurement error: %.4f, %.4f'%\
                    (rms_observed, 2.*(1.-rms_observed**2)))
                print('  estimated using distribution to remove measurement error: %.4f, %.4f'%\
                    (est_true_dist, 2.*(1.-est_true_dist**2)))

        print('')
        print('')
        x_vals = np.vstack([np.log(self.mean_snr/20.).flatten(),
                            np.log(self.mean_res/0.5).flatten(),
                            np.ones(len(self.mean_res.flatten()))]).transpose()
        est_log_sig_e = np.log(self.est_sigma_e.flatten())
        result  = np.linalg.lstsq(x_vals, est_log_sig_e)
        #Eq.(18)
        model_sigma_e = np.exp(result[0][2])*((self.mean_snr/20.)**result[0][0])*\
            ((self.mean_res/0.5)**result[0][1])
        print('model sigma_e = %.3f (SNR/20)^%.3f (res/0.5)^%.3f'%\
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
        plt.savefig(self.plot_pref+'sigmae_estimated.'+self.plot_type)


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
        plt.savefig(self.plot_pref+'sigmae_cat.'+self.plot_type)

        # Plot ratio of catalog to estimated sigma_e.
        fig = plt.figure()
        ax = fig.add_subplot(111)
        sc = ax.scatter(self.mean_snr, self.mean_res, c=self.cat_est_sigma_e/self.est_sigma_e, cmap=cm)
        plt.colorbar(sc)
        ax.set_title('Catalog/estimated RMS measurement error')
        plt.xscale('log')
        plt.xlabel('SNR')
        plt.ylabel('resolution')
        plt.savefig(self.plot_pref+'sigmae_cat_estimated_ratio.'+self.plot_type)


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
        plt.savefig(self.plot_pref+'sigmae_model.'+self.plot_type)
        # model / estimated sigma_e as a function of SNR, resolution.  Should be close to 1.
        fig = plt.figure()
        ax = fig.add_subplot(111)
        sc = ax.scatter(self.mean_snr, self.mean_res, c=self.est_sigma_e/model_sigma_e, cmap=cm)
        plt.colorbar(sc)
        ax.set_title('Estimated/model RMS')
        plt.xscale('log')
        plt.xlabel('SNR')
        plt.ylabel('resolution')
        plt.savefig(self.plot_pref+'sigmae_ratio.'+self.plot_type)

        # save data table for model/actual
        out_data = (self.mean_snr.flatten(), self.mean_res.flatten(),
                    (self.est_sigma_e/model_sigma_e).flatten())
        out_data = np.column_stack(out_data)
        np.savetxt(self.plot_pref+'sigmae_ratio.dat', out_data)
