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

def m_func(x, b, c, d, e):
    """
    Empirically-motivated model we are trying to fit for m(SNR, res).

    x[0,:] is SNR, x[1,:] is resolution.
    """
    model = (x[0,:]/20.)**d
    model *= b*((x[1,:]/0.5)**c)
    return model+e


def a_func(x, b, c, d):
    """
    Empirically-motivated model we are trying to fit for a(SNR, res).

    x[0,:] is SNR, x[1,:] is resolution.
    """
    model = b*(x[1,:]-c)
    model *= (x[0,:]/20.)**d
    return model

class reGaus_shearCalibConfig(pexConfig.Config):
    anaDir      =   pexConfig.Field(dtype=str, default='anaCat', doc = 'simulation catalogs')
    obsDir      =   pexConfig.Field(dtype=str, default='/work/xiangchong.li/work/S16AFPFS/final/noBg/field/', doc = 'observed catalogs')
    pltDir      =   pexConfig.Field(dtype=str, default='plot/shearCalib/', doc = 'plot directory')
    nbin_1d     =   pexConfig.Field(dtype=int, default=8, doc = 'number of bins for 1-D parameter space')
    ngal_bin    =   pexConfig.Field(dtype=int, default=150000, doc = 'number of bins for 1-D parameter space')
    use_model   =   pexConfig.Field(dtype=bool, default=True, doc = 'whether to use model to estimate weight')
    cancel_shape_noise=  pexConfig.Field(dtype=bool, default=True, doc = 'whether to use orthogonal galaxies to cancel shape noise')
    separate_components= pexConfig.Field(dtype=bool, default=False, doc = 'whether to seperate g1 g2 components')
    round_trip  =   pexConfig.Field(dtype=bool, default=False, doc = 'do a round trip?')
    test_weight_bias  =   pexConfig.Field(dtype=bool, default=False, doc = 'test weight bias?')
    bin_plot    =   pexConfig.Field(dtype=bool, default=False, doc = 'whether to seperate g1 g2 components')
    
    max_snr =   pexConfig.Field(dtype=float, default=100000., doc = 'Maximum value of signal to noise ratio')
    max_e   =   pexConfig.Field(dtype=float, default=2., doc = 'Maximum value of distortion |e|')
    max_sigma_e =   pexConfig.Field(dtype=float, default=0.4, doc = 'Maximum value of sigmaE')
    lower_psf_size  =   pexConfig.Field(dtype=float, default=0., doc = 'Minimum of PSF size')
    upper_psf_size  =   pexConfig.Field(dtype=float, default=100., doc = 'Maximum of PSF size')
    
    def setDefaults(self):
        pexConfig.Config.setDefaults(self)

    def validate(self):
        pexConfig.Config.validate(self)
        if self.config.test_weight_bias and (not self.config.round_trip):
            raise RuntimeError("Cannot test weight bias when not doing round trip")
        if self.config.test_weight_bias and self.config.cancel_shape_noise:
            raise RuntimeError("Cannot test weight bias when cancelling shape noise")
        # Checks that the paths provided are valid
        if not os.path.isdir(self.anaDir):
            raise RuntimeError("Simulation catalog directory %s does not exist "
                         "or is not a directory" % (self.anaDir))
        if not os.path.isdir(self.obsDir):
            raise RuntimeError("Data catalog directory %s does not exist "
                         "or is not a directory" % (self.anaDir))
        if not os.path.exists(self.pltDir):
            os.mkdir(self.pltDir)

class reGaus_shearCalibRunner(TaskRunner):
    @staticmethod
    def getTargetList(parsedCmd, **kwargs):
        return [(ref, kwargs) for ref in range(1)] 

def unpickle(factory, args, kwargs):
    """Unpickle something by calling a factory"""
    return factory(*args, **kwargs)

class reGaus_shearCalibTask(BatchPoolTask):
    ConfigClass     =   reGaus_shearCalibConfig
    RunnerClass     =   reGaus_shearCalibRunner
    _DefaultName    =   "reGaus_shearCalib"
    
    def __init__(self,**kwargs):
        BatchPoolTask.__init__(self, **kwargs)
        
        if self.config.cancel_shape_noise and not self.config.round_trip:
            self.prefix='shear'
        if self.config.cancel_shape_noise and self.config.round_trip:
            self.prefix='test_shear'
        if not self.config.cancel_shape_noise and self.config.round_trip:
            self.prefix='weightBias'
        
        # Initialize array to store multiplicative bias and associated quantities.
        self.snr_bins = np.zeros((self.config.nbin_1d, self.config.nbin_1d))
        self.res_bins = np.zeros((self.config.nbin_1d, self.config.nbin_1d))
        self.m1_bins = np.zeros((self.config.nbin_1d, self.config.nbin_1d))
        self.a1_bins = np.zeros((self.config.nbin_1d, self.config.nbin_1d))
        self.m2_bins = np.zeros((self.config.nbin_1d, self.config.nbin_1d))
        self.a2_bins = np.zeros((self.config.nbin_1d, self.config.nbin_1d))
        self.sigm1_bins = np.zeros((self.config.nbin_1d, self.config.nbin_1d))
        self.siga1_bins = np.zeros((self.config.nbin_1d, self.config.nbin_1d))
        
        #True g1,g2
        self.g1_true = []
        self.g2_true = []
        self.psf_e1 = []
        self.psf_e2 = []
        self.sigma_g = []
        self.weight_s16a = []

        self.add_m_val  =   0.
        self.add_a_val  =   0.
        
        self.log.info('Beginning shearCalib analysis for the following scenario:')
        self.log.info('   %d bins in both SNR and resolution'%self.config.nbin_1d)
        self.log.info('   Percentile ranges for PSF size: %.2f-%.2f'%\
            (self.config.lower_psf_size, self.config.upper_psf_size))
        self.log.info('')
        
    def __reduce__(self):
        """Pickler"""
        return unpickle, (self.__class__, [], dict(config=self.config, name=self._name,
                parentTask=self._parentTask, log=self.log))
    
    @abortOnError
    def run(self,Id):
        #Prepare the pool
        rows    =   Table.read('noiVar_rPsf_s18a.csv')
        pool    =   Pool("reGaus_shearCalib")
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
                                         max_snr=self.config.max_snr,
                                         cancel_shape_noise=self.config.cancel_shape_noise,
                                         max_e=self.config.max_e,
                                         max_sigma_e=self.config.max_sigma_e)
                                         
        return data[mask]
        
    def finalise(self,catalog,reweight_s16a=True):
        """
        Compute shear biases and make plots
        """
        # Get mask for imaging conditions.  Note this will always preserve shape noise cancellation
        # due to how the sims are set up, though we can do an assert just to check.
        psf_size = utilities.get_psf_size(catalog)
        cond_mask = np.ones_like(catalog).astype(bool)
        if self.config.lower_psf_size > 0.:
            cut_val = np.percentile(psf_size, self.config.lower_psf_size)
            cond_mask &= (psf_size >= cut_val)
        if self.config.upper_psf_size < 100.:
            cut_val = np.percentile(psf_size, self.config.upper_psf_size)
            cond_mask &= (psf_size < cut_val)
        assert np.all(cond_mask[0::2]==cond_mask[1::2]),\
                "Mask for imaging conditions does not obey shape noise cancellation sanity check"
        catalog = catalog[cond_mask]
        
        # Set up sliding window definition and make a mask for each SNR, resolution bin.
        snr = utilities.get_snr(catalog)
        res = utilities.get_res(catalog)
        
        if self.config.round_trip:
            if self.config.test_weight_bias:
                weight_bias = True
            else:
                weight_bias = False
            model_m = utilities.get_m_model(catalog,weight_bias=weight_bias)
            model_a = utilities.get_a_model(catalog, data_version=weight_bias=weight_bias)
            psf_e1, psf_e2 = utilities.get_psf_ellip(catalog)
            model_c1 = model_a*psf_e1
            model_c2 = model_a*psf_e2
        
        
        # We need to get the weights here, to use when averaging over galaxies to get mean effective
        # SNR and resolution.
        # Note, it might be worth making a utility for this combined shape + S16A weight, since it's
        # done here and in estimate_subfield_shear.
        if self.config.use_model: 
            shape_weights = utilities.get_weight_model(catalog)
        else:
            shape_weights = catalog['ishape_hsm_regauss_derived_shape_weight']
        if reweight_s16a:
            shape_weights = shape_weights*catalog['reweight_s16a']
        if self.config.cancel_shape_noise:
            shape_weights[1::2] = shape_weights[0::2]
            assert np.all(catalog['reweight_s16a'][0::2]==catalog['reweight_s16a'][1::2]), \
                'Galaxies in pairs do not have the same reweight_s16a factors'
        
        # Sort out binning. 
        snr_lower, snr_upper, res_lower, res_upper = \
            utilities.sliding_window_def(snr, res,n_gal=self.config.ngal_bin,
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
            
            # Get the estimate in bins in SNR/res as well.
            for j in range(self.config.nbin_1d):
                for k in range(self.config.nbin_1d):
                    mm = bin_mask['%d_%d'%(j,k)] & (subfield_id == i)
                    # Some weird bins might have no objects in particular subfields.
                    # Don't use estimate_subfield_shear on those!
                    if mm.astype(int).sum() > 2:
                        # Another criterion: make sure the objects have nonzero weight
                        if np.sum(shape_weights[mm])>0:
                            if self.self.config.round_trip:
                                roundtrip_data = (model_m[mm], model_c1[mm], model_c2[mm])
                            else:
                                roundtrip_data = None
                            tmp_g1, tmp_g2, tmp_resp = utilities.estimate_subfield_shear(
                                catalog[mm], use_model=self.config.use_model,
                                cancel_shape_noise=self.config.cancel_shape_noise,
                                simple_resp=True, report_resp=True,
                                reweight_s16a=reweight_s16a,
                                roundtrip_data=roundtrip_data,
                                test_weight_bias=self.config.test_weight_bias)
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

        self.weight_s16a = np.array(self.weight_s16a)
        self.g1_true = np.array(self.g1_true)
        self.g2_true = np.array(self.g2_true)
        self.psf_e1 = np.array(self.psf_e1)
        self.psf_e2 = np.array(self.psf_e2)

        #For each of the sliding bins in SNR/resolution.
        for j in range(self.config.nbin_1d):
            for k in range(self.config.nbin_1d):
                sigma_vals = self.sigma_g_bins[j,k,:]
                if reweight_s16a:
                    division_val = np.sqrt(self.weight_s16a)
                    division_val[division_val<1.e-6] = 1.e-6
                    sigma_vals /= division_val
                tmp_use = sigma_vals < 1.e10
                if self.config.separate_components:
                    self.m1_bins[j,k], self.a1_bins[j,k], \
                        self.sigm1_bins[j,k], self.siga1_bins[j,k], \
                        self.m2_bins[j,k], self.a2_bins[j,k], sig_m2, sig_a2, \
                        fig1, fig2, fig3, fig4 = \
                        utilities.get_m_c_a(self.g1_bins[j,k,tmp_use],
                                            self.g2_bins[j,k,tmp_use],
                                            self.g1_true[tmp_use],
                                            self.g2_true[tmp_use],
                                            self.psf_e1[tmp_use],
                                            self.psf_e2[tmp_use],
                                            sigma_vals[tmp_use],
                                            separate_components=self.config.separate_components)
                else:
                    self.m1_bins[j,k], self.a1_bins[j,k], \
                        self.sigm1_bins[j,k], self.siga1_bins[j,k], \
                        fig1, fig2 = \
                        utilities.get_m_c_a(self.g1_bins[j,k,tmp_use],
                                            self.g2_bins[j,k,tmp_use],
                                            self.g1_true[tmp_use], self.g2_true[tmp_use],
                                            self.psf_e1[tmp_use], self.psf_e2[tmp_use],
                                            sigma_vals[tmp_use],
                                            separate_components=self.config.separate_components)

                self.log.info('Bin %d %d bias analysis:'%(j,k))
                if self.config.separate_components:
                    self.log.info('m1 = %.3f +/- %.3f' % (self.m1_bins[j,k], self.sigm1_bins[j,k]))
                    self.log.info('a1 = %.3f +/- %.3f' % (self.a1_bins[j,k], self.siga1_bins[j,k]))
                    self.log.info('m2 = %.3f +/- %.3f' % (self.m2_bins[j,k], sig_m2))
                    self.log.info('a2 = %.3f +/- %.3f' % (self.a2_bins[j,k], sig_a2))
                    if self.config.bin_plot:
                        fig1.savefig(self.config.pltDir+\
                                         '%s_bias_bin%d_%d_1.png'%(j,k))
                        fig2.savefig(self.config.pltDir+\
                                         '%s_bias_bin%d_%d_2.png'%(j,k))
                        fig3.savefig(self.config.pltDir+\
                                         '%s_additive_residual_bin%d_%d_1.png'%(j,k) )
                        fig4.savefig(self.config.pltDir+\
                                         '%s_additive_residual_bin%d_%d_2.png'%(j,k) )
                else:
                    self.log.info('m = %.3f +/- %.3f' % (self.m1_bins[j,k], self.sigm1_bins[j,k]))
                    self.log.info('a = %.3f +/- %.3f' % (self.a1_bins[j,k], self.siga1_bins[j,k]))
                    if self.config.bin_plot:
                        fig1.savefig(self.config.pltDir+\
                                         '%s_bias_bin%d_%d.png'%(j,k))
                        fig2.savefig(self.config.pltDir+\
                                         '%s_additive_residual_bin%d_%d.png'%(j,k) )
        # Do some clever summary plots and model fitting.
        cm = plt.cm.get_cmap('RdYlBu')
        if self.config.separate_components:
            m_vmin = min(np.min(self.m1_bins+self.add_m_val), np.min(self.m2_bins+self.add_m_val))
            m_vmax = max(np.max(self.m1_bins+self.add_m_val), np.max(self.m2_bins+self.add_m_val))
        else:
            m_vmin = np.min(self.m1_bins+self.add_m_val)
            m_vmax = np.max(self.m1_bins+self.add_m_val)
        fig = plt.figure()
        ax = fig.add_subplot(111)
        sc = ax.scatter(self.snr_bins, self.res_bins, c=self.m1_bins+self.add_m_val, cmap=cm,
                    vmin=m_vmin, vmax=m_vmax, s=80, edgecolors='none')
        cbar = plt.colorbar(sc)
        if self.config.separate_components:
            cbar.ax.set_title(r'$m_1$')
        else:
            if self.add_m_val > 0:
                cbar.ax.set_title(r'$m+$constant')
            else:
                cbar.ax.set_title(r'$m$')
        plt.xscale('log')
        plt.xlabel('S/N')
        plt.ylabel('resolution')
        if self.config.separate_components:
            plt.savefig(self.config.pltDir+'%s_m1_2d.eps'%(self.prefix))
        else:
            plt.savefig(self.config.pltDir+'%s_m_2d.eps'%(self.prefix))

            # Also plot errors in this case.
            fig = plt.figure()
            ax = fig.add_subplot(111)
            sc = ax.scatter(self.snr_bins, self.res_bins, c=self.sigm1_bins, cmap=cm, s=80, edgecolors='none')
            cbar = plt.colorbar(sc)
            cbar.ax.set_title(r'$\sigma_m$')
            plt.xscale('log')
            plt.xlabel('S/N')
            plt.ylabel('resolution')
            plt.savefig(self.config.pltDir+'%s_sigm_2d.eps'%(self.prefix))

        if self.config.separate_components:
            fig = plt.figure()
            ax = fig.add_subplot(111)
            sc = ax.scatter(self.snr_bins, self.res_bins, c=self.m2_bins+self.add_m_val, cmap=cm,
                            vmin=m_vmin, vmax=m_vmax, s=80, edgecolors='none')
            cbar = plt.colorbar(sc)
            cbar.ax.set_title(r'$m_2$')
            plt.xscale('log')
            plt.xlabel('S/N')
            plt.ylabel('resolution')
            plt.savefig(self.config.pltDir+'%s_m2_2d.eps'%(self.prefix))
        
        if self.config.separate_components:
            a_vmin = min(np.min(self.a1_bins+self.add_a_val), np.min(self.a2_bins+self.add_a_val))
            a_vmax = max(np.max(self.a1_bins+self.add_a_val), np.max(self.a2_bins+self.add_a_val))
        else:
            a_vmin = np.min(self.a1_bins+self.add_a_val)
            a_vmax = np.max(self.a1_bins+self.add_a_val)
        fig = plt.figure()
        ax = fig.add_subplot(111)
        sc = ax.scatter(self.snr_bins, self.res_bins, c=self.a1_bins+self.add_a_val, cmap=cm,
                        vmin=a_vmin, vmax=a_vmax, s=80, edgecolors='none')
        cbar = plt.colorbar(sc)
        if self.config.separate_components:
            cbar.ax.set_title(r'$a_1$')
        else:
            if self.add_m_val > 0:
                cbar.ax.set_title(r'$a$')
            else:
                cbar.ax.set_title(r'$a$')
        plt.xscale('log')
        plt.xlabel('S/N')
        plt.ylabel('resolution')
        if self.config.separate_components:
            plt.savefig(self.config.pltDir+'%s_a1_2d.eps'%(self.prefix))
        else:
            plt.savefig(self.config.pltDir+'%s_a_2d.eps'%(self.prefix))

            # Also plot errors in this case.
            fig = plt.figure()
            ax = fig.add_subplot(111)
            sc = ax.scatter(self.snr_bins, self.res_bins, c=self.siga1_bins, cmap=cm, s=80, edgecolors='none')
            cbar = plt.colorbar(sc)
            cbar.ax.set_title(r'$\sigma_a$')
            plt.xscale('log')
            plt.xlabel('S/N')
            plt.ylabel('resolution')
            plt.savefig(self.config.pltDir+'%s_siga_2d.eps'%(self.prefix))

        if self.config.separate_components:
            fig = plt.figure()
            ax = fig.add_subplot(111)
            sc = ax.scatter(self.snr_bins, self.res_bins, c=self.a2_bins+self.add_a_val, cmap=cm,
                            vmin=a_vmin, vmax=a_vmax, s=80, edgecolors='none')
            cbar = plt.colorbar(sc)
            cbar.ax.set_title(r'$a_2$')
            plt.xscale('log')
            plt.xlabel('S/N')
            plt.ylabel('resolution')
            plt.savefig(self.config.pltDir+'%s_a2_2d.eps'%(self.prefix))

        if self.config.separate_components:
            return

        # Also do model fitting in this case.
        out_data = (self.snr_bins.flatten(), self.res_bins.flatten(),
                    self.m1_bins.flatten(), self.sigm1_bins.flatten(),
                    self.a1_bins.flatten(), self.siga1_bins.flatten())
        out_data = np.column_stack(out_data)
        np.savetxt(self.config.pltDir+'%s_m_a_2d.dat'%(self.prefix), out_data)

        if self.config.cancel_shape_noise:
            fake_x = np.vstack((self.snr_bins.flatten(), self.res_bins.flatten()))
            m_opt, m_cov = opt.curve_fit(
                m_func, fake_x, self.m1_bins.flatten(), sigma=self.sigm1_bins.flatten(),
                bounds=([-np.inf,0,-np.inf,-np.inf],[np.inf,np.inf,0,np.inf]))
        else:
            fake_x = np.vstack((self.snr_bins[1:-1,1:-1].flatten(),
                                self.res_bins[1:-1,1:-1].flatten()))
            m_opt, m_cov = opt.curve_fit(
                m_func, fake_x, self.m1_bins[1:-1,1:-1].flatten(),
                sigma=self.sigm1_bins[1:-1,1:-1].flatten(),
                bounds=([-np.inf,0,-np.inf,-np.inf],[np.inf,np.inf,0,np.inf]))
            fake_x = np.vstack((self.snr_bins.flatten(), self.res_bins.flatten()))
        self.log.info('m model:')
        self.log.info('m=b*(res/0.5)^c*(snr/20)^d+e')
        for ind in range(len(m_opt)):
            self.log.info('%f +/- %f'%(m_opt[ind],m_cov[ind,ind])) 
        m_model = m_func(fake_x, *m_opt).reshape(self.config.nbin_1d, self.config.nbin_1d)
        # And plot the model, residuals, and their significance
        fig = plt.figure()
        ax = fig.add_subplot(111)
        sc = ax.scatter(self.snr_bins, self.res_bins, c=m_model, cmap=cm,
                        vmin=m_vmin, vmax=m_vmax, s=80, edgecolors='none')
        cbar = plt.colorbar(sc)
        cbar.ax.set_title(r'$m$ (model)')
        plt.xscale('log')
        plt.xlabel('S/N')
        plt.ylabel('resolution')
        plt.savefig(self.config.pltDir+'%s_m_model_2d.eps'%(self.prefix))

        fig = plt.figure()
        ax = fig.add_subplot(111)
        sc = ax.scatter(self.snr_bins, self.res_bins, c=self.m1_bins-m_model,
                        cmap=cm, s=80, edgecolors='none')
        cbar = plt.colorbar(sc)
        cbar.ax.set_title(r'$m$ - model')
        plt.xscale('log')
        plt.xlabel('S/N')
        plt.ylabel('resolution')
        plt.savefig(self.config.pltDir+'%s_Delta_m_model_2d.eps'%(self.prefix))

        fig = plt.figure()
        ax = fig.add_subplot(111)
        sc = ax.scatter(self.snr_bins, self.res_bins,
                        c=(self.m1_bins-m_model)/self.sigm1_bins, cmap=cm, s=80, edgecolors='none')
        cbar = plt.colorbar(sc)
        cbar.ax.set_title(r'($m$ - model)/$\sigma_m$')
        plt.xscale('log')
        plt.xlabel('S/N')
        plt.ylabel('resolution')
        plt.savefig(self.config.pltDir+'%s_m_model_signif_2d.eps'%(self.prefix))

        #For additive bias
        if self.config.cancel_shape_noise:
            fake_x = np.vstack((self.snr_bins.flatten(), self.res_bins.flatten()))
            a_opt, a_cov = opt.curve_fit(
                a_func, fake_x, self.a1_bins.flatten(), sigma=self.siga1_bins.flatten(),
                bounds=([-np.inf,-np.inf,-np.inf],[np.inf,np.inf,0.]))
        else:
            fake_x = np.vstack((self.snr_bins[1:-1,1:-1].flatten(),
                                self.res_bins[1:-1,1:-1].flatten()))
            a_opt, a_cov = opt.curve_fit(
                a_func, fake_x, self.a1_bins[1:-1,1:-1].flatten(),
                sigma=self.siga1_bins[1:-1,1:-1].flatten(),
                bounds=([-np.inf,-np.inf,-np.inf],[np.inf,np.inf,0.]))
            fake_x = np.vstack((self.snr_bins.flatten(), self.res_bins.flatten()))
        self.log.info('a model:')
        self.log.info('a=b*(res-c)*(snr/20)^d')
        for ind in range(len(a_opt)):
            self.log.info('%f +/- %f'%(a_opt[ind],a_cov[ind,ind]))
        a_model = a_func(fake_x, *a_opt).reshape(self.config.nbin_1d, self.config.nbin_1d)
        # Plot model, residuals, significance
        fig = plt.figure()
        ax = fig.add_subplot(111)
        sc = ax.scatter(self.snr_bins, self.res_bins, c=a_model, cmap=cm,
                        vmin=a_vmin, vmax=a_vmax, s=80, edgecolors='none')
        cbar = plt.colorbar(sc)
        cbar.ax.set_title(r'$a$ (model)')
        plt.xscale('log')
        plt.xlabel('S/N')
        plt.ylabel('resolution')
        plt.savefig(self.config.pltDir+'%s_a_model_2d.eps'%(self.prefix))

        fig = plt.figure()
        ax = fig.add_subplot(111)
        sc = ax.scatter(self.snr_bins, self.res_bins,
                        c=self.a1_bins-a_model, cmap=cm, s=80, edgecolors='none')
        cbar = plt.colorbar(sc)
        cbar.ax.set_title(r'$a$ - model')
        plt.xscale('log')
        plt.xlabel('S/N')
        plt.ylabel('resolution')
        plt.savefig(self.config.pltDir+'%s_Delta_a_model_2d.eps'%(self.prefix))

        fig = plt.figure()
        ax = fig.add_subplot(111)
        sc = ax.scatter(self.snr_bins, self.res_bins,
                        c=(self.a1_bins-a_model)/self.siga1_bins, cmap=cm, s=80, edgecolors='none')
        cbar = plt.colorbar(sc)
        cbar.ax.set_title(r'($a$ - model)/$\sigma_a$')
        plt.xscale('log')
        plt.xlabel('S/N')
        plt.ylabel('resolution')
        plt.savefig(self.config.pltDir+'%s_a_model_signif_2d.eps'%(self.prefix))
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
