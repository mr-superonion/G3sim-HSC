#!/usr/bin/env python
import os
import sys
import matplotlib
matplotlib.use('agg')
matplotlib.rcParams['ps.useafm'] = True
matplotlib.rcParams['pdf.use14corefonts'] = True
matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt
import sim_analysis_utilities as utilities

import numpy as np
import scipy.stats as stats
from astropy.table import Table, vstack

import lsst.daf.base as dafBase
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase

from lsst.pipe.base import TaskRunner
from lsst.ctrl.pool.parallel import BatchPoolTask
from lsst.ctrl.pool.pool import Pool, abortOnError

class reGaus_eRMSConfig(pexConfig.Config):
    anaDir  =   pexConfig.Field(dtype=str, default='anaCat', doc = 'simulation catalogs')
    obsDir  =   pexConfig.Field(dtype=str, default='../database/s18/fields', doc = 'observed catalogs')
    pltDir  =   pexConfig.Field(dtype=str, default='plot/eRMS/', doc = 'plot directory')
    nbin_1d =   pexConfig.Field(dtype=int, default=20, doc = 'number of bins for 1-D parameter space')
    max_snr =   pexConfig.Field(dtype=float, default=200., doc = 'Maximum value of signal to noise ratio')
    max_e   =   pexConfig.Field(dtype=float, default=2., doc = 'Maximum value of distortion |e|')
    max_sigma_e =   pexConfig.Field(dtype=float, default=0.4, doc = 'Maximum value of eRMS')
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

class reGaus_eRMSRunner(TaskRunner):
    @staticmethod
    def getTargetList(parsedCmd, **kwargs):
        return [(ref, kwargs) for ref in range(1)]
def unpickle(factory, args, kwargs):
    """Unpickle something by calling a factory"""
    return factory(*args, **kwargs)
class reGaus_eRMSTask(BatchPoolTask):
    ConfigClass     =   reGaus_eRMSConfig
    RunnerClass     =   reGaus_eRMSRunner
    _DefaultName    =   "reGaus_eRMS"

    def __init__(self,**kwargs):
        BatchPoolTask.__init__(self, **kwargs)
        self.max_snr = self.config.max_snr
        self.max_e = self.config.max_e
        self.max_sigma_e = self.config.max_sigma_e
        self.lower_psf_size = self.config.lower_psf_size
        self.upper_psf_size = self.config.upper_psf_size

        # Initialise lists to store data
        self.names  =   ('e1','e2','sigma_e','model_sigma_e','model_cc','resolution','snr','mag','weight')

    def __reduce__(self):
        """Pickler"""
        return unpickle, (self.__class__, [], dict(config=self.config, name=self._name,
                parentTask=self._parentTask, log=self.log))

    @abortOnError
    def run(self,Id):
        #Prepare the pool
        pool    =   Pool("reGaus_eRMS")
        pool.cacheClear()
        outcomes=   pool.map(self.process,utilities.field_keys)
        outcomes=   [x for x in outcomes if x is not None]
        outtable=   vstack(outcomes)
        self.finalise(outtable)
        return

    def process(self,cache,fieldKey,reweight_s16a=True,use_wl_flag=True):
        fieldname=  utilities.field_keys[fieldKey]
        if fieldname==  'AEGIS':
            return None
        self.log.info('Processing field %s' %(fieldname))
        fname   =   os.path.join(self.config.anaDir,fieldname,'%s.fits'%(fieldname))
        data    =   Table.read(fname)
        """
        Get estimated and true shape information, lensing shear, and other information.
        """
        # Define a mask with some set of cuts.  Two notes:
        # cancel_shape_noise MUST be false since this is observed
        # data
        mask    =   utilities.galaxy_selector(catalog=data,
                                         use_wl_flag=use_wl_flag,
                                         bright_star_mask=False,
                                         max_snr=self.max_snr,
                                         force_90_pair=False)
        dataTab =   []
        # Store NumPy arrays for e1,e2,sigma_e,model_sigma_e,
        # model_cc,resolution,snr,weight
        dataTab.append(data[mask]['ishape_hsm_regauss_e1'])
        dataTab.append(data[mask]['ishape_hsm_regauss_e2'])
        dataTab.append(utilities.get_sigma_e(data[mask]))
        dataTab.append(utilities.get_sigma_e_model(data[mask]))
        dataTab.append(np.zeros_like(data[mask]).astype(float))
        dataTab.append(data[mask]['ishape_hsm_regauss_resolution'])
        dataTab.append(utilities.get_snr(data[mask]))
        dataTab.append(utilities.get_imag(data[mask]))
        if reweight_s16a:
            dataTab.append(data[mask]['reweight_s16a'])
        else:
            dataTab.append(np.ones_like(dataTab[0]).astype(float))
        # Store subfield true shears.
        return Table(data=dataTab,names=self.names)

    def finalise(self,inTable):
        """
        Make distributions of the quantities that we need to calculate
        in order to estimate the distribution of shape measurement
        errors.  Everything here is unweighted; the only exception is
        weighting to account for the S16A conditions.
        """
        # First: assemble the arrays for individual subfields into one big master catalog.
        master_e1 = inTable['e1']
        master_e2 = inTable['e2']
        master_sigma_e = inTable['sigma_e']
        master_model_sigma_e = inTable['model_sigma_e']
        master_resolution = inTable['resolution']
        master_snr  = inTable['snr']
        master_mag  =   inTable['mag']
        master_weight = inTable['weight']

        # Start by estimating the RMS as a function of mag, resolution, and SNR in 1d.  Do this
        # using catalog and model sigma_e's.
        quantities = ['resolution', 'snr', 'mag']
        for x in range(len(quantities)):
            vals = eval('master_'+quantities[x])
            bin_edges = np.percentile(vals, np.linspace(0., 100., self.config.nbin_1d+1))
            bin_lower = bin_edges[:-1]
            bin_upper = bin_edges[1:]
            cat_rms = np.zeros(self.config.nbin_1d)
            model_rms = np.zeros(self.config.nbin_1d)
            cat_rms_corr = np.zeros(self.config.nbin_1d)
            model_rms_corr = np.zeros(self.config.nbin_1d)
            mean_val = np.zeros(self.config.nbin_1d)
            for ind in range(self.config.nbin_1d):
                in_bin = (vals >= bin_lower[ind]) & (vals < bin_upper[ind])
                # Calculate a standard per-component RMS using the catalog and model measurement error values.
                cat_rms[ind] = np.sqrt(0.5*np.sum(
                        master_weight[in_bin]*\
                            (master_e1[in_bin]**2+master_e2[in_bin]**2-2*master_sigma_e[in_bin]**2)) / \
                                      np.sum(master_weight[in_bin]))
                model_rms[ind] = np.sqrt(0.5*np.sum(
                        master_weight[in_bin]*\
                            (master_e1[in_bin]**2+master_e2[in_bin]**2-2*master_model_sigma_e[in_bin]**2)) / \
                                        np.sum(master_weight[in_bin]))
                # Also do a version that includes correction for the effect of correlations between shape and measurement error.
                mean_val[ind] = \
                    np.sum(master_weight[in_bin]*vals[in_bin])/np.sum(master_weight[in_bin])
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.plot(mean_val, cat_rms, 'r-', label='With catalog sigmas')
            ax.plot(mean_val, model_rms, 'k-', label='With model sigmas')
            plt.legend()
            if 'snr' in quantities[x]: plt.xscale('log')
            plt.xlabel(quantities[x])
            plt.ylabel(r'$\sigma_{SN}$')
            plt.savefig(os.path.join(self.config.pltDir,'intrinsicshape_'+quantities[x]+'.eps'))

        # Then do the 2d estimation as a function of SNR and resolution.
        snr_lower, snr_upper, res_lower, res_upper = \
            utilities.sliding_window_def(master_snr, master_resolution,
                                         n_bin_1d=self.config.nbin_1d)
        cat_rms = np.zeros((self.config.nbin_1d, self.config.nbin_1d))
        model_rms = np.zeros((self.config.nbin_1d, self.config.nbin_1d))
        mean_snr = np.zeros((self.config.nbin_1d, self.config.nbin_1d))
        mean_res = np.zeros((self.config.nbin_1d, self.config.nbin_1d))
        for i in range(len(snr_lower)):
            for j in range(len(snr_lower)):
                in_bin = \
                    (master_snr >= snr_lower[i]) & \
                    (master_snr < snr_upper[i]) & \
                    (master_resolution >= res_lower[i,j]) & \
                    (master_resolution < res_upper[i,j])
                cat_rms[i,j] = np.sqrt(0.5*np.sum(
                        master_weight[in_bin]*\
                            (master_e1[in_bin]**2+master_e2[in_bin]**2-2*master_sigma_e[in_bin]**2)) / \
                                           np.sum(master_weight[in_bin]))
                model_rms[i,j] = np.sqrt(0.5*np.sum(
                        master_weight[in_bin]*\
                            (master_e1[in_bin]**2+master_e2[in_bin]**2-2*master_model_sigma_e[in_bin]**2)) / \
                                             np.sum(master_weight[in_bin]))
                mean_snr[i,j] = \
                    np.sum(master_weight[in_bin]*master_snr[in_bin])/np.sum(master_weight[in_bin])
                mean_res[i,j] = \
                    np.sum(master_weight[in_bin]*master_resolution[in_bin])/np.sum(master_weight[in_bin])
        # Plot these 2d maps -- all four of them.
        cm = plt.cm.get_cmap('RdYlBu')
        vmin = min((np.min(cat_rms), np.min(model_rms)))
        vmax = max((np.max(cat_rms), np.max(model_rms)))
        fig = plt.figure()
        ax = fig.add_subplot(111)
        sc = ax.scatter(mean_snr, mean_res, c=cat_rms, cmap=cm,
                        vmin=vmin, vmax=vmax, edgecolors='none')
        ax.set_xlim(left=10., right=200.)
        plt.colorbar(sc)
        ax.set_title('Using catalog RMS')
        plt.xscale('log')
        plt.xlabel('S/N')
        plt.ylabel('resolution')
        plt.savefig(os.path.join(self.config.pltDir,'intrinsicshape_cat_2d.eps'))

        fig = plt.figure()
        ax = fig.add_subplot(111)
        sc = ax.scatter(mean_snr, mean_res, c=model_rms, cmap=cm,
                        vmin=vmin, vmax=vmax, edgecolors='none')
        ax.set_xlim(left=10., right=200.)
        cbar = plt.colorbar(sc)
        cbar.ax.set_title(r'$e_{\mathrm{rms}}$')
        plt.xscale('log')
        plt.xlabel('S/N')
        plt.ylabel('resolution')
        plt.savefig(os.path.join(self.config.pltDir,'intrinsicshape_2d.eps'))

        fig = plt.figure()
        ax = fig.add_subplot(111)
        sc = ax.scatter(mean_snr, mean_res, c=model_rms, cmap=cm,
                        vmin=np.min(model_rms), vmax=np.max(model_rms), edgecolors='none')
        ax.set_xlim(left=10., right=200.)
        cbar = plt.colorbar(sc)
        cbar.ax.set_title(r'$e_{\mathrm{rms}}$')
        #ax.set_title('Using model RMS')
        plt.xscale('log')
        plt.xlabel('S/N')
        plt.ylabel('resolution')
        plt.savefig(os.path.join(self.config.pltDir,'intrinsicshape_2d_limcolorrange.eps'))


        # save data table for model/actual
        out_data = (mean_snr.flatten(), mean_res.flatten(),
                    model_rms.flatten())
        out_data = np.column_stack(out_data)
        np.savetxt(os.path.join(self.config.pltDir,'intrinsicshape_2d.dat'), out_data)
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
