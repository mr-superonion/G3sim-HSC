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


class reGaus_shearCalibConfig(pexConfig.Config):
    anaDir      =   pexConfig.Field(dtype=str, default='anaCat', doc = 'simulation catalogs')
    pltDir      =   pexConfig.Field(dtype=str, default='plot/reGausCalib/', doc = 'plot directory')
    nbin_1d     =   pexConfig.Field(dtype=int, default=8, doc = 'number of bins for 1-D parameter space')
    ngal_bin    =   pexConfig.Field(dtype=int, default=150000, doc = 'number of bins for 1-D parameter space')
    use_model   =   pexConfig.Field(dtype=bool, default=True, doc = 'whether to use model to estimate weight')
    reweight_s16a=  pexConfig.Field(dtype=bool, default=True, doc = 'reweigting?')
    force_90_pair=  pexConfig.Field(dtype=bool, default=True, doc = 'whether to use orthogonal galaxies to cancel shape noise')
    force_weight=   pexConfig.Field(dtype=bool, default=True, doc = 'test weight bias?')
    round_trip  =   pexConfig.Field(dtype=bool, default=False, doc = 'do a round trip?')
    test_weight_bias  =   pexConfig.Field(dtype=bool, default=False, doc = 'test weight bias?')
    
    
    def setDefaults(self):
        pexConfig.Config.setDefaults(self)

    def validate(self):
        pexConfig.Config.validate(self)
        if self.test_weight_bias and (not self.round_trip):
            raise RuntimeError("Cannot test weight bias when not doing round trip")
        if self.test_weight_bias and self.force_weight:
            raise RuntimeError("Cannot test weight bias when forcing weight")
        # Checks that the paths provided are valid
        if not os.path.isdir(self.anaDir):
            raise RuntimeError("Simulation catalog directory %s does not exist "
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
        
        if self.config.force_weight:
            self.prefix='shear'
            self.pltApd=''
        else:
            self.prefix='weightBias'
            self.pltApd='wt'
        if self.config.round_trip:
            self.prefix='test_'+self.prefix
        
        self.add_m_val  =   0.
        self.add_a_val  =   0.

        # Initialize array to store multiplicative bias and associated quantities.
        self.dshape     =   (self.config.nbin_1d,self.config.nbin_1d)
        self.snr_bins   =   np.zeros(self.dshape)
        self.res_bins   =   np.zeros(self.dshape)

        self.fake_x     =   None
        self.m1_bins    =   np.zeros(self.dshape)
        self.a1_bins    =   np.zeros(self.dshape)
        self.sigm1_bins =   np.zeros(self.dshape)
        self.siga1_bins =   np.zeros(self.dshape)
        
        
        self.log.info('Beginning shearCalib analysis for the following scenario:')
        self.log.info('  %d bins in both SNR and resolution'%self.config.nbin_1d)
        self.log.info('')
        
    def __reduce__(self):
        """Pickler"""
        return unpickle, (self.__class__, [], dict(config=self.config, name=self._name,
                parentTask=self._parentTask, log=self.log))
    
    @abortOnError
    def run(self,Id):
        #Prepare the pool
        slide_cat   =   utilities.prepare_slide_catalog(self.config.anaDir) 
        self.separate_slide(slide_cat)
        # For Bias each bins
        self.getBiasBin()
        self.plotBiasBin()
        # For Bias model
        biasModel   =   self.fitModel()
        self.plotModel(biasModel)
        return
        
    def separate_slide(self,slide_cat):
        self.log.info('Detemining the sliding window')
        mask    =   utilities.galaxy_selector(catalog=slide_cat,
                    max_snr = 1000.,max_sigma_e=0.35,max_e =2.,
                    use_wl_flag=True,force_90_pair=self.config.force_90_pair)
        slide_cat=  slide_cat[mask]
        # Set up sliding window definition and make a mask for each SNR, resolution bin.
        snr = utilities.get_snr(slide_cat)
        res = utilities.get_res(slide_cat)
        # Sort out binning. 
        if self.config.reweight_s16a:
            weights =   slide_cat['reweight_s16a']
        else:
            weights =   None
        self.snr_lower, self.snr_upper, self.res_lower, self.res_upper = \
            utilities.sliding_window_def(snr, res,weights,
            n_gal=self.config.ngal_bin,n_bin_1d=self.config.nbin_1d)
        
        # We need to get the weights here, to use when averaging over galaxies to get mean effective
        # SNR and resolution.
        # Note, it might be worth making a utility for this combined shape + S16A weight, since it's
        # done here and in estimate_subfield_shear.
        if self.config.use_model: 
            shape_weights = utilities.get_weight_model(slide_cat)
        else:
            shape_weights = slide_cat['ishape_hsm_regauss_derived_shape_weight']
        if self.config.reweight_s16a:
            shape_weights = shape_weights*slide_cat['reweight_s16a']
        if self.config.force_weight:
            shape_weights[1::2] = shape_weights[0::2]
            assert np.all(slide_cat['reweight_s16a'][0::2]==slide_cat['reweight_s16a'][1::2]), \
                'Galaxies in pairs do not have the same reweight_s16a factors'
        
        # Record the mean (snr,res) for each bin
        for j in range(self.config.nbin_1d):
            for k in range(self.config.nbin_1d):
                # This line below doesn't preserve shape noise cancellation, since each galaxy in
                # the pair has the cut imposed separately.
                tmp_mask = (snr >= self.snr_lower[j]) & (snr < self.snr_upper[j]) & \
                (res >= self.res_lower[j,k]) & (res < self.res_upper[j,k])
                # This line will re-impose the shape noise cancellation (by only using the first in
                # the pair for the cut).
                if self.config.force_90_pair:
                    tmp_mask[1::2] = tmp_mask[0::2]
                mean_snr = np.sum(shape_weights[tmp_mask]*snr[tmp_mask])/ \
                    np.sum(shape_weights[tmp_mask])
                mean_res = np.sum(shape_weights[tmp_mask]*res[tmp_mask])/ \
                    np.sum(shape_weights[tmp_mask])
                self.snr_bins[j,k] = mean_snr
                self.res_bins[j,k] = mean_res
        self.fake_x = np.vstack((self.snr_bins.flatten(), self.res_bins.flatten()))
        return
        
        
    
    def process(self,cache,row):
        isim    =   row['ifield']
        fN      =   row['field_name']
        self.log.info('Beginning for field: %s, sim: %s' %(fN,isim))
        catFname=   os.path.join(self.config.anaDir,'%s'%fN,'simCat','%d.fits' %isim)
        catalog =   Table.read(catFname)
        #Estimate shear for each isim
        mask    =   utilities.galaxy_selector(catalog=catalog,
                    max_snr = 1000.,max_sigma_e=0.35,max_e =2.,
                     use_wl_flag=True,force_90_pair=self.config.force_90_pair)
        catalog =  catalog[mask]
        if len(catalog)<2:
            return None
        # Set up sliding window definition and make a mask for each SNR, resolution bin.
        snr = utilities.get_snr(catalog)
        res = utilities.get_res(catalog)
        
        if self.config.use_model: 
            shape_weights = utilities.get_weight_model(catalog)
        else:
            shape_weights =catalog['ishape_hsm_regauss_derived_shape_weight']
        if self.config.reweight_s16a:
            shape_weights = shape_weights*catalog['reweight_s16a']
        if self.config.force_weight:
            shape_weights[1::2] = shape_weights[0::2]
            assert np.all(catalog['reweight_s16a'][0::2]==catalog['reweight_s16a'][1::2]), \
                'Galaxies in pairs do not have the same reweight_s16a factors'
        if self.config.round_trip:
            model_m = utilities.get_m_model(catalog,weight_bias=self.config.test_weight_bias)
            model_c1,model_c2 = utilities.get_c_model(catalog,weight_bias=self.config.test_weight_bias)
        # The shear estimation routine is built to work on simulation subfields
        # We need to iterate on all independent subfields
        g1_bins = np.zeros(self.dshape)
        g2_bins = np.zeros(self.dshape)
        sigma_g_bins= np.zeros(self.dshape)
        # Get true shear
        g1_true, g2_true= utilities.get_true_shear(catalog[0])
        # Define the PSF shape and rotation angle to apply to get into that frame
        psf_e1, psf_e2  =   utilities.get_psf_ellip(catalog[0])
        weight_s16a     =   catalog['reweight_s16a'][0]
            
        # Get the estimate in bins in SNR/res as well.
        for j in range(self.config.nbin_1d):
            for k in range(self.config.nbin_1d):
                mm= (snr >= cache.snr_lower[j])&(snr < cache.snr_upper[j])&\
                (res>=cache.res_lower[j,k])&(res < cache.res_upper[j,k])
                if self.config.force_90_pair:
                    mm[1::2] = mm[0::2]
                # Some weird bins might have no objects in particular subfields.
                # Don't use estimate_subfield_shear on those!
                if mm.astype(int).sum()>2 and np.sum(shape_weights[mm])>0:
                    if self.config.round_trip:
                        roundtrip_data = (model_m[mm], model_c1[mm], model_c2[mm])
                    else:
                        roundtrip_data = None
                    tmp_g1, tmp_g2, tmp_resp = utilities.estimate_subfield_shear(
                        catalog[mm], use_model=self.config.use_model,
                        simple_resp=True, report_resp=True,
                        reweight_s16a=self.config.reweight_s16a,
                        roundtrip_data=roundtrip_data,
                        force_weight=self.config.force_weight)
                    g1_bins[j,k] = tmp_g1
                    g2_bins[j,k] = tmp_g2
                    sigma_g_bins[j,k] = 0.25/np.sqrt(mm.astype(int).sum())
                else:
                    g1_bins[j,k] = 0.0
                    g2_bins[j,k] = 0.0
                    sigma_g_bins[j,k] = 1.e12
        
        oDat    =   (g1_true,g2_true,g1_bins,g2_bins,sigma_g_bins
                    ,psf_e1,psf_e2,weight_s16a)
        return oDat 
    
    def getBiasBin(self):
        biasBinFname    =   self.config.pltDir+'%s_m_a_2d.csv'%(self.prefix)
        if os.path.exists(biasBinFname):
            biasBinTab  =   Table.read(biasBinFname)
            self.snr_bins=  np.array(biasBinTab['snr']).reshape(self.dshape)
            self.res_bins=  np.array(biasBinTab['res']).reshape(self.dshape)
            self.m1_bins=  np.array(biasBinTab['m']).reshape(self.dshape)
            self.sigm1_bins=  np.array(biasBinTab['mSigma']).reshape(self.dshape)
            self.a1_bins=  np.array(biasBinTab['a']).reshape(self.dshape)
            self.siga1_bins=  np.array(biasBinTab['aSigma']).reshape(self.dshape)
            return
        
        self.log.info('Estimating bias for each bin')
        pool    =   Pool("reGaus_shearCalib")
        pool.cacheClear()
        pool.storeSet(snr_lower=self.snr_lower,snr_upper=self.snr_upper,
                    res_lower=self.res_lower,res_upper=self.res_upper)
        rows    =   Table.read('noiVar_rPsf_sim.csv')
        outcomes=   pool.map(self.process,rows)
        outcomes=   [x for x in outcomes if x is not None]
        nsim    =   len(outcomes)
        self.g1_true=   np.zeros(nsim)
        self.g2_true=   np.zeros(nsim)
        self.g1_bins=   np.zeros((self.config.nbin_1d,self.config.nbin_1d,nsim))
        self.g2_bins=   np.zeros((self.config.nbin_1d,self.config.nbin_1d,nsim))
        self.sigma_g_bins=   np.zeros((self.config.nbin_1d,self.config.nbin_1d,nsim))
        self.psf_e1 =   np.zeros(nsim)
        self.psf_e2 =   np.zeros(nsim)
        self.weight_s16a=   np.zeros(nsim)
        for iout,outcome in enumerate(outcomes):
            self.g1_true[iout]  =   outcome[0]
            self.g2_true[iout]  =   outcome[1]
            self.g1_bins[:,:,iout]  =   outcome[2]
            self.g2_bins[:,:,iout]  =   outcome[3]
            self.sigma_g_bins[:,:,iout]  =   outcome[4]
            self.psf_e1[iout]  =   outcome[5]
            self.psf_e2[iout]  =   outcome[6]
            self.weight_s16a[iout]  =   outcome[7]
        
        #For each of the sliding bins in SNR/resolution.
        for j in range(self.config.nbin_1d):
            for k in range(self.config.nbin_1d):
                sigma_vals = self.sigma_g_bins[j,k,:]
                if self.config.reweight_s16a:
                    division_val = np.sqrt(self.weight_s16a)
                    division_val[division_val<1.e-6] = 1.e-6
                    sigma_vals /= division_val
                tmp_use = sigma_vals < 1.e10
                self.m1_bins[j,k],self.a1_bins[j,k], \
                self.sigm1_bins[j,k],self.siga1_bins[j,k]= \
                    utilities.get_m_c_a(self.g1_bins[j,k,tmp_use],
                                        self.g2_bins[j,k,tmp_use],
                                        self.g1_true[tmp_use], 
                                        self.g2_true[tmp_use],
                                        self.psf_e1[tmp_use],
                                        self.psf_e2[tmp_use],
                                        sigma_vals[tmp_use],
                                        separate_components=False)
                self.log.info('Bin %d %d bias analysis:'%(j,k))
                self.log.info('m = %.3f +/- %.3f' % (self.m1_bins[j,k], self.sigm1_bins[j,k]))
                self.log.info('a = %.3f +/- %.3f' % (self.a1_bins[j,k], self.siga1_bins[j,k]))
        if not self.config.force_weight:
            shearName   =   self.config.pltDir+'shear_m_a_2d.csv'
            shearTab    =   Table.read(shearName)
            self.m1_bins-= np.array(shearTab['m']).reshape(self.dshape)
            self.a1_bins-= np.array(shearTab['a']).reshape(self.dshape)

        #Write data 
        out_data    =   (self.snr_bins.flatten(), self.res_bins.flatten(),
                self.m1_bins.flatten(), self.sigm1_bins.flatten(),
                self.a1_bins.flatten(), self.siga1_bins.flatten())

        mcBinTab    =   Table(data=out_data,names=('snr','res','m','mSigma','a','aSigma'))
        mcBinTab.write(biasBinFname)
        return

    def plotBiasBin(self):
        self.log.info('ploting bias for each bin')
        cm = plt.cm.get_cmap('RdYlBu')

        m_vmin = np.min(self.m1_bins+self.add_m_val)
        m_vmax = np.max(self.m1_bins+self.add_m_val)
        fig = plt.figure(figsize=(5,8))
        ax = fig.add_subplot(211)
        sc = ax.scatter(self.snr_bins, self.res_bins, c=self.m1_bins+self.add_m_val, cmap=cm,
                    vmin=m_vmin, vmax=m_vmax, s=80, edgecolors='none')
        cbar = plt.colorbar(sc)
        if self.add_m_val > 0:
            cbar.ax.set_title(r'$m+$constant')
        else:
            cbar.ax.set_title(r'$m_{%s}$'%(self.pltApd))
        plt.xscale('log')
        plt.xlabel('S/N')
        plt.ylabel('resolution')
        ax = fig.add_subplot(212)
        sc = ax.scatter(self.snr_bins, self.res_bins, c=self.m1_bins/self.sigm1_bins, cmap=cm, s=80, edgecolors='none')
        cbar = plt.colorbar(sc)
        cbar.ax.set_title(r'$m_{%s}/\sigma_m$'%(self.pltApd) )
        plt.xscale('log')
        plt.xlabel('S/N')
        plt.ylabel('resolution')
        plt.tight_layout()
        plt.savefig(self.config.pltDir+'%s_m_2d.eps'%(self.prefix))

        a_vmin = np.min(self.a1_bins+self.add_a_val)
        a_vmax = np.max(self.a1_bins+self.add_a_val)
        fig = plt.figure(figsize=(5,8))
        ax = fig.add_subplot(211)
        sc = ax.scatter(self.snr_bins, self.res_bins, c=self.a1_bins+self.add_a_val, cmap=cm,vmin=a_vmin, vmax=a_vmax, s=80, edgecolors='none')
        cbar = plt.colorbar(sc)
        plt.xscale('log')
        plt.xlabel('S/N')
        plt.ylabel('resolution')
        if self.add_a_val > 0:
            cbar.ax.set_title(r'$m+$constant')
        else:
            cbar.ax.set_title(r'$a_{%s}$' %(self.pltApd))
        # Also plot errors in this case.
        ax = fig.add_subplot(212)
        sc = ax.scatter(self.snr_bins, self.res_bins, c=self.a1_bins/self.siga1_bins, cmap=cm, s=80, edgecolors='none')
        cbar = plt.colorbar(sc)
        cbar.ax.set_title(r'$a_{%s}/\sigma_a$' %(self.pltApd))
        plt.xscale('log')
        plt.xlabel('S/N')
        plt.ylabel('resolution')
        plt.tight_layout()
        plt.savefig(self.config.pltDir+'%s_a_2d.eps'%(self.prefix))
        return

    def fitModel(self):
        maFname  =   self.config.pltDir+'%s_m_a_model_par.npy'%(self.prefix)
        if os.path.exists(maFname):
            return np.load(maFname).item()
        outMA   =   {}
        # Also do model fitting.
        m_opt, m_cov = opt.curve_fit(
            utilities.m_func, self.fake_x, self.m1_bins.flatten(), sigma=self.sigm1_bins.flatten(),
            bounds=([-np.inf,0,-np.inf,-np.inf],[np.inf,np.inf,0,np.inf]))

        outMA.update({'m_opt':m_opt,'m_cov':m_cov})
        self.log.info('m model:')
        self.log.info('m=b*(res/0.5)^c*(snr/20)^d+e')
        for ind in range(len(m_opt)):
            self.log.info('%f +/- %f'%(m_opt[ind],m_cov[ind,ind])) 
        
        #For additive bias
        a_opt, a_cov = opt.curve_fit(
            utilities.a_func, self.fake_x, self.a1_bins.flatten(), sigma=self.siga1_bins.flatten(),
            bounds=([-np.inf,-np.inf,-np.inf],[np.inf,np.inf,0.]))
        
        outMA.update({'a_opt':a_opt,'a_cov':a_cov})
        self.log.info('a model:')
        self.log.info('a=b*(res-c)*(snr/20)^d')
        for ind in range(len(a_opt)):
            self.log.info('%f +/- %f'%(a_opt[ind],a_cov[ind,ind]))
        np.save(maFname,outMA)
        return outMA

    def plotModel(self,fitModel):
        cm = plt.cm.get_cmap('RdYlBu')
        m_vmin = np.min(self.m1_bins+self.add_m_val)
        m_vmax = np.max(self.m1_bins+self.add_m_val)
        a_vmin = np.min(self.a1_bins+self.add_a_val)
        a_vmax = np.max(self.a1_bins+self.add_a_val)
        m_opt   =   fitModel['m_opt']
        a_opt   =   fitModel['a_opt']
        # plot the model, residuals, and their significance
        m_model =   utilities.m_func(self.fake_x, *m_opt).reshape(self.dshape)
        deltaM  =   self.m1_bins-m_model
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
        sc = ax.scatter(self.snr_bins, self.res_bins, c=deltaM,
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
    
        # Plot model, residuals, significance
        a_model = utilities.a_func(self.fake_x, *a_opt).reshape(self.config.nbin_1d, self.config.nbin_1d)
        deltaA  =   self.a1_bins-a_model
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
        sc = ax.scatter(self.snr_bins,self.res_bins,c=deltaA, 
                        cmap=cm, s=80, edgecolors='none')
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
        out_data=(self.snr_bins.flatten(), self.res_bins.flatten(),
                deltaM.flatten(), self.sigm1_bins.flatten(),
                deltaA.flatten(), self.siga1_bins.flatten())
        mcDeltaBinTab    =   Table(data=out_data,names=('snr','res','dm','mSigma','da','aSigma'))
        mcDeltaBinTab.write(self.config.pltDir+'%s_dm_da_2d.csv'%(self.prefix))
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
