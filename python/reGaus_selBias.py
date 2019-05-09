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


class reGaus_selBiasConfig(pexConfig.Config):
    anaDir      =   pexConfig.Field(dtype=str, default='anaCat', doc = 'simulation catalogs')
    pltDir      =   pexConfig.Field(dtype=str, default='plot/reGausSelBias/', doc = 'plot directory')
    use_model   =   pexConfig.Field(dtype=bool, default=True, doc = 'whether to use model to estimate weight')
    reweight_s16a=  pexConfig.Field(dtype=bool, default=True, doc = 'reweigting?')
    force_90_pair=  pexConfig.Field(dtype=bool, default=True, doc = 'whether to use orthogonal galaxies to cancel shape noise')
    force_weight=   pexConfig.Field(dtype=bool, default=True, doc = 'test weight bias?')
    
    
    def setDefaults(self):
        pexConfig.Config.setDefaults(self)

    def validate(self):
        pexConfig.Config.validate(self)
        # Checks that the paths provided are valid
        if not os.path.isdir(self.anaDir):
            raise RuntimeError("Simulation catalog directory %s does not exist "
                         "or is not a directory" % (self.anaDir))
        if not os.path.exists(self.pltDir):
            os.mkdir(self.pltDir)

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
        
        if self.config.force_weight:
            self.prefix='shear'
            self.pltApd=''
        else:
            self.prefix='weightBias'
            self.pltApd='wt'
        # Initialize array to store multiplicative bias and associated quantities.
        self.log.info('Beginning selBias analysis:')
        self.log.info('')
        
    def __reduce__(self):
        """Pickler"""
        return unpickle, (self.__class__, [], dict(config=self.config, name=self._name,
                parentTask=self._parentTask, log=self.log))
    
    @abortOnError
    def run(self,Id):
        # For Bias each bins
        self.getBiasBin()
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
        if self.config.force_weight:
            weight_bias=False
        else:
            weight_bias=True
        model_m = utilities.get_m_model(catalog,weight_bias=weight_bias)
        model_c1,model_c2 = utilities.get_c_model(catalog,weight_bias=weight_bias)
        # Get true shear
        g1_true, g2_true= utilities.get_true_shear(catalog[0])
        # Define the PSF shape and rotation angle to apply to get into that frame
        psf_e1, psf_e2  =   utilities.get_psf_ellip(catalog[0])
        weight_s16a     =   catalog['reweight_s16a'][0]
            
        # Get the estimate.
        if np.sum(shape_weights)>0:
            roundtrip_data = (model_m, model_c1, model_c2)
            g1, g2, resp = utilities.estimate_subfield_shear(
                catalog, use_model=self.config.use_model,
                simple_resp=True, report_resp=True,
                reweight_s16a=self.config.reweight_s16a,
                roundtrip_data=roundtrip_data,
                force_weight=self.config.force_weight)
            sigma_g = 0.25/np.sqrt(len(shape_weights))
        else:
            g1= 0.0
            g2= 0.0
            sigma_g= 1.e12
       
        oDat    =   (g1_true,g2_true,g1,g2,sigma_g,psf_e1,psf_e2,weight_s16a)
        return oDat 
    
    def getBiasBin(self):
        self.log.info('Estimating bias for each bin')
        pool    =   Pool("reGaus_selBias")
        pool.cacheClear()
        rows    =   Table.read('noiVar_rPsf_sim.csv')
        outcomes=   pool.map(self.process,rows)
        outcomes=   [x for x in outcomes if x is not None]
        nsim    =   len(outcomes)
        self.g1_true=   np.zeros(nsim);self.g2_true=   np.zeros(nsim)
        self.g1_bins=   np.zeros(nsim);self.g2_bins=   np.zeros(nsim)
        self.sigma_g_bins=   np.zeros(nsim)
        self.psf_e1 =   np.zeros(nsim);self.psf_e2 =   np.zeros(nsim)
        self.weight_s16a=   np.zeros(nsim)
        for iout,outcome in enumerate(outcomes):
            self.g1_true[iout]  =   outcome[0];self.g2_true[iout]  =   outcome[1]
            self.g1_bins[iout]  =   outcome[2];self.g2_bins[iout]  =   outcome[3]
            self.sigma_g_bins[iout]  =   outcome[4]
            self.psf_e1[iout]  =   outcome[5];self.psf_e2[iout]  =   outcome[6]
            self.weight_s16a[iout]  =   outcome[7]
        
        #For each of the sliding bins in SNR/resolution.
        if self.config.reweight_s16a:
            division_val = np.sqrt(self.weight_s16a)
            division_val[division_val<1.e-6] = 1.e-6
            self.sigma_g_bins/= division_val
        tmp_use = self.sigma_g_bins < 1.e8
        self.m1_bins,self.a1_bins, \
        self.sigm1_bins,self.siga1_bins= \
            utilities.get_m_c_a(self.g1_bins[tmp_use],self.g2_bins[tmp_use],
                                self.g1_true[tmp_use], self.g2_true[tmp_use],
                                self.psf_e1[tmp_use],self.psf_e2[tmp_use],
                                self.sigma_g_bins[tmp_use],separate_components=False)
        self.log.info('m = %.3f +/- %.3f' % (self.m1_bins, self.sigm1_bins))
        self.log.info('a = %.3f +/- %.3f' % (self.a1_bins, self.siga1_bins))
        return

    def plotBiasBin(self):
        self.log.info('ploting bias for each bin')
        cm = plt.cm.get_cmap('RdYlBu')

        m_vmin = np.min(self.m1_bins)
        m_vmax = np.max(self.m1_bins)
        fig = plt.figure(figsize=(5,8))
        ax = fig.add_subplot(211)
        sc = ax.scatter(self.snr_bins, self.res_bins, c=self.m1_bins, cmap=cm,
                    vmin=m_vmin, vmax=m_vmax, s=80, edgecolors='none')
        cbar = plt.colorbar(sc)
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

        a_vmin = np.min(self.a1_bins)
        a_vmax = np.max(self.a1_bins)
        fig = plt.figure(figsize=(5,8))
        ax = fig.add_subplot(211)
        sc = ax.scatter(self.snr_bins, self.res_bins, c=self.a1_bins, cmap=cm,vmin=a_vmin, vmax=a_vmax, s=80, edgecolors='none')
        cbar = plt.colorbar(sc)
        plt.xscale('log')
        plt.xlabel('S/N')
        plt.ylabel('resolution')
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
