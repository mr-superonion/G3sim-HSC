#!/usr/bin/env python
import os
import smatch
import fitsio
import numpy as np
import catSelect
from astropy.table import Table, vstack
import sim_analysis_utilities as utilities

import matplotlib
matplotlib.use('agg')
matplotlib.rcParams['ps.useafm'] = True
matplotlib.rcParams['pdf.use14corefonts'] = True
matplotlib.rcParams['text.usetex'] = True

import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import matplotlib.lines as mlines
from matplotlib.colors import LogNorm
import lsst.daf.base as dafBase
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase

from lsst.pipe.base import TaskRunner
from lsst.ctrl.pool.parallel import BatchPoolTask
from lsst.ctrl.pool.pool import Pool, abortOnError

class galDisAnalysis2DSizeConfig(pexConfig.Config):
    anaDir  =   pexConfig.Field(dtype=str, default='anaCat', doc = 'simulation catalogs directory')
    obsDir  =   pexConfig.Field(dtype=str, default='obsCat', doc = 'observed catalogs directory')
    pltDir  =   pexConfig.Field(dtype=str, default='plot/galDis/', doc = 'plot directory')
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

class galDisAnalysis2DSizeRunner(TaskRunner):
    @staticmethod
    def getTargetList(parsedCmd, **kwargs):
        return [(ref, kwargs) for ref in range(1)] 
def unpickle(factory, args, kwargs):
    """Unpickle something by calling a factory"""
    return factory(*args, **kwargs)
class galDisAnalysis2DSizeTask(BatchPoolTask):
    ConfigClass     =   galDisAnalysis2DSizeConfig
    RunnerClass     =   galDisAnalysis2DSizeRunner
    _DefaultName    =   "galDisAnalysis2DSize"
    def __init__(self,**kwargs):
        BatchPoolTask.__init__(self, **kwargs)
        self.datanames   =   ['res','snr','mag','traceR','logb','weight']
    def __reduce__(self):
        """Pickler"""
        return unpickle, (self.__class__, [], dict(config=self.config, name=self._name,
                parentTask=self._parentTask, log=self.log))

    @abortOnError
    def run(self,Id):
        if not os.path.exists(os.path.join(self.config.anaDir,'allObs4.fits')):
            self.log.info('canot find input files')
            return
        else:
            self.log.info('already have prepared catalog')
            allSim=Table.read(os.path.join(self.config.anaDir,'allSim4.fits'))
            allObs=Table.read(os.path.join(self.config.anaDir,'allObs4.fits'))
        #self.finalPlot_gal(allSim,allObs)
        self.finalPlot_psf2(allSim,allObs)
        return

    def finalPlot_gal(self,sims,data):
        """
        # Plot the snr - (res, traceR, detR) Joint distributions 
        """

        # Do mask
        mask    =   data['inS16L']
        mask    &=  (~data['inS16'])
        data    =   data[mask]

        fig2 = plt.figure(figsize=(15, 5))

        """
        # The First Panel
        """
        ax = fig2.add_subplot(1,3,1)
        counts1,xbins1,ybins1,image = ax.hist2d(data['snr'], data['res'], 32, weights=data['weight'], range=[[10,80],[0.3,1]], cmap='gray_r', normed=True,vmax=0.08);
        c = ax.contour(counts1.T,extent=[10,80,0.3,1],cmap='viridis')
        ax.hist2d(data['snr'], data['res'], 128, weights=data['weight'], range=[[10,80],[0.3,1]], cmap='gray_r', normed=True,vmax=0.08,alpha=0.8);
        ax.set_xlabel('CModel SNR')
        ax.set_ylabel('reGaus resolution')
        #ax.set_xscale('log')
        ax.xaxis.set_major_formatter(mtick.FormatStrFormatter('%d'))
        ax.xaxis.set_minor_formatter(mtick.FormatStrFormatter('%d'))

        """
        # The Second Panel
        """
        ax = fig2.add_subplot(1,3,2)
        counts1,xbins1,ybins1,image = ax.hist2d(data['snr'], data['traceR'], 32, weights=data['weight'], range=[[10,80],[0.,5]], cmap='gray_r', normed=True,vmax=0.045);
        c = ax.contour(counts1.T,extent=[10,80,0.,5.],cmap='viridis')
        ax.hist2d(data['snr'], data['traceR'], 128, weights=data['weight'], range=[[10,80],[0.,5]], cmap='gray_r', normed=True,alpha=0.8,norm=LogNorm(vmin=0.0002,vmax=0.2));
        ax.set_xlabel('CModel SNR')
        ax.set_ylabel('traceR (arcmin)')
        #ax.set_xscale('log')
        ax.xaxis.set_major_formatter(mtick.FormatStrFormatter('%d'))
        ax.xaxis.set_minor_formatter(mtick.FormatStrFormatter('%d'))

        """
        #  Adding legend
        """
        data_line = mlines.Line2D([], [], color='black', ls='-', label='Data')
        sims_line = mlines.Line2D([], [], color='black', ls='--',label='Simulation')
        ax.legend(handles=[data_line,sims_line])

        """
        # The third panel
        """
        ax = fig2.add_subplot(1,3,3)
        counts1,xbins1,ybins1,image = ax.hist2d(data['snr'], data['detR'], 32, weights=data['weight'], range=[[10,80],[0.,5.]], cmap='gray_r', normed=True,vmax=0.1);
        c = ax.contour(counts1.T,extent=[10,80,0.,5.],cmap='viridis')
        ax.hist2d(data['snr'], data['detR'], 128, weights=data['weight'], range=[[10,80],[0.,5.]], cmap='gray_r', normed=True,vmax=0.1,alpha=0.8,norm=LogNorm(vmin=0.0002,vmax=0.2));
        ax.set_xlabel('CModel SNR')
        ax.set_ylabel('detR (arcsec)')
        #ax.set_xscale('log')
        ax.xaxis.set_major_formatter(mtick.FormatStrFormatter('%d'))
        ax.xaxis.set_minor_formatter(mtick.FormatStrFormatter('%d'))

        """
        # Save figure
        """
        plt.tight_layout()
        fig2.savefig(os.path.join(self.config.pltDir,'galDis2D_good_size.png'))
        return

    def finalPlot_psf(self,sims,data):
        """
        # Plot the res - (absEP, traceRP, detRP) Joint distributions 
        """

        # Do mask
        mask    =   data['inS16L']
        mask    &=  (data['inS16'])
        data    =   data[mask]

        fig2 = plt.figure(figsize=(15, 5))

        """
        # The First Panel
        """
        ax = fig2.add_subplot(1,3,1)
        counts1,xbins1,ybins1,image = ax.hist2d(data['res'], data['absEP'], 32, weights=data['weight'], range=[[0.5,1.],[0.,0.25]], cmap='gray_r', normed=True)
        c = ax.contour(counts1.T,extent=[0.5,1,0.,0.25], cmap='viridis')
        ax.hist2d(data['res'], data['absEP'], 128, weights=data['weight'], range=[[0.5,1], [0.,0.25]], cmap='gray_r', normed=True, alpha=0.8)
        ax.set_xlabel('reGaus resolution')
        ax.set_ylabel(r'$|e|$ PSF')

        """
        # The Second Panel
        """
        ax = fig2.add_subplot(1,3,2)
        counts1,xbins1,ybins1,image = ax.hist2d(data['res'], data['traceRP'], 32, weights=data['weight'], range=[[0.5,1], [0.2, 0.6]], cmap='gray_r', normed=True)
        c = ax.contour(counts1.T,extent=[0.5,1,0.2,0.6], cmap='viridis')
        ax.hist2d(data['res'], data['traceRP'], 128, weights=data['weight'], range=[[0.5,1], [0.2,0.6]], cmap='gray_r', normed=True, alpha=0.8)
        ax.set_xlabel('reGaus resolution')
        ax.set_ylabel('traceR PSF (arcmin)')

        """
        #  Adding legend
        """
        data_line = mlines.Line2D([], [], color='black', ls='-', label='Data')
        sims_line = mlines.Line2D([], [], color='black', ls='--',label='Simulation')
        ax.legend(handles=[data_line,sims_line])

        """
        # The third panel
        """
        ax = fig2.add_subplot(1,3,3)
        counts1,xbins1,ybins1,image = ax.hist2d(data['res'], data['detRP'], 32, weights=data['weight'], range=[[0.5,1.], [0.1, 0.5]], cmap='gray_r', normed=True)
        c = ax.contour(counts1.T,extent=[0.5,1.,0.1,0.5],cmap='viridis')
        ax.hist2d(data['res'], data['detRP'], 128, weights=data['weight'], range=[[0.5,1], [0.1,0.5]], cmap='gray_r', normed=True, alpha=0.8)
        ax.set_xlabel('reGaus resolution')
        ax.set_ylabel('detR PSF (arcmin)')

        """
        # Save figure
        """
        plt.tight_layout()
        fig2.savefig(os.path.join(self.config.pltDir,'galDis2D_good_size_psf_largeRes.png'))
        return

    def finalPlot_psf2(self,sims,data):
        """
        # Plot the res - (absEP, traceRP, detRP) Joint distributions 
        """

        # Do mask
        mask    =   data['inS16L']
        mask    &=  (data['inS16'])
        data    =   data[mask]

        fig2 = plt.figure(figsize=(15, 5))

        """
        # The First Panel
        """
        ax = fig2.add_subplot(1,3,1)
        counts1,xbins1,ybins1,image = ax.hist2d(data['res'], data['absEP'], 32, weights=data['weight'], range=[[0.5,1.],[0.,0.25]], cmap='gray_r', normed=True)
        c = ax.contour(counts1.T,extent=[0.5,1,0.,0.25], cmap='viridis')
        ax.hist2d(data['res'], data['absEP'], 128, weights=data['weight'], range=[[0.5,1], [0.,0.25]], cmap='gray_r', normed=True, alpha=0.8)
        ax.set_xlabel('reGaus resolution')
        ax.set_ylabel(r'$|e|$ PSF')

        """
        # The Second Panel
        """
        ax = fig2.add_subplot(1,3,2)
        counts1,xbins1,ybins1,image = ax.hist2d(data['absEP'], data['traceRP'], 32, weights=data['weight'], range=[[0.,0.25], [0.2, 0.6]], cmap='gray_r', normed=True)
        c = ax.contour(counts1.T,extent=[0.,0.25,0.2,0.6], cmap='viridis')
        ax.hist2d(data['res'], data['traceRP'], 128, weights=data['weight'], range=[[0.,0.25], [0.2,0.6]], cmap='gray_r', normed=True, alpha=0.8)
        ax.set_xlabel(r'$|e|$ PSF')
        ax.set_ylabel('traceR PSF (arcmin)')

        """
        #  Adding legend
        """
        data_line = mlines.Line2D([], [], color='black', ls='-', label='Data')
        sims_line = mlines.Line2D([], [], color='black', ls='--',label='Simulation')
        ax.legend(handles=[data_line,sims_line])

        """
        # The third panel
        """
        ax = fig2.add_subplot(1,3,3)
        counts1,xbins1,ybins1,image = ax.hist2d(data['absEP'], data['detRP'], 32, weights=data['weight'], range=[[0.,0.25], [0.1, 0.5]], cmap='gray_r', normed=True)
        c = ax.contour(counts1.T,extent=[0.,0.25,0.1,0.5],cmap='viridis')
        ax.hist2d(data['res'], data['detRP'], 128, weights=data['weight'], range=[[0.,0.25], [0.1,0.5]], cmap='gray_r', normed=True, alpha=0.8)
        ax.set_xlabel(r'$|e|$ PSF')
        ax.set_ylabel('detR PSF (arcmin)')

        """
        # Save figure
        """
        plt.tight_layout()
        fig2.savefig(os.path.join(self.config.pltDir,'galDis2D_good_elli_psf_largeRes.png'))
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
