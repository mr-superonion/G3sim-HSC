#!/usr/bin/env python

import os
import sys
import glob
import fitsio
import optparse
import sim_analysis_utilities as utilities
import numpy as np
from astropy.table import Table, vstack

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


class galDisAnalysis2DConfig(pexConfig.Config):
    anaDir  =   pexConfig.Field(dtype=str, default='anaCat', doc = 'simulation catalogs directory')
    obsDir  =   pexConfig.Field(dtype=str, default='/work/xiangchong.li/work/S16AFPFS/final/noBg/field/', doc = 'observed catalogs directory')
    pltDir  =   pexConfig.Field(dtype=str, default='plot2/galDis/', doc = 'plot directory')
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

class galDisAnalysis2DRunner(TaskRunner):
    @staticmethod
    def getTargetList(parsedCmd, **kwargs):
        return [(ref, kwargs) for ref in range(1)] 

def unpickle(factory, args, kwargs):
    """Unpickle something by calling a factory"""
    return factory(*args, **kwargs)

class galDisAnalysis2DTask(BatchPoolTask):
    ConfigClass     =   galDisAnalysis2DConfig
    RunnerClass     =   galDisAnalysis2DRunner
    _DefaultName    =   "galDisAnalysis2D"
    
    def __init__(self,**kwargs):
        BatchPoolTask.__init__(self, **kwargs)
        self.datanames   =   ['res','snr','mag','absE','weight']
    def __reduce__(self):
        """Pickler"""
        return unpickle, (self.__class__, [], dict(config=self.config, name=self._name,
                parentTask=self._parentTask, log=self.log))
    
    @abortOnError
    def run(self,Id):
        #Prepare the pool
        pool        =   Pool("galDisAnalysis2D")
        pool.cacheClear()
        outcomes    =   pool.map(self.process,utilities.field_keys.keys())
        outcomes    =   [x for x in outcomes if x is not None]
        allObs      =   []
        allSim      =   []
        for outcome in outcomes:
            allObs.append(outcome[0])
            allSim.append(outcome[1])
        
        self.finalPlot(vstack(allSim),vstack(allObs))
        return
    
    def process(self,cache,fieldKey):
        if fieldKey=='aegis':
            return None
        fieldname   =   utilities.field_keys[fieldKey]
        self.log.info('Processing field: %s' %(fieldname))
        anaFieldDir =   os.path.join(self.config.anaDir,'%s' %(fieldname))
        
        self.log.info('Achieving data from observations')
        obsFname=  os.path.join(self.config.obsDir,'%s.fits' %(fieldname))
        obsDat  =   self.getObsGalData(obsFname)
        
        self.log.info('Achieving data from simulations')
        simFname    =   os.path.join(anaFieldDir,'%s.fits' %(fieldname))
        simDat  =   self.getSimGalData(simFname)
        return [obsDat,simDat]
    
    def getObsGalData(self,obsFname):
        data    =   fitsio.read(obsFname)
        mask    =   np.ones(len(data)).astype(bool)
        mask    =   mask&(~np.isnan(data['reGaus_resolution']))
        mask    =   mask&(~np.isnan(data['CModel_snr']))
        mask    =   mask&(~np.isnan(data['CModel_mag']))
        mask    =   mask&(~np.isnan(data['reGaus_sigmaE']))
        mask    =   mask&((data['CModel_mag'])<24.5)
        mask    =   mask&(data['extendness']>0)
        mask    =   mask&(data['reGaus_resolution']>0.3)
        mask    =   mask&(data['CModel_snr']>10.)
        mask    =   mask&(data['reGaus_sigmaE']<2.)
        mask    =   mask&(data['fdfc'])
        data    =   data[mask]
        absE    =   data['reGaus_sigmaE']
        res     =   data['reGaus_resolution']
        mag     =   data['CModel_mag']
        snr     =   data['CModel_snr']
        weight  =   np.ones(len(data))
        dataTab =   [res,snr,mag,absE,weight]
        return Table(data=dataTab,names=self.datanames)
        
        
    def getSimGalData(self,simFname):
        #Get the noise variance and psf information
        data    =   fitsio.read(simFname)
        mask    =   data['weak_lensing_flag']
        data    =   data[mask]
        res     =   data['ext_shapeHSM_HsmShapeRegauss_resolution']
        snr     =   data['modelfit_CModel_flux']/data['modelfit_CModel_fluxSigma']
        absE    =   data['ext_shapeHSM_HsmShapeRegauss_e1']**2.+data['ext_shapeHSM_HsmShapeRegauss_e2']**2.
        absE    =   np.sqrt(absE)
        mag     =   data['Mag']
        weight  =   data['reweight_s16a']
        dataTab =   [res,snr,mag,absE,weight]
        return Table(data=dataTab,names=self.datanames)
    
    
    
    def finalPlot(self,sims,data):
        # Second plot with
        fig2 = plt.figure(figsize=(15, 5))
        ax = fig2.add_subplot(1,3,1)
        counts2,xbins2,ybins2,image = ax.hist2d(sims['snr'], sims['res'], 32, weights=sims['weight'], range=[[10,80],[0.3,1]], cmap='gray', normed=True,vmax=0.08);
        counts1,xbins1,ybins1,image = ax.hist2d(data['snr'], data['res'], 32, weights=data['weight'], range=[[10,80],[0.3,1]], cmap='gray_r', normed=True,vmax=0.08);
        ax.hist2d(data['snr'], data['res'], 128, weights=data['weight'], range=[[10,80],[0.3,1]], cmap='gray_r', normed=True,vmax=0.08,alpha=0.8);

        c = ax.contour(counts1.T,extent=[10,80,0.3,1],cmap='viridis', levels=[0.01, 0.025, 0.05, 0.075],label='data')
        c = ax.contour(counts2.T,extent=[10,80,0.3,1], cmap='viridis', levels=[0.01, 0.025, 0.05, 0.075], linestyles='dashed',label='simulation')
        ax.set_xlabel('CModel S/N')
        ax.set_ylabel('reGaus resolution')
        #ax.set_xscale('log')
        ax.xaxis.set_major_formatter(mtick.FormatStrFormatter('%d'))
        ax.xaxis.set_minor_formatter(mtick.FormatStrFormatter('%d'))

        ax = fig2.add_subplot(1,3,2)

        counts2,xbins2,ybins2,image = ax.hist2d(sims['snr'], sims['absE'], 32, weights=sims['weight'], range=[[10,80],[0.,2]], cmap='gray', normed=True,vmax=0.045);
        counts1,xbins1,ybins1,image = ax.hist2d(data['snr'], data['absE'], 32, weights=data['weight'], range=[[10,80],[0.,2]], cmap='gray_r', normed=True,vmax=0.045);
        ax.hist2d(data['snr'], data['absE'], 128, weights=data['weight'], range=[[10,80],[0.,2]], cmap='gray_r', normed=True,alpha=0.8,norm=LogNorm(vmin=0.0002,vmax=0.2));

        c = ax.contour(counts1.T,extent=[10,80,0.,2],cmap='viridis', levels=[0.0065, 0.012, 0.025, 0.035],label='data')
        c = ax.contour(counts2.T,extent=[10,80,0.,2], cmap='viridis', levels=[0.0065, 0.012, 0.025, 0.035], linestyles='dashed',label='simulation')
        ax.set_xlabel('CModel SNR')
        ax.set_ylabel(r'reGaus $|e|$')
        #ax.set_xscale('log')
        ax.xaxis.set_major_formatter(mtick.FormatStrFormatter('%d'))
        ax.xaxis.set_minor_formatter(mtick.FormatStrFormatter('%d'))

        #  Adding legend
        data_line = mlines.Line2D([], [], color='black', ls='-', label='Data')
        sims_line = mlines.Line2D([], [], color='black', ls='--',label='Simulation')
        ax.legend(handles=[data_line,sims_line])

        ax = fig2.add_subplot(1,3,3)

        counts2,xbins2,ybins2,image = ax.hist2d(sims['snr'], sims['mag'], 32, weights=sims['weight'], range=[[10,80],[22,24.5]], cmap='gray', normed=True,vmax=0.1);
        counts1,xbins1,ybins1,image = ax.hist2d(data['snr'], data['mag'], 32, weights=data['weight'], range=[[10,80],[22,24.5]], cmap='gray_r', normed=True,vmax=0.1);
        ax.hist2d(data['snr'], data['mag'], 128, weights=data['weight'], range=[[10,80],[22,24.5]], cmap='gray_r', normed=True,vmax=0.1,alpha=0.8,norm=LogNorm(vmin=0.005,vmax=0.1));

        c = ax.contour(counts1.T,extent=[10,80,22,24.5],cmap='viridis', levels=[0.01, 0.025, 0.05, 0.075],label='data')
        c = ax.contour(counts2.T,extent=[10,80,22,24.5], cmap='viridis', levels=[0.01, 0.025, 0.05, 0.075], linestyles='dashed',label='simulation')
        ax.set_xlabel('CModel SNR')
        ax.set_ylabel('CModel mag')
        ax.set_ylim(24.5, 22)
        #ax.set_xscale('log')
        ax.xaxis.set_major_formatter(mtick.FormatStrFormatter('%d'))
        ax.xaxis.set_minor_formatter(mtick.FormatStrFormatter('%d'))
        plt.tight_layout()
        fig2.savefig(os.path.join(self.config.pltDir,'galDis2D.eps'))
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
