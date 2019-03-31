#!/usr/bin/env python

import os
import sys
import glob
import optparse
import astropy.io.fits as pyfits
import sim_analysis_utilities as utilities
import numpy as np
from astropy.table import Table, vstack

import matplotlib
matplotlib.use('agg')
matplotlib.rcParams['ps.useafm'] = True
matplotlib.rcParams['pdf.use14corefonts'] = True
matplotlib.rcParams['text.usetex'] = True

import matplotlib.pyplot as plt
import lsst.daf.base as dafBase
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase

from lsst.pipe.base import TaskRunner
from lsst.ctrl.pool.parallel import BatchPoolTask
from lsst.ctrl.pool.pool import Pool, abortOnError


class psfNoiDisPlotConfig(pexConfig.Config):
    anaDir  =   pexConfig.Field(dtype=str, default='anaCat', doc = 'simulation catalogs directory')
    obsDir  =   pexConfig.Field(dtype=str, default='/work/xiangchong.li/work/S16AFPFS/final/noBg/field-matched/', doc = 'observed catalogs directory')
    pltDir  =   pexConfig.Field(dtype=str, default='plot/psfNoi', doc = 'plot directory')
    pltType =   pexConfig.Field(dtype=str, default='eps', doc = 'plot directory')
    version =   pexConfig.Field(dtype=str, default='v1', doc = 'input catalogs')
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

class psfNoiDisPlotRunner(TaskRunner):
    @staticmethod
    def getTargetList(parsedCmd, **kwargs):
        return [(ref, kwargs) for ref in range(1)] 

def unpickle(factory, args, kwargs):
    """Unpickle something by calling a factory"""
    return factory(*args, **kwargs)

class psfNoiDisPlotTask(BatchPoolTask):
    ConfigClass     =   psfNoiDisPlotConfig
    RunnerClass     =   psfNoiDisPlotRunner
    _DefaultName    =   "psfNoiDisPlot"
    
    def __init__(self,**kwargs):
        BatchPoolTask.__init__(self, **kwargs)
        self.nbin1      =   20
        self.fwhm_bin   =   np.linspace(0.3,0.9,self.nbin1)
        self.fwhm_bin_0 =   self.fwhm_bin[:-1]
        self.fwhm_bin_1 =   self.fwhm_bin[1:]
        self.gridPsf    =   (self.fwhm_bin_0+self.fwhm_bin_1)/2.
        self.var_bin    =   np.linspace(0.001,0.01,self.nbin1)
        self.var_bin_0  =   self.var_bin[:-1]
        self.var_bin_1  =   self.var_bin[1:]
        self.gridNoi    =   (self.var_bin_0+self.var_bin_1)/2.
        
    def __reduce__(self):
        """Pickler"""
        return unpickle, (self.__class__, [], dict(config=self.config, name=self._name,
                parentTask=self._parentTask, log=self.log))
    
    @abortOnError
    def run(self,Id):
        tableF      =   []
        obsDis_var  =   np.zeros(self.nbin1-1)
        obsDis_fwhm =   np.zeros(self.nbin1-1)
        obsNum      =   np.zeros(self.nbin1-1)
        
        simDis_var  =   np.zeros(self.nbin1-1)
        simDis_fwhm =   np.zeros(self.nbin1-1)
        simNum      =   np.zeros(self.nbin1-1)
        
        simRDis_var  =   np.zeros(self.nbin1-1)
        simRDis_fwhm =   np.zeros(self.nbin1-1)
        simRNum      =   np.zeros(self.nbin1-1)
        
        
        figN=   plt.figure()
        figP=   plt.figure()
        for ii,fieldKey in enumerate(utilities.field_keys):
            if fieldKey=='aegis':
                continue
            
            fieldname   =   utilities.field_keys[fieldKey]
            anaFieldDir =   os.path.join(self.config.anaDir,'%s' %(fieldname))
            tablename   =   os.path.join(anaFieldDir,'varPSF.csv')
            t           =   Table.read(tablename)
            t['field_name']=    [fieldname]*len(t)
            t.sort('ifield')
            tableF.append(t)
            
            disnameOri  =   os.path.join(anaFieldDir,'varPSF_dis_sim.csv')
            t_dis       =   Table.read(disnameOri)
            simDis_var  +=  t_dis['histNoi']*t_dis['num']
            simDis_fwhm +=  t_dis['histPsf']*t_dis['num']
            simNum      +=  t_dis['num'] 
            
            disnameObs  =    os.path.join(anaFieldDir,'varPSF_dis_obs.csv')
            t_dis_obs   =   Table.read(disnameObs)
            obsDis_var  +=  t_dis_obs['histNoi']*t_dis_obs['num']
            obsDis_fwhm +=  t_dis_obs['histPsf']*t_dis_obs['num']
            obsNum      +=  t_dis_obs['num'] 
            
            disnameRw   =   os.path.join(anaFieldDir,'varPSF_dis_rw.csv')
            t_dis_rw    =   Table.read(disnameRw)
            simRDis_var +=  t_dis_rw['histNoi']*t_dis_rw['num']
            simRDis_fwhm+=  t_dis_rw['histPsf']*t_dis_rw['num']
            simRNum     +=  t_dis_rw['num']
            
            
            axN =   figN.add_subplot(3,3,ii)
            axN.plot(t_dis_obs['gridNoi'],t_dis_obs['histNoi'],'-',label='obs',c='b')
            axN.plot(t_dis_rw['gridNoi'],t_dis_rw['histNoi'],'--',label='reweight',c='r')
            axN.set_xlabel('noise variance')
            axN.set_title(fieldname)
            
            axP =   figP.add_subplot(3,3,ii)
            axP.plot(t_dis_obs['gridPsf'],t_dis_obs['histPsf'],'-',label='obs',c='b')
            axP.plot(t_dis_rw['gridPsf'],t_dis_rw['histPsf'],'--',label='reweight',c='r')
            axP.set_xlabel('PSF FWHM [arcsec]')
            axP.set_title(fieldname)
            
        obsDis_var /=obsNum
        obsDis_fwhm/=obsNum
        
        simDis_var /=simNum
        simDis_fwhm/=simNum
        
        simRDis_var /=simRNum
        simRDis_fwhm/=simRNum
        
        axN =   figN.add_subplot(3,3,7)
        axN.plot(self.gridNoi,obsDis_var,'-',label='obs',c='b')
        axN.plot(self.gridNoi,simRDis_var,'--',label='reweight',c='r')
        axN.set_xlabel('noise variance')
        axN.set_title('All')
        axN.legend(loc='center left', bbox_to_anchor=(1.5, 0.5),fontsize=15)
        figN.tight_layout()
        figN.savefig(os.path.join(self.config.pltDir,'noiVar.eps'))
            
        axP =   figP.add_subplot(3,3,7)
        axP.plot(self.gridPsf,obsDis_fwhm,'-',label='obs',c='b')
        axP.plot(self.gridPsf,simRDis_fwhm,'--',label='reweight',c='r')
        axP.set_xlabel('PSF FWHM [arcsec]')
        axP.set_title('All')
        axP.legend(loc='center left', bbox_to_anchor=(1.5, 0.5),fontsize=15)
        figP.tight_layout()
        figP.savefig(os.path.join(self.config.pltDir,'psfFWHM.eps'))
        
        tableF =    vstack(tableF)
        tableF.sort('ifield')
        tableF.write('noiVar_rPsf_sim.csv',overwrite=True)
        #Prepare the pool
        pool        =   Pool("psfNoiDisPlot")
        pool.cacheClear()
        pool.map(self.writeWeight,tableF)
        return
    
    def writeWeight(self,cache,row):
        ifield  =   row['ifield']
        self.log.info('writing weight to catalog %04d' %(ifield))
        fieldName=  row['field_name']
        anaFname    =   os.path.join(self.config.anaDir,fieldName,'src-%04d.fits'%(ifield))
        catalog     =   Table.read(anaFname)
        catalog['reweight_s16a']=row['weight']
        catalog.write(anaFname,overwrite=True)
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
