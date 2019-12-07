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

class galDisAnalysis2DConfig(pexConfig.Config):
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
        self.datanames   =   ['res','snr','mag','absE','logb','weight']
    def __reduce__(self):
        """Pickler"""
        return unpickle, (self.__class__, [], dict(config=self.config, name=self._name,
                parentTask=self._parentTask, log=self.log))

    @abortOnError
    def run(self,Id):
        if not os.path.exists(os.path.join(self.config.anaDir,'allObs4.fits')):
            #Prepare the pool
            pool        =   Pool("galDisAnalysis2D")
            pool.cacheClear()
            outcomes    =   pool.map(self.process,utilities.field_keys.keys())
            outcomes    =   [x for x in outcomes if x is not None]
            allObsList  =   []
            allSimList  =   []
            for outcome in outcomes:
                allObsList.append(outcome[0])
                allSimList.append(outcome[1])
            allSim      =   vstack(allSimList)
            allObs      =   vstack(allObsList)
            del outcomes
            del allObsList
            del allSimList
            allSim.write(os.path.join(self.config.anaDir,'allSim4.fits'))
            allObs.write(os.path.join(self.config.anaDir,'allObs4.fits'))
        else:
            self.log.info('already have prepared catalog')
            allSim=Table.read(os.path.join(self.config.anaDir,'allSim4.fits'))
            allObs=Table.read(os.path.join(self.config.anaDir,'allObs4.fits'))
        """
        self.logbNumPlot(allObs)
        self.logbPlot2D(allObs)
        self.compareS16(allObs)
        """
        self.finalPlot(allSim,allObs)
        return

    def process(self,cache,fieldKey):
        fieldname   =   utilities.field_keys[fieldKey]
        self.log.info('Processing field: %s' %(fieldname))
        anaFieldDir =   os.path.join(self.config.anaDir,'%s' %(fieldname))
 
        self.log.info('Achieving galaxy data from observations')
        obsFname    =  os.path.join(self.config.obsDir,'%s.fits' %(utilities.field_names[fieldname]))
        obsDat      =   self.getObsGalData(obsFname,fieldname)
 
        self.log.info('Achieving galaxy data from simulations')
        simFname    =   os.path.join(anaFieldDir,'%s.fits' %(fieldname))
        simDat      =   self.getSimGalData(simFname)
        return [obsDat,simDat]

    def getObsGalData(self,obsFname,fieldname):
        """
        #load the necessary data from catalog (after FDFC, WL cut)
        """
        data    =   fitsio.read(obsFname)
        """
        # Match to S16A
        """
        obsFname16= os.path.join(os.environ['homeWrk'],\
                'S16AStandard/S16AStandardCalibrated/field',\
                '%s_calibrated.fits' %fieldname)
        
        data16  =   fitsio.read(obsFname16)
        ra      =   data['i_ra']
        dec     =   data['i_dec']
        """
        # Match to S16A
        """
        matchL  =   smatch.match(data16['ira'],data16['idec'],20./3600.,\
                ra,dec,maxmatch=1)
        inS16L  =   np.zeros(len(data)).astype(bool)
        inS16L[matchL['i2']]  =   True
        mDist   =   np.arccos(matchL['cosdist'])*180./np.pi

        matchM  =  matchL[mDist<10./3600.]
        inS16M  =   np.zeros(len(data)).astype(bool)
        inS16M[matchM['i2']]  =   True

        matchS  =   matchL[mDist<5./3600.]
        inS16S  =   np.zeros(len(data)).astype(bool)
        inS16S[matchS['i2']]  =   True

        match   =   matchL[mDist<0.5/3600.]
        inS16   =   np.zeros(len(data)).astype(bool)
        inS16[match['i2']]  =   True

        """
        # Main Column for S16A
        """
        res16   =   np.zeros(len(data))
        res16[match['i2']]  =   utilities.get_res(data16[match['i1']])
        snr16   =   np.zeros(len(data))
        snr16[match['i2']]  =   utilities.get_snr(data16[match['i1']])
        mag16   =   np.zeros(len(data))
        mag16[match['i2']]  =   utilities.get_imag(data16[match['i1']])
        logb16   =   np.zeros(len(data))
        logb16[match['i2']]  =   utilities.get_logb(data16[match['i1']])
        absE16  =   np.zeros(len(data))
        absE16[match['i2']] =   utilities.get_abs_ellip(data16[match['i1']])
        fwhm16   =   np.zeros(len(data))
        fwhm16[match['i2']]  =   utilities.get_psf_size(data16[match['i1']])

        """
        # FDFC Column for S16A
        """
        inS16FDFC   =   catSelect.getFDFCS16(data,'i_ra','i_dec')

        """
        # Main Column for S18A
        """
        res     =   utilities.get_res(data)
        snr     =   utilities.get_snr(data)
        mag     =   utilities.get_imag(data)
        absE    =   utilities.get_abs_ellip(data)
        logb    =   utilities.get_logb(data)
        fwhm    =   utilities.get_psf_size(data,'fwhm')
        traceR  =   utilities.get_gal_size(data,type='trace')
        detR    =   utilities.get_gal_size(data,type='det')
        traceRP =   utilities.get_psf_size(data,type='trace')
        detRP   =   utilities.get_psf_size(data,type='det')
        absEP   =   utilities.get_abs_ellip_psf(data)
        fieldnameA=[fieldname]*len(data)

        """
        # Flags for S18A
        """
        badPsf  =   data['i_pixelflags_inexact_psf']
        badPsfC =   data['i_pixelflags_inexact_psfcenter']
        bO_rep  =   data['i_pixelflags_bright_object_ireprocess']
        bO_repC =   data['i_pixelflags_bright_objectcenter_ireprocess']
        bO_S18  =   data['mask_briObj_S18']
        objId   =   data['object_id']
        weight  =   np.ones(len(data))
        datanmA =   ['fwhm', 'badPsf', 'badPsfC', 'bO_rep', 'bO_repC', 'bO_S18_r', 'inS16', 'inS16L', 'inS16FDFC', \
                'res16', 'snr16', 'mag16', 'absE16', 'logb16', 'fwhm16', 'fieldname', 'ra', 'dec', 'object_id', \
                'traceR', 'detR', 'traceRP', 'detRP', 'absEP', 'inS16M', 'inS16S']
        dataTab =   [res, snr, mag, absE, logb, weight, fwhm, badPsf, badPsfC, bO_rep, \
                bO_repC, bO_S18, inS16, inS16L, inS16FDFC, \
                res16, snr16, mag16, absE16, logb16, fwhm16, fieldnameA, ra, dec, objId, \
                traceR, detR, traceRP, detRP, absEP, inS16M, inS16S]
        return Table(data=dataTab,names=self.datanames+datanmA)

    def getSimGalData(self,simFname):
        data    =   fitsio.read(simFname)
        data    =   utilities.update_wl_flags_sim(data)
        mask    =   data['weak_lensing_flag']
        """
        mask    =   mask&(utilities.get_snr(data) < 80.)
        """
        data    =   data[mask]
        absE    =   utilities.get_abs_ellip(data)
        res     =   utilities.get_res(data)
        mag     =   utilities.get_imag(data)
        snr     =   utilities.get_snr(data)
        logb    =   utilities.get_logb(data)
        weight  =   data['reweight_s16a']
        dataTab =   [res,snr,mag,absE,logb,weight]
        return Table(data=dataTab,names=self.datanames)

    def logbNumPlot(self,data):
        """
        # Plot the galaxy  number loss vs. logb
        """
        fig =   plt.figure(figsize=(6, 6))
        ax  =   fig.add_subplot(1,1,1)

        self.log.info('beginning logbNumPlot')
        data    =   data[data['bO_S18_r'].astype(bool)]
        """
        # Prepare 
        """
        numAll  =   len(data)
        thres   =   np.arange(-1.2,-0.375,0.05)
        frac    =   np.array([np.sum(data['logb']<t) for t in thres])/numAll
        ax.plot(thres,frac)
        ax.set_xlabel('log10(b)')
        ax.set_ylabel('fraction')
        ax.grid()
        """
        # Plot
        """
        plt.tight_layout()
        fig.savefig(os.path.join(self.config.pltDir,'frac-logb_newBriObjMsk.png'))
        return

    def logbPlot2D(self,data):
        """
        # Plot the logb vs. resolution
        """
        fig =   plt.figure(figsize=(6, 6))
        ax  =   fig.add_subplot(1,1,1)

        data    =   data[data['bO_S18_r'].astype(bool)]

        """
        # Problematic data (extended faint galaxies)
        """
        mask=   (data['res']>0.9)&(data['snr']<15.)
        data2=  data[mask]
        counts,xbins,ybins,image = ax.hist2d(data2['logb'], data2['res'],8, weights=data2['weight'], range=[[-2,-.375],[0.9,1]], cmap='gray_r', normed=True)
        c   =   ax.contour(counts.T,extent=[-2,-0.375,0.9,1], cmap='gray_r',linestyles='dashed')
        print(c.levels)

        """
        # All of the data
        """
        hd  =   ax.hist2d(data['logb'], data['res'], 64, weights=data['weight'], range=[[-2,-0.375],[0.3,1]], cmap='jet', normed=True,alpha=0.8)
        ax.set_xlabel('log10(b)')
        ax.set_ylabel('reGaus resolution')
        cbar=   fig.colorbar(hd[3],ax=ax)

        """
        # Plot
        """
        plt.tight_layout()
        fig.savefig(os.path.join(self.config.pltDir,'logb-res_newBriObjMsk.png'))
        return

    def compareS16(self,data):
        mask    =   data['bO_S18_r'].astype(bool)
        mask    &=  data['inS16']
        data    =   data[mask]
        fig2    =   plt.figure(figsize=(15, 10))

        """
        # SNR scatter
        """
        ax = fig2.add_subplot(2,3,1)
        ax.scatter(data['snr'],data['snr16'],s=0.1)
        ax.set_xlabel('SNR S18')
        ax.set_ylabel('SNR S16')

        """
        # Mag scatter
        """
        ax = fig2.add_subplot(2,3,2)
        ax.scatter(data['mag'],data['mag16'],s=0.1)
        ax.set_xlabel('mag S18')
        ax.set_ylabel('mag S16')

        """
        # logb scatter
        """
        ax = fig2.add_subplot(2,3,3)
        ax.scatter(data['logb'],data['logb16'],s=0.1)
        ax.set_xlabel('logb S18')
        ax.set_ylabel('logb S16')

        """
        # Resolution scatter
        """
        ax = fig2.add_subplot(2,3,4)
        ax.scatter(data['res'],data['res16'],s=0.1)
        ax.set_xlabel(r'$R_2$ S18')
        ax.set_ylabel(r'$R_2$ S16')

        """
        # FWHM scatter
        """
        ax = fig2.add_subplot(2,3,5)
        ax.scatter(data['fwhm'],data['fwhm16'],s=0.1)
        ax.set_xlabel('fwhm S18 (arcsec)')
        ax.set_ylabel('fwhm S16 (arcsec)')

        """
        # absE scatter
        """
        ax = fig2.add_subplot(2,3,6)
        ax.scatter(data['absE'],data['absE16'],s=0.1)
        ax.set_xlabel(r'$|e|$ S18')
        ax.set_ylabel(r'$|e|$ S16')
        """
        # Save figure
        """
        plt.tight_layout()
        fig2.savefig(os.path.join(self.config.pltDir,'compareS16_6Obs_newBriObjMsk.png'))
        return

    def finalPlot(self,sims,data):
        """
        # Plot the snr - (res,mag.|e|) Joint distributions 
        """

        """
        # Do mask
        data    =   data[~data['badPsfC']]
        maxlogb =   -0.4
        sims    =   sims[sims['logb']<maxlogb]
        mask    =   data['logb']<maxlogb
        mask    &=   data['bO_S18_r'].astype(bool)
        mask    &=  (data['absEP']<0.1)
        mask    &=  (data['fwhm']>0.47)
        """
        mask    =  (data['inS16'])
        data    =   data[mask]

        fig2 = plt.figure(figsize=(15, 5))

        """
        # The First Panel
        """
        ax = fig2.add_subplot(1,3,1)
        counts1,xbins1,ybins1,image = ax.hist2d(data['snr'], data['res'], 32, weights=data['weight'], range=[[10,80],[0.3,1]], cmap='gray_r', normed=True,vmax=0.08);
        c = ax.contour(counts1.T,extent=[10,80,0.3,1],cmap='viridis')#, levels=[0.01, 0.025, 0.05, 0.075])
        """
        counts2,xbins2,ybins2,image = ax.hist2d(sims['snr'], sims['res'], 32, weights=sims['weight'], range=[[10,80],[0.3,1]], cmap='gray', normed=True,vmax=0.08);
        c = ax.contour(counts2.T,extent=[10,80,0.3,1],cmap='viridis',linestyles='dashed')#, levels=[0.01, 0.025, 0.05, 0.075])
        """
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
        counts1,xbins1,ybins1,image = ax.hist2d(data['snr'], data['absE'], 32, weights=data['weight'], range=[[10,80],[0.,2]], cmap='gray_r', normed=True,vmax=0.045);
        c = ax.contour(counts1.T,extent=[10,80,0.,2],cmap='viridis')#, levels=[0.0065, 0.012, 0.025, 0.035])
        """
        counts2,xbins2,ybins2,image = ax.hist2d(sims['snr'], sims['absE'], 32, weights=sims['weight'], range=[[10,80],[0.,2]], cmap='gray', normed=True,vmax=0.045);
        c = ax.contour(counts2.T,extent=[10,80,0.,2], cmap='viridis', linestyles='dashed')#, levels=[0.0065,0.012,0.025,0.035])
        """
        ax.hist2d(data['snr'], data['absE'], 128, weights=data['weight'], range=[[10,80],[0.,2]], cmap='gray_r', normed=True,alpha=0.8,norm=LogNorm(vmin=0.0002,vmax=0.2));
        ax.set_xlabel('CModel SNR')
        ax.set_ylabel(r'reGaus $|e|$')
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
        counts1,xbins1,ybins1,image = ax.hist2d(data['snr'], data['mag'], 32, weights=data['weight'], range=[[10,80],[22,24.5]], cmap='gray_r', normed=True,vmax=0.1);
        c = ax.contour(counts1.T,extent=[10,80,22,24.5],cmap='viridis')#, levels=[0.01, 0.025, 0.05, 0.075])
        """
        counts2,xbins2,ybins2,image = ax.hist2d(sims['snr'], sims['mag'], 32, weights=sims['weight'], range=[[10,80],[22,24.5]], cmap='gray', normed=True,vmax=0.1);
        c = ax.contour(counts2.T,extent=[10,80,22,24.5], cmap='viridis', linestyles='dashed')#, levels=[0.01,0.025,0.05,0.075])
        """
        ax.hist2d(data['snr'], data['mag'], 128, weights=data['weight'], range=[[10,80],[22,24.5]], cmap='gray_r', normed=True,vmax=0.1,alpha=0.8,norm=LogNorm(vmin=0.005,vmax=0.1));
        ax.set_xlabel('CModel SNR')
        ax.set_ylabel('CModel mag')
        ax.set_ylim(24.5,22)
        #ax.set_xscale('log')
        ax.xaxis.set_major_formatter(mtick.FormatStrFormatter('%d'))
        ax.xaxis.set_minor_formatter(mtick.FormatStrFormatter('%d'))

        """
        # Save figure
        """
        plt.tight_layout()
        fig2.savefig(os.path.join(self.config.pltDir,'galDis2D_S16.png'))
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
