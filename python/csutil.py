import os
import gc
import fitsio
import numpy as np
import astropy.io.ascii as ascii
import astropy.table as astTab


default_datDir='wrkDir/s19/s19-1/anaCat_newS19Mask_fdeltacut/cosmicShear/covariance-wopz/'

class cs_covariance():
    def __init__(self,datDir=default_datDir):
        self.datDir=datDir

    # pend in ['G','E0','E1','GE0']
    def getAve(self,pend):
        # Estimate the average of the correlation functions
        outFname=os.path.join(self.datDir,'cor%s_ave.csv'%pend)
        if not os.path.isfile(outFname):
            xip_ave=np.zeros(30)
            xim_ave=np.zeros(30)
            n   =   0
            for isim in range(108):
                for irot in range(13):
                    fname=os.path.join(self.datDir,'sim-%03d' %isim,'rotmat%d_cor%s.csv' %(irot,pend))
                    assert os.path.isfile(fname), print(isim,irot)
                    data=ascii.read(fname)
                    xip_ave+= data['xip']
                    xim_ave+= data['xim']
                    n+=1
                    del data
                    gc.collect()
            xip_ave=xip_ave/n
            xim_ave=xim_ave/n
            data=ascii.read(fname)
            rnom=data['R_nom']
            out=astTab.Table([rnom,xip_ave,xim_ave],names=['R_nom','xip','xim'])
            out.write(outFname)
        else:
            out=ascii.read(outFname)
            rnom    =   out['R_nom']
            xip_ave =   out['xip']
            xim_ave =   out['xim']
        return rnom,xip_ave,xim_ave

    # check the histgram of c-c_ave
    def getPlotHist(self,pend):
        xip_ave,xim_ave=getAve(pend)
        xipAll=[]
        ximAll=[]
        for isim in range(108):
            for irot in range(13):
                fname=os.path.join(self.datDir,'sim-%03d' %isim,'rotmat%d_cor%s.csv' %(irot,pend))
                assert os.path.isfile(fname), print(isim,irot)
                data=ascii.read(fname)
                xipAll.append(data['xip']-xip_ave)
                ximAll.append(data['xim']-xim_ave)
        plt.close()
        plt.hist(np.hstack(xipAll),bins=50)
        return

    def getCovMat(self,pend):
        # Estiamte the Covariance Matrix
        outFnameP=os.path.join(self.datDir,'%s%s_cov_pp.npy' %(pend,pend))
        outFnameM=os.path.join(self.datDir,'%s%s_cov_mm.npy' %(pend,pend))
        outFnameC=os.path.join(self.datDir,'%s%s_cov_pm.npy' %(pend,pend))
        if not os.path.isfile(outFnameC):
            xip_cov=np.zeros((30,30))
            xim_cov=np.zeros((30,30))
            xic_cov=np.zeros((30,30))

            xip_ave,xim_ave=getAve(pend)
            xip_ave=xip_ave.data
            xim_ave=xim_ave.data
            n=0
            for isim in range(108):
                for irot in range(13):
                    fname=os.path.join(self.datDir,'sim-%03d' %isim,'rotmat%d_cor%s.csv' %(irot,pend))
                    data=ascii.read(fname)
                    mp=np.matrix(data['xip'].data-xip_ave)
                    xip_cov=xip_cov+mp.T.dot(mp)
                    mm=np.matrix(data['xim'].data-xim_ave)
                    xim_cov=xim_cov+mm.T.dot(mm)
                    xic_cov=xic_cov+mp.T.dot(mm)
                    n+=1
                    del data
                    gc.collect()
            xip_cov=xip_cov/n
            xim_cov=xim_cov/n
            xic_cov=xic_cov/n
            np.save(outFnameP,xip_cov)
            np.save(outFnameM,xim_cov)
            np.save(outFnameC,xic_cov)
        else:
            xip_cov=np.load(outFnameP)
            xim_cov=np.load(outFnameM)
            xic_cov=np.load(outFnameC)
        return xip_cov,xim_cov,xic_cov
