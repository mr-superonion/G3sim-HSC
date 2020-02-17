#!/usr/bin/env python
import os
import fitsio
import smatch
from astropy.table import Table,vstack
import numpy as np
import sim_analysis_utilities as utilities

import matplotlib
matplotlib.rcParams['ps.useafm'] = True
matplotlib.rcParams['pdf.use14corefonts'] = True
matplotlib.rcParams['text.usetex'] = True
matplotlib.use('agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

import lsst.afw.image as afwImg
import lsst.afw.detection as afwDet
import lsst.afw.geom as afwGeom
import lsst.afw.coord as afwCoord
import lsst.afw.table as afwTable

from lsst.meas.base.noiseReplacer import NoiseReplacerConfig,NoiseReplacer

noiRep  =   False
inSrc   =   False
ngrid   =   64
npg     =   500
outDir  =   'exfaGalAll'
if noiRep:
    outDir+='-noiRep'
if inSrc:
    outDir+='-srcInser'

expname =   '%s/%s/calexp-HSC-I-%s-%s.fits'
srcname =   '%s/%s/meas-HSC-I-%s-%s.fits'
bnames  =   ['TPid', 'ra', 'dec', 'object_id', 'traceR16', \
            'traceR', 'traceRP', 'mag', 'mag16', 'xPix', 'yPix']

cmap    =   plt.cm.Paired
norm    =   matplotlib.colors.Normalize(vmin=21,vmax=25.)
mappable=   plt.cm.ScalarMappable(norm=norm, cmap=cmap)

def readData(version):
    fname   =   'database/compare/match-%s.fits' %version
    names   =   ['TPid', 'TPid16', 'ra', 'dec', 'snr', 'snr16', \
                'object_id', 'object_id16', 'traceR', 'traceR16', \
                'res', 'res16', 'mag', 'mag16', 'fieldname']
    data    =   Table.read(fname)
    # Do mask
    mask    =   (~np.isnan(data['mag16']))
    mask    &=  (~np.isnan(data['mag']))
    mask    &=  (~np.isnan(data['traceR16']))
    mask    &=  (~np.isnan(data['traceR']))
    # WL cut
    mask    &=  (data['mag']<24.5)
    mask    &=  (data['res']>0.3)
    mask    &=  (data['snr']>10.)
    mask    &=  (data['npass']>=2)
    mask    &=  (data['logb']<=-0.375)
    # Bad cut
    mask    &=  (data['mag16']>24.5)
    mask    &=  (data['magA10']>25.5)
    mask    &=  ((data['traceR']-data['traceR16'])>0.2)
    #
    mask    &=  (data['TPid']==data['TPid16'])
    return data[mask][names]

def prepare():
    fname =   'database/compare/susMatch-s16-s18-s19.fits'
    if not os.path.exists(fname):
        print('Preparing matched catalog')
        data18  =   readData('s18')
        data19  =   readData('s19')

        match   =   smatch.match(data18['ra'],data18['dec'],0.4/3600.,\
            data19['ra'],data19['dec'],maxmatch=1)
        mask    =   data18[match['i1']]['TPid']==data19[match['i2']]['TPid']
        match   =   match[mask]
        index   =   np.unique(match['i1'],return_index=True)[1]
        match   =   match[index]

        data18  =   data18[match['i1']]
        data19  =   data19[match['i2']]
        data18['object_id19']=  data19['object_id']
        data18['snr19']   =   data19['snr']
        data18['mag19']   =   data19['mag']
        data18['res19']   =   data19['res']
        data18['traceR19']=   data19['traceR']
        data18['TPid19']  =   data19['TPid']
        data18.write(fname)
    else:
        print('Already have matched catalog')
        data18  =  Table.read(fname)
    return data18

def cutStampFromExposure(centroid,src,exposure):
    # Prepare
    imDat   =   exposure.getImage().getArray()
    mskDat  =   exposure.getMask().getArray()
    expBox  =   exposure.getBBox()
    psfExp  =   exposure.getPsf()

    centx   =   int(centroid.getX()+0.5-expBox.getBeginX())
    centy   =   int(centroid.getY()+0.5-expBox.getBeginY())
    rCut    =   min(centx,centy,expBox.getWidth()-1-centx,\
                    expBox.getHeight()-1-centy,ngrid//2)
    rCut    =   int(rCut)
    if rCut<2:
        return None
    #get image
    left    =   centx-rCut
    right   =   centx+rCut
    down    =   centy-rCut
    up      =   centy+rCut
    #image, mask and PSF
    galData =   np.zeros((ngrid,ngrid))
    mskData =   np.zeros((ngrid,ngrid))
    psfData =   np.zeros((ngrid,ngrid))
    shift   =   int(ngrid//2-rCut)
    num     =   0
    for j in range(2*rCut):
        for i in range(2*rCut):
            galData[shift+j,shift+i]=   imDat[down+j,left+i]
            mskData[shift+j,shift+i]=   mskDat[down+j,left+i]
            num                     |=  mskDat[down+j,left+i]

    binary  =   []
    idig    =   0
    while num:
        bit     =   num % 2
        if bit:
            binary.append(idig)
        num     =   num // 2
        idig    +=  1
    afwGeom.Point2D(centroid)
    psfData2    =   psfExp.computeImage(afwGeom.Point2D(centroid)).getArray()
    ny          =   psfData2.shape[0]
    nx          =   psfData2.shape[1]
    shifty      =   (ny-ngrid-1)//2
    shiftx      =   (nx-ngrid-1)//2
    if shiftx>=0:
        for j in range(ngrid):
            for i in range(ngrid):
                    psfData[j,i]=psfData2[shifty+j,shiftx+i]
    else:
        shiftx=abs(shiftx)
        shifty=abs(shifty)
        for j in range(ny):
            for i in range(nx):
                    psfData[j+shifty,i+shiftx]=psfData2[j,i]
    return galData,mskData,psfData,binary

def process(data,version):
    nGal    =   len(data)
    if version=='s16':
        tractDirE=  os.path.join(os.environ['s16w'],'deepCoadd/HSC-I/')
        tractDirS=  os.path.join(os.environ['s16w'],'deepCoadd-results/HSC-I/')
        idname  =   'object_id16'
    elif version=='s18':
        tractDirE=  os.path.join(os.environ['s18w'],'deepCoadd-results/HSC-I/')
        tractDirS=  os.path.join(os.environ['s18w'],'deepCoadd-results/HSC-I/')
        idname  =   'object_id'
    elif version=='s19':
        tractDirE=  os.path.join(os.environ['s19w'],'deepCoadd-results/HSC-I/')
        tractDirS=  os.path.join(os.environ['s19w'],'deepCoadd-results/HSC-I/')
        idname  =   'object_id19'
    Dir     =   os.path.join(outDir,version)
    if not os.path.exists(Dir):
        os.system('mkdir -p %s' %Dir)

    igal        =   0
    igroup      =   0
    mskDataA    =   np.zeros((npg,ngrid,ngrid))
    galDataA    =   np.zeros((npg,ngrid,ngrid))
    psfDataA    =   np.zeros((npg,ngrid,ngrid))
    catA        =   []
    tpSet       =   list(set(data['TPid']))
    print(tpSet)
    for tp in tpSet:
        dataU   =   data[data['TPid']==tp]
        # Read Exposure and catalog
        tract   =   tp//1000
        patch   =   '%d,%d' %(tp%1000//100,tp%1000%100)
        # Exposure
        expFname=   os.path.join(tractDirE,expname %(tract,patch,tract,patch))
        exposure=   afwImg.ExposureF.readFits(expFname)
        # Sources
        srcFname=   os.path.join(tractDirS,srcname %(tract,patch,tract,patch))
        sources =   afwTable.SourceCatalog.readFits(srcFname)

        if noiRep:
            # noise replacement
            footprints  = {measRecord.getId(): (measRecord.getParent(), measRecord.getFootprint())
                          for measRecord in sources}
            noiConfig   =   NoiseReplacerConfig()
            noiseReplacer=  NoiseReplacer(noiConfig, exposure, footprints)

        for ss in dataU:
            # Find source
            src =   sources.find(ss[idname])
            if inSrc and noiRep:
                noiseReplacer.insertSource(src.getId())
            #pixP2D= afwGeom.Point2I(ss['xPix'],ss['yPix'])
            pixP2D= afwGeom.Point2I(src.getCentroid())
            galData,mskData,psfData,binary =   cutStampFromExposure(pixP2D,src,exposure)
            if inSrc and noiRep:
                noiseReplacer.removeSource(src.getId())

            mskDataA[igal%npg]=mskData
            galDataA[igal%npg]=galData
            psfDataA[igal%npg]=psfData
            catA.append(ss)
            igal+=1
            if igal%npg==0 or igal==nGal:
                outFnameM=   os.path.join(Dir,'msk-%d.fits' %igroup)
                outFnameG=   os.path.join(Dir,'gal-%d.fits' %igroup)
                outFnameP=   os.path.join(Dir,'psf-%d.fits' %igroup)
                outFnameC=   os.path.join(Dir,'cat-%d.csv'  %igroup)
                fitsio.write(outFnameM,mskDataA)
                """
                fitsio.write(outFnameG,galDataA)
                fitsio.write(outFnameP,psfDataA)
                vstack(catA).write(outFnameC)
                """
                igroup +=1
                mskDataA    =   np.zeros((npg,ngrid,ngrid))
                galDataA    =   np.zeros((npg,ngrid,ngrid))
                psfDataA    =   np.zeros((npg,ngrid,ngrid))
                catA        =   []
    return

def makePlot(gal,cat):
    nl      =   int(np.sqrt(len(cat)))
    if nl<=1:
        return
    figname =   os.path.join(outDir,'fig-%d.png' %tp)
    fig     =   plt.figure(figsize=(3*nl+2, 3*nl+2))
    ax      =   fig.add_subplot(1,1,1)
    bdry    =   np.percentile(gal,np.arange(0,101,10))
    norm2   =   matplotlib.colors.BoundaryNorm(bdry,ncolors=256)
    im      =   ax.imshow(gal,norm=norm2,cmap='gray')
    for ic,cc in enumerate(cat):
        ix  =   ic%nl
        iy  =   ic//nl
        rg16=   cc['traceR16']/0.168
        rg  =   cc['traceR']/0.168
        rg  =   min(rg,ngrid//2)
        rg16=   min(rg16,ngrid//2)
        clr =   cmap(norm(cc['mag']))
        clr16=  cmap(norm(cc['mag16']))
        xy  =   (ix*ngrid+ngrid//2,iy*ngrid+ngrid//2)
        cirG=   mpatches.Circle(xy, radius=rg, ls='-', \
                facecolor='None', edgecolor=clr,lw=3.)
        cirP=   mpatches.Circle(xy, radius=rg16, ls='--', \
                facecolor='None', edgecolor=clr16, lw=2.)
        ax.add_patch(cirG)
        ax.add_patch(cirP)

    mappable.set_array(cat['mag'])
    fig.colorbar(im,orientation='vertical')
    fig.colorbar(mappable=mappable,orientation='horizontal')
    plt.tight_layout()
    fig.savefig(figname)
    plt.close('all')
    del fig
    return

def main():
    dataAll  =   prepare()
    dataAll  =   dataAll[dataAll['fieldname']=='HECTOMAP']
    process(dataAll,'s19')
    return

if __name__=='__main__':
    main()
