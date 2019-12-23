#!/usr/bin/env python
import os
import fitsio
import astropy.io.ascii as ascii
import numpy as np
import sim_analysis_utilities as utilities

import matplotlib
matplotlib.rcParams['ps.useafm'] = True
matplotlib.rcParams['pdf.use14corefonts'] = True
matplotlib.rcParams['text.usetex'] = True
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

import lsst.afw.image as afwImg
import lsst.afw.detection as afwDet
import lsst.afw.geom as afwGeom
import lsst.afw.coord as afwCoord
import lsst.afw.table as afwTable

from lsst.meas.base.noiseReplacer import NoiseReplacerConfig,NoiseReplacer

neiRep  =   True
ngrid   =   32
baseDir =   'exfaGal-src3'
if neiRep:
    baseDir+='-neiRep'
version =   's18-10-12'
obsFname=   'anaCat/allObs-%s-%s.fits'
tractDir    =   os.path.join(os.environ['s18w'],'deepCoadd-results/HSC-I/')
bnames  =   ['TPid', 'ra', 'dec', 'object_id', \
            'traceR', 'traceRP', 'mag']
cmap    =   plt.cm.Paired
norm    =   matplotlib.colors.Normalize(vmin=21,vmax=25.)
mappable=   plt.cm.ScalarMappable(norm=norm, cmap=cmap)

def getObsGalData(fname):
    data    =   fitsio.read(fname)
    mask    =   data['magA10']>25.5
    mask    |=  (data['snrlb']>0.01*(data['snrA10']))
    mask    &=  data['inS16FDFC']
    mask    &=  (data['npass']>=2)
    mask    &=  (data['logb']<=-0.375)
    mask    &=  (~data['inS16'])
    return data[mask][bnames]

def cutStampFromExposure(centroid,ngrid,expBox,imDat,mskDat,psfExp):
    centx       =   int(centroid.getX()+0.5-expBox.getBeginX())
    centy       =   int(centroid.getY()+0.5-expBox.getBeginY())
    rCut        =   min(centx,centy,expBox.getWidth()-1-centx,\
                    expBox.getHeight()-1-centy,ngrid//2)
    rCut        =   int(rCut)
    if rCut<2:
        return None
    #get image
    left        =   centx-rCut
    right       =   centx+rCut
    down        =   centy-rCut
    up          =   centy+rCut
    #image, mask and PSF
    galData     =   np.zeros((ngrid,ngrid))
    mskData     =   np.zeros((ngrid,ngrid))
    psfData     =   np.zeros((ngrid,ngrid))
    shift       =   int(ngrid//2-rCut)
    num         =   0
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
    afwGeom.coordinates.Point2D(centroid)
    psfData2    =   psfExp.computeImage(afwGeom.coordinates.Point2D(centroid)).getArray()
    shift       =   (psfData2.shape[0]-ngrid-1)//2
    for j in range(ngrid):
        for i in range(ngrid):
            psfData[j,i]=psfData2[shift+j,shift+i]
    return galData,mskData,psfData,binary

def process(tp,data,outDir):
    outFnameG=   os.path.join(outDir,'gal-%d.fits' %tp)
    outFnameM=   os.path.join(outDir,'msk-%d.fits' %tp)
    outFnameP=   os.path.join(outDir,'psf-%d.fits' %tp)
    outFnameC=   os.path.join(outDir,'cat-%d.csv' %tp)
    if os.path.exists(outFnameC):
        return fitsio.read(outFnameG),fitsio.read(outFnameM),\
        fitsio.read(outFnameP),ascii.read(outFnameC)

    tract   =   tp//1000
    patch   =   '%d,%d' %(tp%1000//100,tp%1000%100)

    # Exposure
    expFname=   os.path.join(tractDir,'%s/%s/calexp-HSC-I-%s-%s.fits'%(tract,patch,tract,patch))
    exposure=   afwImg.ExposureF.readFits(expFname)
    imDat   =   exposure.getImage().getArray()
    mskDat  =   exposure.getMask().getArray()
    bbox    =   exposure.getBBox()
    psfExp  =   exposure.getPsf()
    wcs     =   exposure.getWcs()

    # Sources
    srcFname=   os.path.join(tractDir,'%s/%s/meas-HSC-I-%s-%s.fits'%(tract,patch,tract,patch))
    sources =   afwTable.SourceCatalog.readFits(srcFname)
    maskU   =   data['TPid']==tp
    dataU   =   data[maskU]

    # Prepare Output data
    ndata   =   len(dataU)
    nl      =   min(int(np.sqrt(ndata)),8)
    nl2     =   nl**2
    galDataA=   np.zeros((nl*ngrid,nl*ngrid))
    mskDataA=   np.zeros((nl*ngrid,nl*ngrid))
    psfDataA=   np.zeros((nl*ngrid,nl*ngrid))
    dataU   =   dataU[0:nl2]

    if neiRep:
        # noise replacement
        footprints  = {measRecord.getId(): (measRecord.getParent(), measRecord.getFootprint())
                              for measRecord in sources}
        noiConfig   =   NoiseReplacerConfig()
        noiseReplacer=  NoiseReplacer(noiConfig, exposure, footprints)

    for i in range(nl2):
        ss  =   dataU[i]
        # Make Sky coordinate
        p2D     =   afwGeom.coordinates.Point2D(ss['ra'],ss['dec'])
        crd2D   =   afwCoord.IcrsCoord(p2D,afwGeom.degrees)
        pixP2D  =   afwGeom.coordinates.Point2I(wcs.skyToPixel(crd2D))
        # Find source
        src     =   sources.find(ss['object_id'])
        if neiRep:
            noiseReplacer.insertSource(src.getId())
        pixP2D2 =   afwGeom.coordinates.Point2I(src.getCentroid())
        galData,mskData,psfData,binary =   cutStampFromExposure(pixP2D2,ngrid,bbox,imDat,mskDat,psfExp)
        # The output data
        ix  =   i%nl
        iy  =   i//nl
        xb  =   ix*ngrid
        yb  =   iy*ngrid
        xe  =   (ix+1)*ngrid
        ye  =   (iy+1)*ngrid
        galDataA[yb:ye,xb:xe]=galData
        mskDataA[yb:ye,xb:xe]=mskData
        psfDataA[yb:ye,xb:xe]=psfData
        if neiRep:
            noiseReplacer.removeSource(src.getId())
    fitsio.write(outFnameG,galDataA)
    fitsio.write(outFnameM,mskDataA)
    fitsio.write(outFnameP,psfDataA)
    ascii.write(dataU,outFnameC)
    return galDataA,mskDataA,psfDataA,dataU

def main():
    fieldnames  =   ['HECTOMAP']#list(utilities.field_names.keys())
    for fn in fieldnames:
        outDir      =   '%s/%s' %(baseDir,fn)
        if not os.path.exists(outDir):
            os.system('mkdir -p %s' %outDir)
        data        =   getObsGalData(obsFname %(version,fn))
        tpSet       =   list(set(data['TPid']))
        for tp in tpSet:
            gal,msk,psf,cat =   process(tp,data,outDir)
            nl      =   int(np.sqrt(len(cat)))
            if nl<=1:
                continue
            figname =   os.path.join(outDir,'fig-%d.png' %tp)
            fig     =   plt.figure(figsize=(3*nl+2, 3*nl+2))
            ax      =   fig.add_subplot(1,1,1)
            bdry    =   np.percentile(gal,np.arange(0,101,10))
            norm2   =   matplotlib.colors.BoundaryNorm(bdry,ncolors=256)
            im      =   ax.imshow(gal,norm=norm2,cmap='gray')
            for ic,cc in enumerate(cat):
                ix  =   ic%nl
                iy  =   ic//nl
                rp  =   cc['traceRP']/0.168
                rg  =   cc['traceR']/0.168
                rg  =   min(rg,ngrid//2)
                clr =   cmap(norm(cc['mag']))
                xy  =   (ix*ngrid+ngrid//2,iy*ngrid+ngrid//2)
                cirG=   mpatches.Circle(xy, radius=rg, ls='-', \
                        facecolor='None', edgecolor=clr,lw=3.)
                cirP=   mpatches.Circle(xy, radius=rp, ls='--', \
                        facecolor='None', edgecolor=clr, lw=2.)
                ax.add_patch(cirG)
                ax.add_patch(cirP)

            mappable.set_array(cat['mag'])
            fig.colorbar(im,orientation='vertical')
            fig.colorbar(mappable=mappable,orientation='horizontal')
            plt.tight_layout()
            fig.savefig(figname)
            del fig
    return

if __name__=='__main__':
    main()
