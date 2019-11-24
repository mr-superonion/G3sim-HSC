#!/usr/bin/env python
import os
import fitsio
import numpy as np
from astropy.table import Table, vstack
import sim_analysis_utilities as utilities

import lsst.afw.image as afwImg
import lsst.afw.detection as afwDet
import lsst.afw.geom as afwGeom
import lsst.afw.coord as afwCoord

def getObsGalData(obsFname):
    data    =   Table.read(obsFname)
    res     =   utilities.get_res(data)
    snr     =   utilities.get_snr(data)
    mask    =   (res>0.9)&(snr<15.)
    return data[mask][['tract','patch','i_ra','i_dec','object_id']]

def cutStampFromExposure(centroid,ngrid,expBox,imDat,mskDat):
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
    shift       =   int(ngrid//2-rCut)
    for j in range(2*rCut):
        for i in range(2*rCut):
            galData[shift+j,shift+i]=   imDat[down+j,left+i]
            mskData[shift+j,shift+i]=   mskDat[down+j,left+i]
    return galData,mskData

def main():
    fieldnames  =   list(utilities.field_names.keys())
    for fn in fieldnames:
        process(fn)
    return

def process(fieldname):
    ngrid       =   32
    outDir      =   'exfaGal/%s' %(fieldname)
    if not os.path.exists(outDir):
        os.mkdir(outDir)
    obsDir      =   '/gpfs02/work/miyatake/S18AL/with_new_star_mask_updated20191117/'
    tractDir    =   os.path.join(os.environ['s18w'],'deepCoadd-results/HSC-I/')
    obsFname    =   os.path.join(obsDir,'%s.fits' %(utilities.field_names[fieldname]))
    data        =   getObsGalData(obsFname)
    TPid        =   data['tract']*1000+data['patch']
    data['TPid']=  TPid
    data.write('exfaGal/%s.csv' %fieldname)
    tpName      =   os.path.join(outDir,'TPset.npy') 
    if not os.path.exists(tpName):
        tpSet   =   list(set(TPid))
        np.save(tpName,tpSet)
    else:
        print('already have TP file')
        tpSet   =   np.load(tpName) 
    for tp in tpSet[0:20]:
        tract   =   tp//1000
        patch   =   '%d,%d' %(tp%1000//100,tp%1000%100)
        expFname=   os.path.join(tractDir,'%s/%s/calexp-HSC-I-%s-%s.fits'%(tract,patch,tract,patch))
        exposure=   afwImg.ExposureF.readFits(expFname)
        imDat   =   exposure.getImage().getArray()
        mskDat  =   exposure.getMask().getArray() 
        bbox    =   exposure.getBBox()
        dataU   =   data[TPid==tp]
        wcs     =   exposure.getWcs()
        for ss in dataU:
            # Make Sky coordinate
            p2D     =   afwGeom.coordinates.Point2D(ss['i_ra'],ss['i_dec'])
            crd2D   =   afwCoord.IcrsCoord(p2D,afwGeom.degrees)
            pixP2D  =   afwGeom.coordinates.Point2I(wcs.skyToPixel(crd2D))
            galData,mskData =   cutStampFromExposure(pixP2D,ngrid,bbox,imDat,mskDat)
            outFname1=   os.path.join(outDir,'gal-%d-%d.fits' %(ss['TPid'],ss['object_id']))
            fitsio.write(outFname1,galData)
            outFname2=   os.path.join(outDir,'msk-%d-%d.fits' %(ss['TPid'],ss['object_id']))
            fitsio.write(outFname2,mskData)
    return

if __name__=='__main__':
    main()
