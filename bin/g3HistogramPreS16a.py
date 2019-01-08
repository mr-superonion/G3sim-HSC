#!/usr/bin/env python
import os
import itertools
import numpy as np
import astropy.io.fits as pyfits
import astropy.io.ascii as ascii
from multiprocessing import Pool
from astropy.table import Table
import lsst.afw.image as afwImg


def histoGram(ifile):
    preFileName =   '/work/xiangchong.li/work/S16ASimHSC/sample4/epoch_catalog-%04d-0.fits' %(ifile)
    if not os.path.exists(preFileName):
        return
    preDat  =   pyfits.getdata(preFileName)
    var     =   preDat['var'][0]
    srcFname=   '/work/xiangchong.li/work/S16ASimHSC/sample4/src-%04d.fits' %(ifile)
    if not os.path.exists(srcFname):
        return
    src     =   pyfits.getdata(srcFname)
    shapePsfs=   src['shape_sdss_psf']
    mask    =   (~np.isnan(shapePsfs[:,0]))
    shapePsfs=  shapePsfs[mask]
    shapePsf=   shapePsfs[0]
    rpsf    =   np.sqrt((shapePsf[0]+shapePsf[1])/2.)
    return (ifile,var,rpsf)

if __name__ == '__main__':
    ndat    =   800
    files   =   range(ndat)
    p       =   Pool(8)

    datarows=   p.map(histoGram,files)
    datarows=   [x for x in datarows if x is not None]
    names   =   ('field','var','rpsf')
    t = Table(rows=datarows, names=names,dtype=('i4','f8','f8'))
    t.write('noiVar_rPsf_s16a.fits')

