#!/usr/bin/env python
import os
import itertools
import numpy as np
import astropy.io.fits as pyfits
import astropy.io.ascii as ascii
from multiprocessing import Pool
from astropy.table import Table


def histoGram(ifile):
    srcFname=   './outcomeFPFS/meas%04d.fits' %(ifile)
    if not os.path.exists(srcFname):
        return
    data    =   pyfits.getdata(srcFname)
    mask    =   np.ones(len(data)).astype(bool)
    mask    =   mask&(~np.isnan(data['base_SdssShape_psf_xx']))
    mask    =   mask&(~np.isnan(data['base_SdssShape_psf_yy']))
    mask    =   mask&(~np.isnan(data['base_Variance_value']))
    data    =   data[mask]
    rpsfs   =   np.sqrt((data['base_SdssShape_psf_xx']+data['base_SdssShape_psf_yy'])/2.)
    rpsf    =   rpsfs[0]
    var     =   data['base_Variance_value'][0]
    return (ifile,var,rpsf)

if __name__ == '__main__':
    ndat    =   1600
    files   =   range(ndat)
    p       =   Pool(8)

    datarows=   p.map(histoGram,files)
    datarows=   [x for x in datarows if x is not None]
    print(datarows)
    names   =   ('field','var','rpsf')
    t = Table(rows=datarows, names=names,dtype=('i4','f8','f8'))
    t.write('noiVar_rPsf_s18a.fits')

