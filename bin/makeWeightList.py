#!/usr/bin/env python
import os
import itertools
import numpy as np
import astropy.io.fits as pyfits
import astropy.io.ascii as ascii
from multiprocessing import Pool
from astropy.table import Table
import lsst.afw.image as afwImg


if __name__ == '__main__':
    ngrid1=25
    ngrid2=20
    histVar=np.linspace(0.1,1.2,ngrid1)
    histR=np.linspace(.8,2.2,ngrid2)
    dvar    =   histVar[1]-histVar[0]
    dR      =   histR[1]-histR[0]
    data    =   Table.read('noiVar_rPsf_s18a.fits')
    ratio2D =   pyfits.getdata('ratio2D.fits')
    weightList=[]
    for idat in data:
        varBin=int((idat['var']*100-0.1)//dvar)
        rBin=int((idat['rpsf']-0.8)//dR)
        if rBin>=0 and rBin<19 and varBin>=0 and varBin<24:
            iratio  =   ratio2D[varBin,rBin]
            weightList.append(iratio)
        else:
            weightList.append(0.)
    t = Table([data['field'],weightList], names=('field','weight'))
    t.write('weightList.fits')


