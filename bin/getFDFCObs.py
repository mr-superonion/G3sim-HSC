#!/usr/bin/env python
import os
import sys
import glob
import catSelect
import astropy.io.fits as pyfits
import sim_analysis_utilities as utilities
import numpy as np
from astropy.table import Table, vstack

def main():
    obsDir   =   '/work/xiangchong.li/work/S16AFPFS/final/noBg/field/'
    for key in utilities.field_keys: 
        fieldname   =   utilities.field_keys[key]
        obsFname    =   os.path.join(obsDir,'%s.fits' %(fieldname))
        data        =   Table.read(obsFname)
        fdfc        =   catSelect.getFDFCS16(data,degree=False)
        data['fdfc']=   fdfc
        data.write(obsFname,overwrite=True)

    return


if __name__=='__main__':
    main()
