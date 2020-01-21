#!/usr/bin/env python
import sys
import os
from astropy.table import Table, vstack
import sim_analysis_utilities as utilities

def main():
    for fieldname in utilities.field_keys.values():
        run(fieldname)
    return

def run(fieldname):
    anaDir  =   './anaCat/'
    print('Processing field: %s' %(fieldname))
    fname   =   os.path.join(anaDir,fieldname,'%s.fits'%(fieldname))
    data    =   Table.read(fname)
    outDir  =   os.path.join(anaDir,fieldname,'simCat')
    if not os.path.exists(outDir):
        os.mkdir(outDir)

    simList =   Table.read('info3.csv')
    numList =   simList[simList['field']==fieldname]['num']
    for isim in numList:
        outname =   os.path.join(outDir,'%s.fits' %(isim))
        dataSub =   data[data['isim']==isim]
        dataSub.write(outname)
        del dataSub
    return

if __name__=='__main__':
    main()
