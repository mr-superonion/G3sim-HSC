#!/usr/bin/env python
import sys
import os
from astropy.table import Table, vstack
import sim_analysis_utilities as utilities


def main():
    for fieldKey in utilities.field_keys:
        run(fieldKey)
    return

    
def run(fieldKey):
    if fieldKey=='aegis':
        return
    anaDir  =   './anaCat/'
    fieldname   =   utilities.field_keys[fieldKey]
    print('Processing field: %s' %(fieldname))
    fname   =   os.path.join(anaDir,fieldname,'%s.fits'%(fieldname))
    data    =   Table.read(fname)
    outDir  =   os.path.join(anaDir,fieldname,'simCat')
    if not os.path.exists(outDir):
        os.mkdir(outDir)
    simList =   Table.read('noiVar_rPsf_s18a.csv')
    numList =   simList[simList['field']==fieldKey]['num']
    for isim in numList:
        outname =   os.path.join(outDir,'%s.fits' %(isim))
        dataSub =   data[data['isim']==isim]
        dataSub.write(outname)
        del dataSub

    return


if __name__=='__main__':
    main()

