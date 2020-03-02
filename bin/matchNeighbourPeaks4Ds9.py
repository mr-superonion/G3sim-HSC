#!/usr/bin/env python
import os
import smatch
import numpy as np
from astropy.table import Table, vstack
import astropy.io.fits as pyfits
import sim_analysis_utilities as utilities
import plot_catalog_utilities as pltUtil

version =   's19'
obsDir  =   'database/compare/'
if not os.path.isdir(obsDir):
    raise RuntimeError("Data catalog directory %s does not exist "
         "or is not a directory" % (obsDir))

def main():
    # Read data
    name1   =   'exfaGalAll/%s/cat-0.csv' %version
    cat1    =   Table.read(name1)
    name2   =   'database/compare/scat-%s.fits' %version
    cat2    =   Table.read(name2)
    # Do match
    match   =   smatch.match(cat1['ra'],cat1['dec'],4./3600.,\
                cat2['ra'],cat2['dec'],maxmatch=-1)
    cat1    =   cat1[match['i1']]
    cat2    =   cat2[match['i2']]
    cat2['id1']=match['i1']
    cat2['id2']=match['i2']
    cat2['dx']=-(cat2['ra']-cat1['ra'])/0.168*3600.
    cat2['dy']=(cat2['dec']-cat1['dec'])/0.168*3600.
    nameO   =   'exfaGalAll/%s/mat-0.csv' %version
    cat2.write(nameO,overwrite=True)
    return

if __name__=='__main__':
    main()

