#!/usr/bin/env python
import os
import smatch
import astropy.table as astTable

mskDir  =   '/work/xiangchong.li/work/S18ASimHSC/s18/brightObjMask/newMask'
if not os.path.isdir(mskDir):
    raise runtimeerror("required directory %s does not exist or is not a directory" % (mskDir))

anaDir  =   'anaCat'
if not os.path.isdir(anaDir):
    raise runtimeerror("required directory %s does not exist or is not a directory" % (anadir))

def main():
    galFname=   os.path.join(anaDir,'allObs-s18-10-12-HECTOMAP.fits')
    catG    =   astTable.Table.read(galFname)
    raG     =   catG['ra']
    decG    =   catG['dec']

    starFname=   os.path.join(mskDir,'hectomap_20191211.csv')
    catS    =   astTable.Table.read(starFname)
    raS     =   catS['ra']
    decS    =   catS['dec']
    rS      =   catS['radius']

    print(len(catS))
    match   =   smatch.match(raS,decS,rS,\
            raG,decG,maxmatch=1)
    print(len(match['i1']))
    print(len(set(match['i1'])))
    print(len(set(match['i2'])))
    catG['mask_briObj_s18'] =   False
    catG['mask_briObj_s18'][match['i2']]=True

    outFname=   os.path.join(anaDir,'allObs-s18-10-12-HECTOMAP-newMask.fits')
    catG.write(outFname)
    return

if __name__=='__main__':
    main()
