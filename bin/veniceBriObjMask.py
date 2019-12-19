#!/usr/bin/env python
import os
import argparse
import sim_analysis_utilities as utilities
regName =   '/work/xiangchong.li/work/S18ASimHSC/s18/brightObjMask/HSC-SSP_brightStarMask_Arcturus/reg/masks_all.reg'

def s18BriObjMask(inCat,outCat,regName,raName='i_ra',decName='i_dec'):
    bName   =   '/work/xiangchong.li/work/S18ASimHSC/s18/brightObjMask/HSC-SSP_brightStarMask_Arcturus/venice-4.0.3/bin/venice'
    flagName=   'mask_briObj_S18'
    command =   '%s -m %s -cat %s -xcol %s -ycol %s -f all -flagName %s -o %s' %(bName,regName,inCat,raName,decName,flagName,outCat)
    print(command)
    os.system(command)
    return

def main(inDir,outDir):
    for fn in utilities.field_names.items():
        inCat   =   os.path.join(inDir,'%s.fits'%fn[1] )
        outCat  =   os.path.join(outDir,'%s.fits'%fn[1])
        utilities.s18BriObjMask(inCat,outCat,regName)
    return


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--inDir',help='directory of input files')
    parser.add_argument('--outDir',help='directory of output files')
    args = parser.parse_args()
    main(args.inDir,args.outDir)
