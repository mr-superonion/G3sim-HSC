#!/usr/bin/env python
"""
Script to loop over simulated subfields and check the noise variance
estimated by the HSC pipeline, as compared to the true one.
"""
import os
import glob
import matplotlib
import numpy as np
matplotlib.use('agg')
matplotlib.rcParams['ps.useafm'] = True
matplotlib.rcParams['pdf.use14corefonts'] = True
matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt

import astropy.table as astTab
import lsst.afw.image as afwImg
from multiprocessing import Pool
outDir  =   './'
outfnameO   =   os.path.join(outDir,'varMeas_vs_varPlane-Obs.csv')
outfnameS   =   os.path.join(outDir,'varMeas_vs_varPlane-Sim.csv')
nprocs      =   8
names       =   ('fname','rescale_factor','var_no_rescale','var_rescale')

def processObs(outfname):
    fnames  =   glob.glob('/work/xiangchong.li/data/hsc_ssp/dr2/'\
            's18a/data/test/s18a_wide_revpsf/deepCoadd-results/HSC-I/85*/*/calexp*.fits')
    pool    =   Pool(nprocs)
    datarows=   pool.map(getPixVariance,fnames)
    datarows=   [x for x in datarows if x is not None]
    table   =   astTab.Table(rows=datarows,names=names)
    table.write(outfname,overwrite=True)
    del pool
    return table

def processSim(outfname):
    fnames  =   glob.glob('expDir/image-*.fits')[0:5]
    pool    =   Pool(nprocs)
    datarows=   pool.map(getPixVariance,fnames)
    datarows=   [x for x in datarows if x is not None]
    table   =   astTab.Table(rows=datarows,names=names)
    table.write(outfname,overwrite=True)
    del pool
    return table

def getPixVariance(fname):
    #[tract,patch]=fname.split('HSC-I')[1].split('/')[1:3]
    exposure=afwImg.ExposureF.readFits(fname)
    md  =   exposure.getMetadata()
    if 'variance_scale' in md.names():
        rescale_factor=md.get('variance_scale')
    else:
        rescale_factor=1.
    msk=exposure.getMask().getArray()[:,:]
    img=exposure.getImage().getArray()[:,:]
    var=exposure.getVariance().getArray()[:,:]
    #blank pixel
    imgBlk=img[msk==0]
    varBlk=var[msk==0]
    #normalized by std
    blkNorm=imgBlk/np.sqrt(varBlk) #unreweighted
    blkNormRescale=imgBlk/np.sqrt(varBlk/rescale_factor) #reweighted
    name=fname.split('/')[-1]
    return (name,rescale_factor,np.var(blkNorm),np.var(blkNormRescale))

def main():
    """
    if not os.path.exists(outfnameO):
        obs =   processObs(outfnameO)
    else:
        obs=astTab.Table.read(outfnameO)
    print('Observation: ave= %.2f' %(np.average(obs['var_no_rescale'])))
    """
    if not os.path.exists(outfnameS):
        sim =   processSim(outfnameS)
    else:
        sim=astTab.Table.read(outfnameS)
    print('Simulation: ave= %.2f' %(np.average(sim['var_no_rescale'])))
    return

if __name__ == '__main__':
    main()
