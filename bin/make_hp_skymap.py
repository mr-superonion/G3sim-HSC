#!/usr/bin/env python
import fitsio
import numpy as np
import healpy as hp

nside=32
npix=hp.nside2npix(nside)
pixAll=np.arange(npix)
vec=hp.boundaries(nside=nside,pix=pixAll,nest=True)
vec=np.swapaxes(vec,1,2)
raBs,decBs=hp.vec2ang(vec,lonlat=True)
raBs[raBs<1e-6]=360.

raBs=raBs.reshape((npix,4))
decBs=decBs.reshape((npix,4))
# neglect data from poles
mask=(decBs[:,0]<85)&(decBs[:,2]>-85)
raBs=raBs[mask]
decBs=decBs[mask]
pixAll=pixAll[mask]
# the healpix cross ra=0.
maskMax=np.abs(raBs[:,3]-np.max(raBs,axis=1))>1e-3
raBs[:,3][maskMax]=raBs[:,3][maskMax]+360.
# maximum and minimum
raMax=raBs[:,3]+0.02;raMin=raBs[:,1]-0.02
decMax=decBs[:,0]+0.02;decMin=decBs[:,2]-0.02
raC=(raMin+raMax)/2.
decC=(decMin+decMax)/2.
dra=raMax-raMin;ddec=decMax-decMin

types=[('pix','i4'),('raMin','>f8'),('raMax','>f8'),('decMin','>f8'),('decMax','>f8'),('dra','>f8'),('ddec','>f8')]
data=np.stack([pixAll,raMin,raMax,decMin,decMax,dra,ddec]).T
data=[tuple(ii) for ii in data]
out=np.array(data,dtype=types)
fitsio.write('/work/xiangchong.li/work/skyMap/healpix-nside%d-nest.fits' %nside,out)
