#!/usr/bin/env python
import os
import numpy as np
import astropy.io.fits as pyfits


"""
#Note
This subroutine is only for FPFS and not for reGauss
"""

def prepareBWeight(ngrid,wFname):
    Bweights=np.zeros((26,ngrid,ngrid))
    for rCut in range(8,32):
        weightArray =   np.zeros((ngrid,ngrid))
        for j in range(ngrid):
            for i in range(ngrid):
                ic  =   i-ngrid/2.
                jc  =   j-ngrid/2.
                rc  =   np.sqrt(ic**2.+jc**2.)
                if rc<rCut:
                    weightArray[j,i]=1.
        weightArray =   np.fft.fftshift(np.fft.ifft2(abs(np.fft.fft2(weightArray))**2.).real) 
        Bweights[rCut-8,:,:]    =   weightArray        
    np.save(wFname,Bweights)
    return

def main():
    #Load the weight caused by boundary
    outFname=   './corPre/noiPows.npy'
    ngrid   =   64
    inWname =   './corPre/boundaryWeight.npy'
    if not os.path.exists(inWname):
        prepareBWeight(ngrid,inWname)
    weights =   np.load(inWname)
    corNoise=   pyfits.getdata('./corPre/correlation.fits')
    ngridC  =   corNoise.shape[0]
    print(ngridC)
    shift   =   int((ngrid-ngridC)//2)+1
    weiCor0 =   np.zeros((ngrid,ngrid))
    weiCor0[shift:shift+ngridC,shift:shift+ngridC]=corNoise[:,:]
    print(np.unravel_index(np.argmax(weiCor0),(64,64)))
    outPows =   np.zeros((24,1,ngrid,ngrid))
    for irad in range(8,32):
        iweiCor =   weiCor0*weights[irad-8,:,:]
        inoiPow =   np.fft.fft2(np.fft.ifftshift(iweiCor)).real
        inoiPow =   np.fft.fftshift(inoiPow)
        itodiv  =   np.sqrt(sum(sum(inoiPow**2.)))
        inoiPow =   inoiPow/itodiv
        outPows[irad-8,0,:,:]=inoiPow
    np.save(outFname,outPows)
    return
    
if __name__ == '__main__':
    main()
