#!/usr/bin/env python
import os
import numpy as np
from astropy.table import Table,vstack
import sim_analysis_utilities as utilities

import matplotlib
matplotlib.use('agg')
matplotlib.rcParams['ps.useafm'] = True
matplotlib.rcParams['pdf.use14corefonts'] = True
matplotlib.rcParams['text.usetex'] = True
import plot_catalog_utilities as pltUtil

import matplotlib.pyplot as plt

nbin    =   64
cc      =   'Blues'
Dir     =   '/work/xiangchong.li/work/S18ASimHSC/prepare2/'
fname   =   os.path.join(Dir,'fieldTractInfo.csv')
fname1  =   os.path.join(Dir,'info.csv')
fname2  =   os.path.join(Dir,'info2.csv')
fname3  =   os.path.join(Dir,'info3.csv')

pltDir  =   os.path.join(Dir,'plot3')
if not os.path.isdir(pltDir):
    os.system('mkdir -p %s' %pltDir)

def getInfo(fname):
    tractList   =   []
    fieldList   =   []
    while True:
        field   =   input('please input field name: \n')
        if not len(field)>0:
            break
        while True:
            tractMin=   input('please input tractMin: \n')
            tractMax=   input('please input tractMax: \n')
            if not len(tractMin)>0:
                break
            tractMin=   eval(tractMin)
            tractMax=   eval(tractMax)
            for i in range(tractMin,tractMax+1):
                tractList.append(i)
                fieldList.append(field)

    tableOut=   Table([tractList,fieldList],names=['tract','field'])
    tableOut.write(fname)
    return tableOut

def goodSample(ti,fd=None):
    if fd!=None:
        mask0   =   ti['field']==fd
    else:
        mask0   =   np.ones(len(ti),dtype=bool)
    ind     =   np.where(mask0)[0]

    ###
    to      =   ti[mask0]
    vMin    =   np.percentile(to['varPlane'],1)
    vMax    =   np.percentile(to['varPlane'],99)
    dv      =   (vMax-vMin)/nbin
    fMin    =   np.percentile(to['fwhm'],1)
    fMax    =   np.percentile(to['fwhm'],99)
    df      =   (fMax-fMin)/nbin
    ###

    nBin    =   np.int_((ti['varPlane']-vMin)/dv)
    pBin    =   np.int_((ti['fwhm']-fMin)/df)
    xBin    =   nBin*100+pBin

    indU    =   np.unique(xBin,return_index=True)[1]
    return np.intersect1d(ind,indU)

def main():
    if not os.path.exists(fname3):
        if not os.path.exists(fname2):
            if not os.path.exists(fname):
                t   =   getInfo(fname)
            else:
                t   =   Table.read(fname)
            to  =   Table.read(fname1)
            to['field']=None
            for row in to:
                row['field']    =   t[t['tract']==row['tract']]['field'][0]
            to.write(fname2)
            del t
        else:
            to=Table.read(fname2)

        to['field2']=   to['field2'].astype('object')
        to['field'] =   to['field'].astype('object')
        to['fBin']  =   np.zeros(len(to),dtype=int)
        for row in to:
            row['field'] =  utilities.field_keys[row['field2']]
            row['fBin']=    eval(row['field2'].split('-')[0][-1])
        t=to
        t.write(fname2)

        nnn     =   int(len(to)/3000)
        ind     =   np.arange(0,len(to),nnn)

        for fd in utilities.field_keys.values():
            ind =   np.unique(np.hstack([ind,goodSample(to,fd)]))
        ind =   np.unique(np.hstack([ind,goodSample(to)]))

        to      =   to[ind]
        to.write(fname3)
    else:
        to= Table.read(fname3)
        t = Table.read(fname2)

    tryList  =  ['all', 'HECTOMAP', 'VVDS',
            'GAMA09H', 'GAMA15H', 'WIDE12H',
            'XMM']

    for itry in tryList:
        if itry != 'all':
            tU  =   t[(t['field'] == itry)]
            toU =   to[(to['field'] == itry)]
        else:
            tU  =   t
            toU =   to
        l1='origin-%s' %itry
        l2='select-%s' %itry

        pltUtil.cols2DPlot(tU,['varPlane'],['fwhm'],pltDir=pltDir,color=cc,label=l1,nfig=1,ngridH=16,ngridC=16,writeFits=True)
        pltUtil.cols2DPlot(toU,['varPlane'],['fwhm'],pltDir=pltDir,color=cc,label=l2,nfig=1,ngridH=16,ngridC=16,writeFits=True)
    return

if __name__=="__main__":
     main()
