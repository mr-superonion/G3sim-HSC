#!/usr/bin/env python
import os
import lsst.afw.image as afwImage
import lsst.daf.persistence as dafPersist

exp=afwImage.ExposureF('./hsc_ssp/dr3/s19a/data/s19a_wide/deepCoadd-results/HSC-I\
        /9314/1,1/calexp-HSC-I-9314-1,1.fits')
info=exp.getInfo()
ccds=info.getCoaddInputs().ccds

dataDirs=[]

addStr='wget --user=XXX --password =XXX \
    --page-requisites --span-hosts \
    -R "*index.html*" -R "*warp*" -r -np -nH -e    \
    robots=off -nc https://hscdata.mtk.nao.ac.jp/'

ff=open('download.sh','w+')

ccd in ccds:
    visit=ccd.get('visit')
    chip=ccd.get('ccd')
    dataId={'visit':visit,'ccd':chip}
    dataDir=butler.get('calexp_filename',dataId)[0].split('li/data/')[1]
    out=addStr+dataDir+'\n'
    ff.write(out)

ff.close()
