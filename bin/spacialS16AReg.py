#!/usr/bin/env python
import os
import numpy as np
import astropy.table as astTable

import matplotlib as mpl
mpl.use('TkAgg')
mpl.rcParams['ps.useafm'] = True
mpl.rcParams['pdf.use14corefonts'] = True
mpl.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt

anaDir  =   'anaCat'
nMapDir =   os.path.join(os.environ['homeWrk'],'S16AFPFS/out/numMap/')
data    =   astTable.Table.read(os.path.join(anaDir,'allObs4.fits'))

def log_10_product(x, pos):
    """The two args are the value and tick position.
    Label ticks with the product of the exponentiation"""
    return '%d' % (x)
formatter = plt.FuncFormatter(log_10_product)
norm=mpl.colors.Normalize(vmin=0,vmax=20.)

fig     =   plt.figure(figsize=(15, 9))
fields  =   ['HECTOMAP','VVDS','WIDE12H','GAMA15H','GAMA09H','XMM']
axes    =   []
for i in range(3):
    f1=fields[i*2]
    f2=fields[i*2+1]
    xlen1=len(np.load('%s/%sNumMap.npy' %(nMapDir,f1))[1])
    xlen2=len(np.load('%s/%sNumMap.npy' %(nMapDir,f2))[1])
    ratio=int((xlen1/(xlen1+xlen2+0.))*100+0.5)
    axes.append(plt.subplot2grid((3, 106), (i, 0), colspan=ratio))
    axes.append(plt.subplot2grid((3, 106), (i, ratio+4), colspan=100-ratio))

for i,field in enumerate(fields):
    print('Plotting for field: %s' %field)
    ax      =   axes[i]
    ref     =   np.load('%s/%sNumMap.npy' %(nMapDir,field))
    xGrid   =   ref[1]
    yGrid   =   ref[2]
    dataU   =   data[data['fieldname']==field]
    """
    mask    =   dataU['inS16']
    mask    &=  dataU['bO_S18_r'].astype(bool)
    data16  =   dataU[mask]
    res     =   dataU['res']
    snr     =   dataU['snr']
    mask    =   (res>0.8)&(res<1.)&(snr>10.)&(snr<15.)
    data16N =   dataU[mask]
    #counts,xbins,ybins,im       =   ax.hist2d(data16['ra'],data16['dec'],[xGrid,yGrid])
    """
    mask0   =   dataU['inS16FDFC']
    mask    =   mask0&(~dataU['inS16'])
    mask    =   mask&(dataU['inS16L'])
    data16N =   dataU[mask]
    counts2,xbins2,ybins2   =   np.histogram2d(data16N['ra'],data16N['dec'],[xGrid,yGrid])
    im2     =   ax.imshow(counts2.T,interpolation='None',extent=[xGrid[0],xGrid[-1],yGrid[-1],yGrid[0]],norm=norm)
    #ax.contour(counts2.T,extent=[xGrid[0],xGrid[-1],yGrid[0],yGrid[-1]])
    ax.set_title(field)
    x0,x1 = ax.get_xlim()
    y0,y1 = ax.get_ylim()
    ax.set_xlim(x1, x0)
    ax.set_ylim(y1, y0)

cbar_ax     =   fig.add_axes([0.95, 0.1, 0.01, 0.8])
cbar        =   fig.colorbar(im2, cax=cbar_ax,format=formatter)#,ticks=[0, 5, 10,15,20])
plt.subplots_adjust(wspace=0.01, hspace=-0.03,left=0.05,right=0.95,bottom=0.05,top=0.95)
pltDir  =   'plot/galDis/'
fig.savefig(os.path.join(pltDir,'spacial_suspect.png'))
