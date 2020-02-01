#!/usr/bin/env python
import os
import numpy as np
import astropy.table as astTable

import matplotlib as mpl
mpl.use('agg')
mpl.rcParams['ps.useafm'] = True
mpl.rcParams['pdf.use14corefonts'] = True
mpl.rcParams['text.usetex'] = True

import matplotlib.pyplot as plt

obsDir  =   'database/'
if not os.path.isdir(obsDir):
    raise runtimeerror("required directory %s does not exist or is not a directory" % (anadir))
version =   'ratio'

def getCounts(data,apertureMask=False):
    mask    =   (data['fieldname']==field)
    mask    &=   data['inS16FDFC']
    mask    &=  (data['logb']<=-0.375)
    if apertureMask:
        mask    &=  (data['magA10']<25.5)
    dataP   =   data[mask]
    return np.histogram2d(dataP['ra'],dataP['dec'],[xGrid,yGrid])[0]

def int_product(x, pos):
    """The two args are the value and tick position.
    Label ticks with the product of the exponentiation"""
    return '%d' % (x)
def float_product(x, pos):
    """The two args are the value and tick position.
    Label ticks with the product of the exponentiation"""
    return '%.1f' % (x)
if 'ratio' not in version:
    norm    =   mpl.colors.Normalize(vmin=0,vmax=40)
    formatter = plt.FuncFormatter(int_product)
else:
    norm    =   mpl.colors.Normalize(vmin=0.7,vmax=1.2)
    formatter = plt.FuncFormatter(float_product)

if 'ratio' not in version:
    datFname=   os.path.join(obsDir,'allObs-%s.fits' %(version))
    if not os.path.exists(datFname):
        raise runtimeerror("required catalog %s does not exist" % (datFname))
    data    =   astTable.Table.read(datFname)
else:
    datFname=   os.path.join(obsDir,'allObs-s18.fits')
    if not os.path.exists(datFname):
        raise runtimeerror("required catalog %s does not exist" % (datFname))
    data1   =   astTable.Table.read(datFname)

    datFname=   os.path.join(obsDir,'allObs-s16.fits')
    if not os.path.exists(datFname):
        raise runtimeerror("required catalog %s does not exist" % (datFname))
    data2   =   astTable.Table.read(datFname)


pltDir  =   's18-2/plot/spacial/'
if not os.path.isdir(pltDir):
    os.system('mkdir -p %s'%pltDir)

nMapDir =   os.path.join(os.environ['homeWrk'],'S16AFPFS/out/numMap/')
# Prepare Figure
fig     =   plt.figure(figsize=(15, 9))
# Prepare axies for different fields
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

# Plot for different fields
hDat=[]
for i,field in enumerate(fields):
    print('Plotting field: %s' %field)
    ax      =   axes[i]
    ax.set_title(field)
    ref     =   np.load('%s/%sNumMap.npy' %(nMapDir,field))
    xGrid   =   ref[1]
    yGrid   =   ref[2]

    if 'ratio' not in version:
        counts  =   getCounts(data,version=='s18')
    else:
        counts1 =   getCounts(data1,True)
        counts2 =   getCounts(data2,False)
        counts  =   np.zeros(counts1.shape)
        mask    =   (counts2>0)
        counts[mask]= counts1[mask]/counts2[mask]
    hDat.extend(counts[counts>0])
    im      =   ax.imshow(counts.T,interpolation='None',extent=[xGrid[0],xGrid[-1],yGrid[-1],yGrid[0]],norm=norm)
    x0,x1 = ax.get_xlim();ax.set_xlim(x1, x0)
    y0,y1 = ax.get_ylim();ax.set_ylim(y1, y0)

cbar_ax     =   fig.add_axes([0.95, 0.1, 0.01, 0.8])
cbar        =   fig.colorbar(im, cax=cbar_ax, format=formatter)
plt.subplots_adjust(wspace=0.01, hspace=-0.03,left=0.05,right=0.95,bottom=0.05,top=0.95)
fig.savefig(os.path.join(pltDir,'%s_numMap_magA10lt255.png' %(version)))

fig=    plt.figure(figsize=(6, 6))
ax =    fig.add_subplot(1,1,1)
ax.hist(hDat, 32, range=(0.02,1.5))
fig.savefig(os.path.join(pltDir,'%s_hist_magA10lt255.png' %(version)))
