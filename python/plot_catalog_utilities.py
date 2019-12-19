import os
import numpy as np
import matplotlib
import astropy.io.fits as pyfits
from astropy.table import Table, join,vstack
matplotlib.rcParams['ps.useafm'] = True
matplotlib.rcParams['pdf.use14corefonts'] = True
matplotlib.rcParams['text.usetex'] = True
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import matplotlib.ticker as ticker
import matplotlib.lines as mlines

obsRanges   =   {'snr':[0,80],
                'snr16':[0,80],
                'snrA10':[0,80],
                'snrA15':[0,80],
                'snrA20':[0,80],
                'snrlb':[-0.1,0.5],
                'logb':[-6.,0],
                'logb16':[-6.,0],
                'logbO':[-6.,0],
                'logbA':[-6.,0],
                'logbR':[-6.,0],
                'res':[0.1,1.0],
                'res16':[0.1,1.0],
                'mag':[21,25.],
                'dmagAlb':[5,15.],
                'mag16':[21,25.],
                'maglb':[30,37.],
                'magA10':[21,25.],
                'fwhm':[0.1,0.8],
                'fwhm16':[0.1,0.8],
                'noivar':[0,0.007],
                'noivar16':[0,0.007],
                'absE':[0,1.],
                'absE16':[0,1.],
                'absEP':[0,1.],
                'absEP16':[0,1.],
                'detR':[0,4],
                'traceR':[0,4]}

# Hex keys for the color scheme for shear catalog paper
plot_colors = {
    'GAMA09H' : '#d73027',
    'GAMA15H' : '#fc8d59',
    'HECTOMAP' : '#fee090',
    'VVDS' : '#000000',
    'WIDE12H' : '#91bfdb',
    'XMM' : '#4575b4'
    }


def cols2DPlot(data,colns1,colns2,pltDir='./',color='gray_r',label='data',nfig=3,sims=None):
    """
    # Plot 2D Joint distributions of variables
    in colns1 and colns2
    """
    assert len(colns1)==len(colns2) and len(colns1)==nfig,\
            print('length of coln1 and colns2 must equal %d'%nfig)
    fig = plt.figure(figsize=(5*nfig, 5))
    fmt = ticker.LogFormatterMathtext()
    fmt.create_dummy_axis()
    weight  =   data['weight']
    if sims is not None:
        weight2  =   sims['weight']
    """
    # The First Panel
    """
    for i in range(nfig):
        x           =   data[colns1[i]]
        y           =   data[colns2[i]]
        if colns1[i] in obsRanges.keys():
            xMin,xMax   =   obsRanges[colns1[i]]
        else:
            print('colname %s is not in obsRanges'%colns1[i])
            xMin,xMax   =   np.percentile(x,[5,95])
        if colns2[i] in obsRanges.keys():
            yMin,yMax   =   obsRanges[colns2[i]]
        else:
            print('colname %s is not in obsRanges'%colns2[i])
            yMin,yMax   =   np.percentile(y,[5,95])
        drange  =   [[xMin,xMax],[yMin,yMax]]
        dextent =   [xMin,xMax,yMin,yMax]
        ax      =   fig.add_subplot(1,nfig,i+1)
        counts  =   np.histogram2d(x,y,32,weights=weight,range=drange)[0]
        vmax    =   np.sum(counts)/128/128*5
        #ax.clabel(c, c.levels, fmt=fmt)
        ax.hist2d(x,y,128,weights=weight,range=drange,cmap=color,vmax=vmax,alpha=0.8);
        judge   =   (colns1[i] in colns2[i]) or (colns2[i] in colns1[i])
        judge   =   judge&( ('16' in colns1[i]) or ('16'in colns2[i]) )
        if judge:
            ax.plot(np.linspace(xMin,xMax,30),np.linspace(xMin,xMax,30),color='black')
        else:
            c       =   ax.contour(counts.T,extent=dextent,cmap='viridis')
        if sims is not None:
            x2      =   sims[colns1[i]]
            y2      =   sims[colns2[i]]
            counts  =   np.histogram2d(x2,y2,32,weights=weight2,range=drange)[0]
            ax.contour(counts2.T,dextent,cmap='viridis',linestyles='dashed')


        if 'mag' in colns1[i]:
            ax.invert_xaxis()
        if 'mag' in colns2[i]:
            ax.invert_yaxis()
        ax.grid()
        ax.set_xlabel(colns1[i])
        ax.set_ylabel(colns2[i])
        if i==1:
            """
            #  Adding legend
            """
            data_line = mlines.Line2D([], [], color='black', ls='-', \
                    label=label.replace('_',' '))
            ax.legend(handles=[data_line])
    plt.tight_layout()

    """
    # Save figure
    """
    outFname=   os.path.join(pltDir, 'galDis2D-%s-%s_%s.png' %(colns1[0],colns2[0],label))
    fig.savefig(outFname)
    return
