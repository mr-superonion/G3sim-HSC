#!/usr/bin/env python
import os
import numpy as np
import astropy.table as astTable

import matplotlib
matplotlib.use('TkAgg')
matplotlib.rcParams['ps.useafm'] = True
matplotlib.rcParams['pdf.use14corefonts'] = True
matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt

def main():
    binAll=np.load('./exfaGal-src/binAll.npy')
    binAll2=np.load('./exBrGal-src/binAll.npy')
    pltDir      =   'plot/galDis/'
    fig =   plt.figure()
    ax  =   fig.add_subplot(111)
    ax.hist(binAll,rwidth=0.5,align='left',normed=True,bins=np.arange(17)-0.5,label='exfaint')
    ax.hist(binAll2,rwidth=0.5,align='mid',normed=True,bins=np.arange(17)-0.5,label='normal')
    plt.xlabel('bit')
    plt.ylabel('num')
    plt.tight_layout()
    plt.legend()
    plt.savefig(os.path.join(pltDir,'badBitHist.png'))
    return

if __name__=='__main__':
    main()
