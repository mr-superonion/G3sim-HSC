---
title: S18A Simulation
author: Xiangchong Li;
        HSC Collaboration
date: 2019-11
...


# Image simulation
## Configuration
+ Input
    * 1600 fields (PSF and noise variance)
    * But only have 8 different noise correlation function
+ rescaling noise variance or not?
    * No rescaling:
        * [one patch](http://member.ipmu.jp/xiangchong.li/s18a_G3sim_hsm/var_varPlane_coadd.html)
        * [one tract](http://member.ipmu.jp/xiangchong.li/s18a_G3sim_hsm/var_varPlane_coadd_8524.html)

# Processing
+ hscpipe 6.7
+ magnitude zero point
```python
magnitude   =   
-2.5 * np.log10(catalog['modelfit_CModel_flux']) + 27.
```


# Analysis Preparation 

## Read the catalog data
/work/xiangchong.li/work/S18ASimHSC/s18

## Apply mask for permissive cuts
+ Extendness>0
+ Deblend Nchild==0
+ centDist < 5

## Add Columns   
+ $e_1$, $e_2$ etc.
+ psfInfo,catInfo,matched with ipos



#
## Identify pairs based on their input ipos values
```python
n = 100*100 # Number of grid points
mask_even = np.in1d(catalog["ipos"] % n, 
    2*np.arange(n/2).astype("int"), assume_unique=True)
mask_odd  = np.in1d(catalog["ipos"] % n, 
    2*np.arange(n/2).astype("int")+1, assume_unique=True)

broken_pairs = np.setxor1d(catalog["ipos"][mask_even] % n,
    catalog["ipos"][mask_odd] % n - 1)
mask_pairs = ~np.in1d(catalog["ipos"] % n,
    np.concatenate([broken_pairs, broken_pairs+1]), 
    assume_unique=True)
catalog["paired"] = mask_pairs
```

# Match PSF and Noise
+ Determine the variance distribution and FWHM distribution of data and sim
+ Determine the ratio between data and sim
+ Reweight sim to match FWHM and Variance to the data

#
![variance reweight](work/plot/psfNoi/noiVar2.pdf)

#
![FWHM reweight](work/plot/psfNoi/psfFWHM2.pdf)



# Check galaxy number distribution

## $1$-D
Four observables: magnitude, SNR, resolution, $\left | e \right |$

## $2$-D
SNR-mag, snr-res, snr-$\left | e \right |$

#
![mag](work/plot/galDis/magDis.pdf)

#
![SNR](work/plot/galDis/snrDis.pdf)

#
![resolution](work/plot/galDis/resolutionDis.pdf)

#
![$|e|$](work/plot/galDis/absEDis.pdf)


#
![s16a paper](ref/s16a_v2.0_v4sims_calib_overall_distributions_1d.png){ height=85% }

#
![$2$-D](work/plot/galDis/galDis2D.png)

#
![s16a paper](ref/s16a_v2.0_v4sims_calib_overall_distributions_2d.png)


# Checking Images
![Image](./outcome/gals.png){ width=50% }
![Image](./outcome/msks.png){ width=50% }

+ Blending?
+ Bad pixels (Mask plane)?

# Definiation of blendedness_abs_flux has Changed?
![](./ref/blendedness_cut.png){ width=50% }
![](./ref/blendedness_define_change1.png){ width=50% }

# Blendedness vs. Resolution
![](./work/plot/galDis/logb-res.png){ height=85% }

# Number loss with blendedness cut
![](./work/plot/galDis/frac-logb.png){ height=85% }


# After blendedness cut..
![$logb<-0.4$](./work/plot/galDis/galDis2D_logb_lt_m040.png)


# After blendedness cut..
![$logb<-0.6$](./work/plot/galDis/galDis2D_logb_lt_m060.png)


# After blendedness cut..
![$logb<-0.8$](./work/plot/galDis/galDis2D_logb_lt_m080.png)


# After blendedness cut..
![$logb<-1.0$](./work/plot/galDis/galDis2D_logb_lt_m100.png)

# Bad pixels
BAD 0
SAT 1
INTRP 2
CR 3
EDGE 4

DETECTED 5
DETECTED_NEGATIVE 6
SUSPECT 7
NO_DATA 8

BRIGHT_OBJECT 9
CLIPPED 10
CROSSTALK 11
INEXACT_PSF 12

NOT_DEBLENDED 13
REJECTED 14
SENSOR_EDGE 15
UNMASKEDNAN 16

# Bad pixels
![histogram](./work/plot/galDis/badBitHist.png)


# Match with S16A
## Six Observables
![compare 6 observables](./work/plot/galDis/compareS16_6Obs.png)

#
## Spacial distribution
 <!---
 bin/compareS16Spacial.py  
 -->
![Spacial Distribution](./work/plot/galDis/compareS16_spacial.png){ width=98% }

