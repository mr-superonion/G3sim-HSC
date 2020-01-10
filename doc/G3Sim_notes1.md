---
title: S18A Simulation
author: Xiangchong Li;
        HSC Collaboration
date: 2019-11

output:
  beamer_presentation:
    slide_level: 2
    toc: true
...

<link href="http://thomasf.github.io/solarized-css/solarized-light.min.css" rel="stylesheet"></link>
<style type="text/css" media="screen">#wrapper { margin: auto; max-width:1000px; } p { max-width: 600px; } </style>
<div id="wrapper">


# Image simulation

## Configuration
+ Input
    * 1600 fields (PSF and noise variance)
    * But only have 8 different noise correlation function
+ rescaling noise variance or not?
    * No rescaling:
        * [one patch](http://member.ipmu.jp/xiangchong.li/s18a_G3sim_hsm/var_varPlane_coadd.html)
        * [one tract](http://member.ipmu.jp/xiangchong.li/s18a_G3sim_hsm/var_varPlane_coadd_8524.html)
## Code for rescale_factor
pipe_tasks/6.5.1-hsc+5/python/lsst/pipe/tasks/scaleVariance.py

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


#
![mag](work/plot/galDis/magDis.pdf)

#
![SNR](work/plot/galDis/snrDis.pdf)

#
![resolution](work/plot/galDis/resolutionDis.pdf)

#
![$|e|$](work/plot/galDis/absEDis.pdf)


## $2$-D
SNR-mag, snr-res, snr-$\left | e \right |$

#
![s16a paper](ref/s16a_v2.0_v4sims_calib_overall_distributions_1d.png){ width=85% }

#
![$2$-D](work/plot/galDis/galDis2D.png)

#
![s16a paper](ref/s16a_v2.0_v4sims_calib_overall_distributions_2d.png)


# Checking Images
![Image](./outcome/gals.png){ width=50% }
![Image](./outcome/msks.png){ width=50% }

# Possible reason?
+ Blending?
    * blendedness?
    * toobig?
    * others?
+ Bad pixels (Mask plane)?
+ Bright Star Mask (Mask plane)?
+ Background Subtraction (Mask plane)?

# Matching
# S16A FDFC region
![S16A FDFC](./work/plot/galDis/galDis2D_s16aFDFC.png){ width=98% }

# Match with S16A
![S16A matched 0.5 arcsec](./work/plot/galDis/galDis2D_inS16.png){ width=98% }

# Increase Masking radius
![S16A matched 20 arcsec](./work/plot/galDis/galDis2D_inS16L.png){ width=98% }

# Spacial distribution
 <!---
 bin/spacialS16AReg.py
 -->
```python
(~inS16)
```
![0.5 arcsec](./work/plot/galDis/spacial_noInS16.png){ width=98% }

# Spacial distribution of suspect galaxies
```python
(~inS16)&(inS16L)
```
![suspect](./work/plot/galDis/spacial_suspect.png){ width=98% }

# Suspect galaxies
![Suspect galaxies](./work/plot/galDis/galDis2D_suspect_size.png){width=98% }

# Spacial distribution of extended faint galaxies
```python
(snr<15)&(res>0.9)
```
![Spacial Distribution](./work/plot/galDis/spacial_extFaiGals.png){width=98% }

# S18 Match to S16A V1
97.8% of suspect galaxies can be matched to S16A V1 (0.5 arcsec)
![](./work/plot/galDis-s18-10-12/galDis2D-snr16-res16_S16v1Mat-suspect-HECTOMAP.png){width=98% }

![](./work/plot/galDis-s18-10-12/galDis2D-snr-res_S16v1Mat-suspect-HECTOMAP.png){width=98% }

#
![](./work/plot/galDis-s18-10-12/galDis2D-snr-snr16_S16v1Mat-all-HECTOMAP.png){width=98% }

# S19 Match to S16A V1
![](./work/plot/galDis-s19-10-12/galDis2D-snr-snr16_S16v1Mat-all-HECTOMAP.png){width=98% }

# Bad pixels

#
![histogram](./work/plot/galDis/badBitHist.png)

# Pixels Masks
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

# Is adding pixFlags help?
![INEXACT_PSF_Center](./work/plot/galDis/galDis2D_inexactPSFC.png)

#
![INEXACT_PSF](./work/plot/galDis/galDis2D_inexactPSF.png)

# Bright Star Mask

# Spacial distribution of suspect galaxies
![suspect](./work/plot/galDis/spacial_suspect.png){ width=98% }

# Spacial distribution of extended faint galaxies
![Spacial Distribution](./work/plot/galDis/spacial_extFaiGals.png){width=98% }

# New Bright Star Mask in HECTOMAP
![](./work/plot/galDis-s18-10-12/galDis2D-snr-res_nomatch_newmask-std-HECTOMAP.png){width=98% }

# Blendedness

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

# Isolated galaxies ($logb < -5.0$ ) with new bright star mask
![](./work/plot/galDis-s18-10-12/galDis2D-snr-res_append_newmask-iso-HECTOMAP.png){width=75%}

# Background
# Different Background Subtraction (Check for S19)
![S19](./work/plot/galDis-s19/galDis2D_final_logbltm035_snrA10gtNop.png){width=98% }

# SNR Aperture
```python
    mask    =  (data['npass']>=2)
    mask    &= (data['logb']<=-0.4)
    mask    &= (data['snrA10']>10)
```
5% galaxies loss
![](./work/plot/galDis-s18/galDis2D_final_logbltm040_snrA10gt10_npass2.png){width=98% }

# SNR Aperture
```python
    mask    =  (data['npass']>=2)
    mask    &= (data['logb']<=-0.5)
    mask    &= (data['snrA10']>10)
```
6% galaxies loss
![](./work/plot/galDis-s18/galDis2D_final_logbltm050_snrA10gt10_npass2.png){width=98% }

# nPass
```python
    mask    =  (data['npass']>=3)
    mask    &= (data['logb']<=-0.5)
    mask    &= (data['snrA10']>10)
```
7% galaxies loss
![](./work/plot/galDis-s18/galDis2D_final_logbltm050_snrA10gt10_npass3.png){width=98% }

# snrAperture and snrlocalBackground
![](./work/plot/galDis-s18-10-12/galDis2D-logb-snrlb_append_newmask-std-HECTOMAP.png){width=75%}

# snrA10 snrlb cut
```python
mask    =  (data['npass']>=2)   #99%
mask    &= (data['logb']<=-0.375)#97%
mask    &= (data['snrlb']<0.04*data['snrA10']-0.5) #89%
```
![](./work/plot/galDis-s18-10-12/galDis2D_final_logbltm037_npass2_mag24.5_snrA10_snrlbT1.png){width=95%}



