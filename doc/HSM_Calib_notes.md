---
title: Calibration for reGauss
author: Xiangchong Li
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
![$2$-D](work/plot/galDis/galDis2D.pdf)

#
![s16a paper](ref/s16a_v2.0_v4sims_calib_overall_distributions_2d.png)
