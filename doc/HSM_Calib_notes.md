---
title: Calibration for reGauss
author: Xiangchong Li
date: 2019-11
...


# Preparation 

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

