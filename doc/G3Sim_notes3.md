---
title: S18A Simulation
author: Xiangchong Li;
        HSC Collaboration
date: 2020-01

output:
  beamer_presentation:
    slide_level: 2
    toc: true
...

<link href="http://thomasf.github.io/solarized-css/solarized-light.min.css" rel="stylesheet"></link>
<style type="text/css" media="screen">#wrapper { margin: auto; max-width:1000px; } p { max-width: 600px; } </style>
<div id="wrapper">

# s16 vs. s17 vs. s18
![](./work/s18-2/plot/galDis-s16-s17-s18-append.png){width=75%}


---------------------------------------------------

# exp_niter and dev_niter
![](./work/s18-2/plot/galDis-s18-10-12/galDis2D-bsf-magA10_final-std-s16-append.png)

----------------------------------------------------

# exp_niter and dev_niter cut
```python
(data['exp_niter']<14)&(data['dev_niter']<20)
```
16% galaxy lost
![](./work/s18-2/plot/galDis-s18-10-12-final/galDis2D_final_logbltm037_npass2_mag24.5_niter14-20.png)

----------------------------------------------------

```python
(data['exp_niter']<16)&(data['dev_niter']<25)
```
11% galaxy lost
![](./work/s18-2/plot/galDis-s18-10-12-final/galDis2D_final_logbltm037_npass2_mag24.5_niter16-25.png)

----------------------------------------------------

```python
data['magA10']<25.5
```
8% galaxy lost
![](./work/s18-2/plot/galDis-s18-10-12-final/galDis2D_final_logbltm037_npass2_mag24.5_magA10-lt25.5.png)


-----------------------------------------------------
