---
title: S18A Simulation
author: Xiangchong Li;
        HSC Collaboration
date: 2019-12

output:
  beamer_presentation:
    slide_level: 2
    toc: true
...

<link href="http://thomasf.github.io/solarized-css/solarized-light.min.css" rel="stylesheet"></link>
<style type="text/css" media="screen">#wrapper { margin: auto; max-width:1000px; } p { max-width: 600px; } </style>
<div id="wrapper">

# Suspect vs. Good
![](./work/plot/galDis-s18-10-12-HECTOMAP/galDis2D-magA10-dmagPC_S16v1Mat-append-HECTOMAP.png)

# local background
## DB code
```python
Flux_pipe/10**(-2/5*(48.6+27.0)) #  => (erg s^{-1} cm^{-2} Hz^{-1})
```

## pipeline code
lsst/meas/base/LocalBackground.h
```C
float const innerRadius = _ctrl.annulusInner*psfSigma;
afw::geom::ellipses::Axes const innerCircle{innerRadius, innerRadius};
auto const& inner = afw::geom::SpanSet::fromShape(afw::geom::ellipses::Ellipse(innerCircle, center));

float const outerRadius = _ctrl.annulusOuter*psfSigma;
afw::geom::ellipses::Axes const outerCircle{outerRadius, outerRadius};
auto const& outer = afw::geom::SpanSet::fromShape(afw::geom::ellipses::Ellipse(outerCircle, center));

auto const& annulus = outer->clippedTo(image.getBBox())->intersectNot(*inner);
auto const& imageValues = annulus->flatten(image.getImage()->getArray(), image.getXY0());
auto const& maskValues = annulus->flatten(image.getMask()->getArray(), image.getXY0());
```
------------------

## pipeline code (continue)

```C
badMaskPlanes({"BAD", "SAT", "NO_DATA"}),
annulusInner(7.0),
annulusOuter(15.0),
bgRej(3.0),
bgIter(3)
```
![](./ref/localBackground.png){ height=30% }

-------------------

## Image without noise replacement
![](./work/exfaGal-src3/HECTOMAP/fig-15827006.png){ height=95% }

------------------

## Image with noise replacement
![](./work/exfaGal-src3-neiRep/HECTOMAP/fig-15827006.png){ height=95% }

------------------

## Code for the noise replacement
```python

footprints  = {measRecord.getId(): (measRecord.getParent(), measRecord.getFootprint()) for measRecord in sources}
noiConfig   = NoiseReplacerConfig()
noiseReplacer=  NoiseReplacer(noiConfig, exposure, footprints)

noiseReplacer.insertSource(src.getId())

#Extract Images...

noiseReplacer.removeSource(src.getId())

```

# New Selection Criteria 1
```python
mask    &= (data['snrlb']<0.01*(data['snrA10']))#91%
mask    &=  data['magA10']<25.5                 #89%
```
![](./work/plot/galDis-s18-10-12/galDis2D_final_logbltm037_npass2_mag24.5_snrA10_snrlbT3.png){width=95%}

-------------------

# New Selection Criteria 2
```python
mask    =  (data['npass']>=2)   #99%
mask    &= (data['logb']<=-0.375)#97%
mask    &= (data['snrlb']<0.04*data['snrA10']-0.5) #89%
```
![](./work/plot/galDis-s18-10-12/galDis2D_final_logbltm037_npass2_mag24.5_snrA10_snrlbT2.png){width=95%}


-------------------

# Bright Star Mask (new vs old vs s16A)
![](./work/plot/spacial/spacial_unmatch05_inNewMask.png)


-------------------
![](./work/plot/spacial/spacial_unmatch05_outNewMask.png)

