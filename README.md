# s18G3

## Download catalogs from DB
+ bin/hscReleaseQuery.py
> query data from DB and download to groups

```shell
hscReleaseQuery.py --user lixiangchong@stars $sqlPATH/hscS18WL.sql --format fits --release-version dr2 --nomail --delete-job
```

+ ./bin/makeGalCatFromDBRun.py
> Delete is_null
> FDFC Mask
> Standard WL Mask

## Image simulation
+ ./bin/rgcSimBatchRun.py: simulate rgc exposures
+ ./bin/processSimBatchRun.py: process simulated exposures

## Determine weak-lensing sample
+ ./bin/galDisAnalysisRun.py
+ ./bin/galDisAnalysis2DRun.py
+ ./bin/galDisAnalysis2DSizeRun.py

## reGauss catalog

### Optimization

### Calibration

