# s18G3
For summaries of the work. see ./doc/%s.md

## From DB to catalog
+ ./bin/hscReleaseQuery.py > ./work/database/s18/tracts/
> hscReleaseQuery.py --user lixiangchong@stars $sqlPATH/hscS18WL.sql --format fits --release-version dr2 --nomail --delete-job
> query data from DB and download to groups
> unzip groupd into tracts

+ ./bin/makeGalCatFromDBRun.py > ./work/database/s18/fields/
> Delete is_null
> FDFC Mask
> Standard WL Mask
> Comprize to Fields


## Test the galaxies
+ ./bin/galDisAnalysisRun.py  -> ./work/anaCat/allObs%s.fits 
+ ./bin/galDisAnalysis2DRun.py
+ ./bin/galDisAnalysis2DSizeRun.py

## Image simulation
+ ./bin/rgcSimBatchRun.py: simulate rgc exposures
+ ./bin/processSimBatchRun.py: process simulated exposures

## Match noise variance and PSF size

## Shear Calibration for reGauss

