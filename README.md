# s18G3

## from DB to catalog
+ ./bin/hscReleaseQuery.py > ./work/database/s18/tracts/
> query data from DB and download to tracts

+ ./bin/makeGalCatFromDBRun.py > ./work/database/s18/tracts2/
> Delete is_null
> FDFC Mask
> Standard WL Mask
> Comprize to Fields
 

## Image simulation
+ ./bin/rgcSimBatchRun.py: simulate rgc exposures
+ ./bin/processSimBatchRun.py: process simulated exposures


## Shear Calibration for reGauss


