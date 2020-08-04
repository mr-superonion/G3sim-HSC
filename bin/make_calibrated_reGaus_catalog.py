#!/usr/bin/env python
import gc
import os
import numpy as np
from astropy.table import Table, vstack
import sim_analysis_utilities as utilities

outDir  =   './catalog_obs_reGaus_calibrated/fields/'
if not os.path.exists(outDir):
    os.system('mkdir -p %s' %outDir)

obsDir  =   '../../database/s19-Hironao-2/fields/'
if not os.path.isdir(obsDir):
    raise RuntimeError("Data catalog directory %s does not exist "
        "or is not a directory" % (obsDir))


def main():
    for fieldname in utilities.field_names:
        print('processing field: %s' %fieldname)
        obsFname=   os.path.join(obsDir,'%s.fits' %(fieldname))
        zFname  =   os.path.join(obsDir,'%s_pz.fits' %(fieldname))
        outFname=   os.path.join(outDir,'%s_calibrated.fits' %(fieldname))
        data    =   Table.read(obsFname)
        zdat    =   Table.read(zFname)
        data['mizuki_photoz_best']      =   zdat['mizuki_photoz_best']
        data['mizuki_photoz_conf_best'] =   zdat['mizuki_photoz_conf_best']
        data['mizuki_photoz_conf_best'] =   zdat['mizuki_photoz_conf_best']
        data['mizuki_photoz_err95_min'] =   zdat['mizuki_photoz_err95_min']
        data['mizuki_photoz_err95_max'] =   zdat['mizuki_photoz_err95_max']
        data    =   utilities.update_wl_flags(data)
        data    =   data[data['weak_lensing_flag']]
        data    =   utilities.update_reGaus_calibration(data)
        data.write(outFname)
        del data
        del zdat
        gc.collect()
    return

if __name__=='__main__':
    main()
