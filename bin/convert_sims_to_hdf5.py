#!/usr/bin/env python
# This script converts the HSC fits data files into one HDF-5 file for faster
# access where only columns of interest can be selectively loaded
# Plus, it now also aggregate the noise rescaling factors for each (tract, patch)
import sys
import os
import optparse
from astropy.table import Table, vstack
from sim_analysis_utilities import *

def main(argv):
    """
    Loads the fits simulation files and saves them into and HDF5 fie
    """
    usage = "usage: %prog [options]"
    description = "Exports the GREAT3-like sims files into single HDF5 file."
    parser = optparse.OptionParser(usage=usage, description=description)

    parser.add_option("--sims_dir", dest="sims_dir", type=str,
                      default=default_sim_dir,
                      help=("Directory with sims catalogs; "
                            "default=%s"%default_sim_dir))

    parser.add_option("--sims_version", dest="sims_version", type=str,
                      default=default_sims_version,
                      help=("Version of the GREAT3-like simulations; "
                            "default=%s"%default_sims_version))

    parser.add_option("--rescale_noise", dest="rescale_noise", type=int,
                      default=0,
                      help=("Rescale misestimated noise variances?; default=0"))

    # Parse command line options
    opt, args = parser.parse_args()

    # Create dictionary to store the catalogs per fields
    catalogs = {}

    if opt.rescale_noise==0:
        filename = os.path.join('src_by_fields.hdf5')
    else:
        filename = os.path.join('src_by_fields_rescale.hdf5')
    if os.path.exists(filename):
        raise RuntimeError("File you are trying to write already exists: %s"%filename)

    for i in range(3200):
        print("Processing subfield %d"%i)
        try:
            catalog  = prepare_sim_catalog(opt.sims_dir,i,rescale_est_variance=opt.rescale_noise)
        except IOError:
            # Technically there could be some other IOError here,
            # e.g., missing variance file rather than this subfield
            # being missing.  But for now let's keep it like this.
            print("Subfield %d missing, skipping to the next one."%i)
            continue
        catalog.meta = None
        fieldname = catalog['field_name'][0]

        if fieldname not in catalogs:
            catalogs[fieldname] = []

        catalogs[fieldname].append(catalog)
    print("Data loaded, exporting to hdf5")
    for fieldname in catalogs:
        catalog =   vstack(catalogs[fieldname])
        catalog.write(filename,path='/%s'%(fieldname),append=True)
if __name__ == "__main__":
    main(sys.argv)
