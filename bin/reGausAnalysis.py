#!/usr/bin/env python

import os
import sys
import optparse
import sim_analysis_utilities as utilities
import numpy as np
from sim_analysis_utilities import default_sim_dir, default_data_dir
from  sim_analysis_utilities import default_sims_version, default_data_version
from base_analysis import get_analysis_arguments, get_analysis_description
from reGaus_sigmaE import sigmae_analysis

# Removed AEGIS since the simulations do not contain this field anyway
data_fields = ['gama9', 'gama15', 'HECTOMAP', 'VVDS', 'wide12', 'XMM']

# Put each one on separate line so we can quickly comment out ones we don't want
# when we are testing.
analyses_list = [
    sigmae_analysis,
    ]

default_n_subfields = 6
default_plot_dir = "plots/"
default_plot_pref = ""
default_plot_type = "eps"

# Default lensing cuts
default_use_wl_flag = 1
default_bright_star_mask = 1
default_reweight_s16a = 1
default_rescale_est_variance = 0
#default_rescale_est_variance = 1
default_use_calib = 0

def replace_calibs(catalog, data_version, sims_version, data_dir):
    """
    This utility replaces the calibration corrections in the baseline
    catalogs with those from utilities.py.  The catalog itself is
    updated in place.
    """
    catalog['ishape_hsm_regauss_derived_sigma_e'] = \
        utilities.get_sigma_e_model(catalog, data_version=data_version, sims_version=sims_version)
    catalog['ishape_hsm_regauss_derived_rms_e'] = \
        utilities.get_erms_model(catalog, data_version=data_version, sims_version=sims_version)
    catalog['ishape_hsm_regauss_derived_shape_weight'] = \
        utilities.get_weight_model(catalog, data_version=data_version, sims_version=sims_version)

    # Fix types before assignment.  Apparently these entries were
    # added as ints so we have to fix this before trying to replace
    # the columns directly.  This won't be a problem if we rerun the
    # hdf5 catalog builder.
    catalog.replace_column('ishape_hsm_regauss_derived_shear_bias_m',
                           catalog['ishape_hsm_regauss_derived_shear_bias_m'].astype(float))
    catalog.replace_column('ishape_hsm_regauss_derived_shear_bias_c1',
                           catalog['ishape_hsm_regauss_derived_shear_bias_c1'].astype(float))
    catalog.replace_column('ishape_hsm_regauss_derived_shear_bias_c2',
                           catalog['ishape_hsm_regauss_derived_shear_bias_c2'].astype(float))
    catalog['ishape_hsm_regauss_derived_shear_bias_m'] = \
        utilities.get_m_model(catalog, data_version=data_version,
                              sims_version=sims_version,
                              dir=data_dir)
    model_a = utilities.get_a_model(catalog, data_version=data_version, sims_version=sims_version)
    psf_e1, psf_e2 = utilities.get_psf_ellip(catalog)
    catalog['ishape_hsm_regauss_derived_shear_bias_c1'] = model_a*psf_e1
    catalog['ishape_hsm_regauss_derived_shear_bias_c2'] = model_a*psf_e2
    return

def main(argv):
    """
    Main function to run the suite of analysis scripts or a subset of them.
    """
    usage = "usage: %prog [options] [analysis name]"
    description = ("Runs a suite of validation analyses on the "
                    "great3_like simulations.")
    parser = optparse.OptionParser(usage=usage, description=description)

    # General options not specific to a given analysis
    parser.add_option("--sims_dir", dest="sims_dir", type=str,
                      default=default_sim_dir,
                      help=("Directory with simulation catalogs; "
                            "default=%s"%default_sim_dir))

    parser.add_option("--data_dir", dest="data_dir", type=str,
                      default=default_data_dir,
                      help=("Directory with data catalogs; "
                            "default=%s"%default_data_dir))

    parser.add_option("--data_version", dest="data_version", type=str,
                      default=default_data_version,
                      help=("Version of the HSC data; "
                            "default=%s"%default_data_version))

    parser.add_option("--sims_version", dest="sims_version", type=str,
                      default=default_sims_version,
                      help=("Version of the simulations; "
                            "default=%s"%default_sims_version))

    parser.add_option("--n_subfield", dest="n_subfield", type=int,
                      default=default_n_subfields,
                      help=("Number of subfields to process; "
                            "default=%d"%default_n_subfields))

    parser.add_option("--plot_dir", dest="plot_dir", type=str,
                      default=default_plot_dir,
                      help=("Path to use when saving plots; "
                            "default=%s"%default_plot_dir))

    parser.add_option("--plot_pref", dest="plot_pref", type=str,
                      default=default_plot_pref,
                      help=("Prefix to use when saving plots; "
                            "default=%s"%default_plot_pref))

    parser.add_option("--plot_type", dest="plot_type", type=str,
                      default=default_plot_type,
                      help=("Suffix to use when saving plots; "
                            "default=%s"%default_plot_type))

    # Add group for default cuts
    group = optparse.OptionGroup(parser,
                                 "LENSING CUTS",
                                 "Default lensing cuts applied to all analyses")

    group.add_option("--use_wl_flag", dest="use_wl_flag",
                     type=int, default=default_use_wl_flag,
                     help=("Use WL flag?; "
                           "default="+str(default_use_wl_flag)))
    group.add_option("--bright_star_mask", dest="bright_star_mask",
                     type=int, default=default_bright_star_mask,
                     help=("Use bright star mask?; "
                           "default="+str(default_bright_star_mask)))
    group.add_option("--reweight_s16a", dest="reweight_s16a",
                     type=int, default=default_reweight_s16a,
                     help=("Reweight sims to S16A conditions?; "
                          "default="+str(default_reweight_s16a)))
    group.add_option("--rescale_est_variance",dest="rescale_est_variance",
                     type=int, default=default_rescale_est_variance,
                     help=("Rescale estimated variance in v4 sims?; "
                           "default="+str(default_rescale_est_variance)))
    group.add_option("--use_calib",dest="use_calib",
                     type=int, default=default_use_calib,
                     help=("Use calibration script for sims, data?; "
                           "default="+str(default_use_calib)))

    parser.add_option_group(group)

    # Add parser options for each analysis
    for a in analyses_list:
        # Create a group per analysis with appropriate name
        name, desc = get_analysis_description(a)
        group = optparse.OptionGroup(parser, name.upper(), desc)

        # Retrieve options for given analysis
        args, vals, desc = get_analysis_arguments(a)

        # Create a prefix for the name of the options
        prefix = name.split(' ')[0] + "_"

        if vals is not None:
            # Add each option to the group
            for (a, v, d) in zip(args, vals, desc):
                group.add_option("--"+prefix+a, dest=prefix+a,
                                 type=type(v),
                                 default=v,
                                 help=d+'; default='+str(v))
        parser.add_option_group(group)

    # Parse command line options
    opts, args = parser.parse_args()

    # Checks that the paths provided are valid
    if not os.path.isdir(opts.sims_dir):
        parser.error("Simulation catalog directory %s does not exist "
                     "or is not a directory" % opts.sims_dir)
    if not os.path.isdir(opts.data_dir):
        parser.error("Data catalog directory %s does not exist "
                     "or is not a directory" % opts.data_dir)
    if not os.path.exists(opts.plot_dir):
        try:
            os.makedirs(opts.plot_dir)
        except:
            parser.error("Plot directory %s does not exist "
                         "and could not be created."% opts.plot_dir)
    if opts.plot_type not in ['pdf','eps','png']:
        raise RuntimeError("Unrecognized plot type: %s"%opts.plot_type)

    # Create complete plot prefix by concatenation of path with prefix
    plot_pref = os.path.join(opts.plot_dir, opts.plot_pref)

    # Initialise each analysis with the arguments
    analyses = []
    for a in analyses_list:
        name, desc = get_analysis_description(a)
        # If only some specific analysis has been requested, only run these ones
        if (len(args) == 0) or (len(args) > 0 and name.split(' ')[0] in args):
            prefix = name.split(' ')[0] + "_"
            print('Initializing analysis: ',name.split(' ')[0])
            inst_args, vals, desc = get_analysis_arguments(a)
            inst_a = a(plot_pref=plot_pref,
                       data_version=opts.data_version,
                       sims_version=opts.sims_version,
                       use_wl_flag=bool(opts.use_wl_flag),
                       bright_star_mask=bool(opts.bright_star_mask),
                       reweight_s16a=bool(opts.reweight_s16a),
                       rescale_est_variance=bool(opts.rescale_est_variance),
                       plot_type=opts.plot_type,
                       **{key: value for key,
                          value in zip(inst_args,
                                       [getattr(opts, prefix + s)
                                        for s in inst_args])})
            inst_a.select_data_sims()
            analyses.append(inst_a)
    use_data = [a.use_data for a in analyses]
    use_sims = [a.use_sims for a in analyses]
    global_use_data = any(use_data)
    global_use_sims = any(use_sims)

    # Cycle through each simulation subfield and feed them to each analysis
    for i in range(opts.n_subfield):
        fieldname = data_fields[i]

        if global_use_sims:
            catalog = utilities.load_sims_catalog_HDF5(dir=opts.sims_dir,
                                                       sims_version=opts.sims_version,
                                                       subfield=data_fields[i],
                                                       rescale_est_variance=opts.rescale_est_variance)
            print("Processing sims fields for %s" % fieldname)

            # use calibrations?
            if bool(opts.use_calib):
                replace_calibs(catalog, opts.data_version, opts.sims_version, opts.data_dir)

            for a in analyses: 
                a.process_sims(catalog, fieldname=fieldname)

        if global_use_data:
            # Load the catalog
            catalog = utilities.load_data_catalog_HDF5(dir=opts.data_dir,
                                                       data_version=opts.data_version,
                                                       subfield=data_fields[i])
            print("Processing data field " + data_fields[i])

            # use calibrations?
            if bool(opts.use_calib):
                replace_calibs(catalog, opts.data_version, opts.sims_version, opts.data_dir)

            for a in analyses:
                a.process_data(catalog, fieldname=fieldname)

    # Finalise and produce plots
    for a in analyses:
        a.finalise()


if __name__ == "__main__":
    main(sys.argv)
