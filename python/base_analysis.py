from abc import ABCMeta, abstractmethod
import re


class base_analysis:
    """
    Abstract class to run an analysis on the simulations
    """
    __metaclass__ = ABCMeta

    def __init__(self, data_version="", sims_version="", plot_pref="plot",
                 use_wl_flag=None, bright_star_mask=None, reweight_s16a=None,
                 rescale_est_variance=None, use_calib=None,
                 plot_type="pdf"):
        """
        Standard initialization of the analysis.

        Parameters
        ----------

        data_version: string
            Version of the data being analyzed

        sims_version: string
            Version of the simulations being analyzed

        plot_pref: string
            Prefix to use when saving plots

        use_wl_flag: int
            use default WL cuts?

        bright_star_mask: int
            apply bright star mask?

        reweight_s16a: int
            reweight sims to S16A observing conditions

        rescale_est_variance: int
            rescale sim outputs to correct for misestimated noise variance

        use_calib: int
            Use calibration script for real data and sims?

        plot_type: string
            suffix to use to specify the output plot type, e.g., "png", "pdf", "eps"
        """
        self.data_version = data_version
        self.sims_version = sims_version
        self.plot_pref = plot_pref
        self.use_wl_flag = use_wl_flag
        self.bright_star_mask = bright_star_mask
        self.reweight_s16a = reweight_s16a
        self.rescale_est_variance = rescale_est_variance
        self.use_calib = use_calib
        self.plot_type = plot_type

    @abstractmethod
    def select_data_sims(self):
        """
        Decide whether the analysis is going to use sims, data, or both.
        All analyses should define this method, and have it set variables
        self.use_data and self.use_sims to True/False.  Then analysis.py 
        will use them decide whether or not to read in data.
        """
        pass

    @abstractmethod
    def process_sims(self, data, fieldname):
        """
        Process a batch of simulated data.

        Parameters
        ----------

        data: Recarray
            Record array from fits file
        """
        pass

    @abstractmethod
    def process_data(self, data, fieldname):
        """
        Process a batch of real data.

        Parameters
        ----------

        data: Recarray
            Record array from fits file

        fieldname: string
            Name of the field
        """
        pass

    @abstractmethod
    def finalise(self):
        """
        Make the plots and exports the results of the analysis
        """
        pass


def get_analysis_arguments(analysis_class):
    """
    Returns list of argument names, default values and description
    """
    # Get list of arguments and default values
    args = analysis_class.__init__.__code__.co_varnames[1:-1] # Expects  the last entry to be additional keyword arguments
    vals = analysis_class.__init__.__defaults__
    doc = analysis_class.__init__.__doc__

    # Extracts the documentation for each argument
    desc = ()
    doc_list = [s.strip() for s in re.split("\n|:", doc)]
    for a in args:
        try:
            ind = doc_list.index(a) + 2
            desc += (doc_list[ind],)
        except:
            desc += ("Argument not documented",)

    return args, vals, desc


def get_analysis_description(analysis_class):
    """
    Returns the name and description of the analysis
    """
    name = analysis_class.__name__.replace('_', ' ')
    description = analysis_class.__doc__.strip()
    return name, description
