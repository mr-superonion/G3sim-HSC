import numpy as np
from utilities import get_abs_ellip, get_res, get_snr, galaxy_selector, get_psf_ellip, get_psf_size
from base_analysis import base_analysis
import pyfits as fits
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['ps.useafm'] = True
matplotlib.rcParams['pdf.use14corefonts'] = True
matplotlib.rcParams['text.usetex'] = True
matplotlib.rcParams.update({'font.size': 16})
import matplotlib.ticker as ticker

class psf_analysis(base_analysis):
    """
    Compares PSF properties in the simulation to the actual data.
    """

    def __init__(self, **kwargs):
        """
        Initialise the analysis with default parameters.
        Applies permissive cuts to simulations and weight each subfield psf by the
        number of galaxies after the cuts.

        Parameters
        ----------

        """
        super(psf_analysis, self).__init__(**kwargs)

        # Range of histograms
        self.hist_sims_nbins = 30
        self.hist_data_nbins = 30
        self.e_range=[-0.2, 0.2]
        self.s_range=[0.4, 0.85]

        # Arrays to store the results
        self.sims_e1 = {}
        self.sims_e2 = {}
        self.sims_size = {}
        self.sims_weight = {}
        self.data_e1 = {}
        self.data_e2 = {}
        self.data_size = {}
        self.data_weight = {}

    def select_data_sims(self):
        """
        Decide whether the analysis is going to use sims, data, or both.
        """
        self.use_data = True
        self.use_sims = True

    def process_data(self, data, fieldname):
        """
        Process the data
        """
        mask = galaxy_selector(data, cancel_shape_noise=False,
                               bright_star_mask=self.bright_star_mask,
                               use_wl_flag=self.use_wl_flag)

        e1, e2 = get_psf_ellip(data[mask])
        s = get_psf_size(data[mask], type='fwhm-arcsec')

        # Check if we are reweighting to match conditions in data and sims.
        if self.reweight_s16a:
            weight = data[mask]['reweight_s16a']
        else:
            weight = np.ones_like(e1).astype(float)

        self.data_e1[fieldname] = e1
        self.data_e2[fieldname] = e2
        self.data_size[fieldname] = s
        self.data_weight[fieldname] = weight

    def process_sims(self, data, fieldname):
        """
        Process sims
        """
        # Compute the quantities of interest for all galaxies
        mask = galaxy_selector(data, cancel_shape_noise=True,
                               bright_star_mask=False,
                               use_wl_flag=self.use_wl_flag)

        e1, e2 = get_psf_ellip(data[mask])
        s = get_psf_size(data[mask], type='fwhm-arcsec')

        weight = np.ones_like(s)
        if self.reweight_s16a:
            weight *= data['reweight_s16a'][mask]

        self.sims_e1[fieldname] = e1
        self.sims_e2[fieldname] = e2
        self.sims_size[fieldname] = s
        self.sims_weight[fieldname] = weight

    def plot_size(self):
        """
        Plot distribution of the psf sizes
        """
        # Extract which fields have been analysed
        fields = list(set(self.data_e1.keys()) | set(self.sims_e1.keys()))
        fields.sort()
        nfields = len(fields)

        fig = plt.figure(figsize=(10,10))
        self.ratio = {}
        if self.reweight_s16a: nplots = nfields+2
        else: nplots = nfields+1
        for i in range(nplots):
            ax = fig.add_subplot(3,3,i+1)
            if i < nfields:
                fieldname = fields[i]

                if fieldname in self.data_e1:
                    h_data, ed_data = np.histogram(self.data_size[fieldname],
                                                   bins=self.hist_data_nbins, range=self.s_range,density=True,
                                                   weights=self.data_weight[fieldname])
                    size_vals = self.data_size[fieldname]
                    weight_vals = self.data_weight[fieldname]
                    m_data = np.sum(size_vals*weight_vals)/np.sum(weight_vals)
                    x_data = 0.5*(ed_data[:-1] + ed_data[1:])
                    ax.axvline(m_data, color='b')
                    ax.plot( x_data, h_data, 'b')
                    ax.text(m_data+0.02,19., "$%0.3f$"%m_data, color='b', ha='left',va='top')
                    #ax.set_ylim(0, h_data.max())
                if fieldname in self.sims_e1:
                    h_sims, ed_sims = np.histogram(self.sims_size[fieldname],
                                                   bins=self.hist_sims_nbins, range=self.s_range,density=True,
                                                   weights=self.sims_weight[fieldname])
                    size_vals = self.sims_size[fieldname]
                    weight_vals = self.sims_weight[fieldname]
                    m_sims = np.sum(size_vals*weight_vals)/np.sum(weight_vals)
                    x_sims = 0.5*(ed_sims[:-1] + ed_sims[1:])
                    ax.axvline(m_sims, color='m', ls='--')
                    ax.plot( x_sims, h_sims, color='m', ls='--')
                    ax.set_ylim(0, 25)
                    ax.text(m_sims-0.02,19., "$%0.3f$"%m_sims, color='m', ha='right',va='top')
                print '%s: Mean and weight in sims= %f, %f'%\
                    (fieldname, m_sims, np.sum(weight_vals))
                print '%s: Mean and weight in data= %f, %f\n'%\
                    (fieldname, m_data, np.sum(weight_vals))
                self.ratio[fieldname] = (ed_sims, h_data.astype(float)/h_sims)
                ax.set_title(fields[i])
                ax.set_xlabel("FWHM [arcsec]")
            elif i == nfields:
                dat = np.concatenate([self.data_size[d] for d in self.data_size]).flatten()
                dat_weight = np.concatenate(
                    [self.data_weight[d] for d in self.data_weight]).flatten()
                sim = np.concatenate([self.sims_size[d] for d in self.sims_size]).flatten()
                sim_weight = np.concatenate(
                    [self.sims_weight[d] for d in self.sims_weight]).flatten()
                h_data, ed_data = np.histogram(dat, bins=self.hist_data_nbins, range=self.s_range, density=True,
                                               weights=dat_weight)
                h_sims, ed_sims = np.histogram(sim, bins=self.hist_sims_nbins, range=self.s_range, density=True,
                                               weights=sim_weight)
                m_data = np.sum(dat*dat_weight)/np.sum(dat_weight)
                m_sims = np.sum(sim*sim_weight)/np.sum(sim_weight)
                x_data = 0.5*(ed_data[:-1] + ed_data[1:])
                x_sims = 0.5*(ed_sims[:-1] + ed_sims[1:])

                ax.plot( x_sims, h_sims, color='m', ls='--', label='Simulations')
                ax.plot( x_data, h_data, color='b', label='S16A')
                self.ratio['all'] = (ed_sims, h_data.astype(float)/h_sims)
                ax.axvline(m_sims, color='m', ls='--')
                ax.axvline(m_data, color='b')
                ax.text(m_sims-0.02,9.6, "$%0.3f$"%m_sims, color='m', ha='right',va='top')
                ax.text(m_data+0.02,9.6, "$%0.3f$"%m_data, color='b', ha='left',va='top')
                ax.set_xlabel("FWHM [arcsec]")
                #ax.text(0.6, 6.5, 'data mean: %.2f\nsims mean: %.2f'%(m_data, m_sims))
                ax.set_ylim(0, 10)
                # Put a legend here only if we're not making another panel.
                if not self.reweight_s16a: ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
                ax.set_title("All fields")
            else:
                # can rely on dat, sim, etc. from last time (all fields), so just remove the reweight_s16a factors
                dat_weight = 1.+0*dat_weight
                sim_weight = 1.+0*sim_weight
                m_data = np.sum(dat*dat_weight)/np.sum(dat_weight)
                m_sims = np.sum(sim*sim_weight)/np.sum(sim_weight)
                h_data, ed_data = np.histogram(dat, bins=self.hist_data_nbins, range=self.s_range, density=True,
                                               weights=dat_weight)
                h_sims, ed_sims = np.histogram(sim, bins=self.hist_sims_nbins, range=self.s_range, density=True,
                                               weights=sim_weight)

                ax.plot( x_sims, h_sims, color='m', ls='--', label='Simulations')
                ax.plot( x_data, h_data, color='b', label='S16A')
                ax.axvline(m_sims, color='m', ls='--')
                ax.axvline(m_data, color='b')
                ax.text(m_sims-0.02,9.6, "$%0.3f$"%m_sims, color='m', ha='right',va='top')
                ax.text(m_data+0.02,9.6, "$%0.3f$"%m_data, color='b', ha='left',va='top')
                ax.set_xlabel("FWHM [arcsec]")
                ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
                ax.set_ylim(0, 10)
                ax.set_title("All fields, unweighted")

        #fig.suptitle("Distribution of PSF FWHM [arcsec] for sims and data", fontsize=14)
        fig.tight_layout()
        fig.subplots_adjust(top=0.9)
        return fig

    def plot_e(self):
        """
        Plot distribution of psf e1 and e2
        """
        # Extract which fields have been analysed
        fields = list(set(self.data_e1.keys()) | set(self.sims_e1.keys()))
        fields.sort()
        nfields = len(fields)

        fig = plt.figure(figsize=(10,10))
        if self.reweight_s16a: nplots = nfields+2
        else: nplots = nfields+1
        for i in range(nplots):
            ax = fig.add_subplot(3,3,i+1)
            if i < nfields:
                fieldname = fields[i]
                if fieldname in self.data_e1:
                    weight = self.data_weight[fieldname]
                    h1_data, ed_data = np.histogram(self.data_e1[fieldname],
                                                    bins=self.hist_data_nbins, range=self.e_range,density=True,
                                                    weights=weight)
                    h2_data, ed_data = np.histogram(self.data_e2[fieldname],
                                                    bins=self.hist_data_nbins, range=self.e_range,density=True,
                                                    weights=weight)
                    x_data = 0.5*(ed_data[:-1] + ed_data[1:])
                    ax.plot( x_data, h1_data, 'b', label=r'$e_1$ %s'%self.data_version)
                    ax.plot( x_data, h2_data, 'r', label=r'$e_2$ %s'%self.data_version)
                if fieldname in self.sims_e1:
                    weight = self.sims_weight[fieldname]
                    h1_sims, ed_sims = np.histogram(self.sims_e1[fieldname],
                                                    bins=self.hist_sims_nbins, range=self.e_range,density=True,
                                                    weights=weight)
                    h2_sims, ed_sims = np.histogram(self.sims_e2[fieldname],
                                                    bins=self.hist_sims_nbins, range=self.e_range,density=True,
                                                    weights=weight)
                    x_sims = 0.5*(ed_sims[:-1] + ed_sims[1:])
                    ax.plot( x_sims, h1_sims, color='c',ls='--', label=r'$e_1$ sims_%s'%self.sims_version)
                    ax.plot( x_sims, h2_sims, color='m',ls=':', label=r'$e_2$ sims_%s'%self.sims_version)
                ax.set_ylim(0,17)
                ax.set_title(fieldname)
                ax.xaxis.set_major_locator(ticker.MultipleLocator(0.1))
            elif i == nfields:
                sim_e1 = np.concatenate([self.sims_e1[d] for d in self.sims_e1]).flatten()
                sim_e2 = np.concatenate([self.sims_e2[d] for d in self.sims_e2]).flatten()
                sim_weight = np.concatenate([self.sims_weight[d] for d in self.sims_weight]).flatten()
                dat_e1 = np.concatenate([self.data_e1[d] for d in self.data_e1]).flatten()
                dat_e2 = np.concatenate([self.data_e2[d] for d in self.data_e2]).flatten()
                dat_weight = np.concatenate([self.data_weight[d] for d in self.data_weight]).flatten()
                h1_sims, ed_sims = np.histogram(sim_e1, bins=self.hist_sims_nbins,
                                                range=self.e_range,density=True, weights=sim_weight)
                h2_sims, ed_sims = np.histogram(sim_e2, bins=self.hist_sims_nbins,
                                                range=self.e_range,density=True, weights=sim_weight)
                h1_data, ed_data = np.histogram(dat_e1, bins=self.hist_data_nbins,
                                                range=self.e_range,density=True, weights=dat_weight)
                h2_data, ed_data = np.histogram(dat_e2, bins=self.hist_data_nbins,
                                                range=self.e_range,density=True, weights=dat_weight)
                x_sims = 0.5*(ed_sims[:-1] + ed_sims[1:])
                x_data = 0.5*(ed_data[:-1] + ed_data[1:])

                ax.plot( x_data, h1_data, color='b', label=r'$e_1$ S16A')
                ax.plot( x_data, h2_data, color='r', label=r'$e_2$ S16A')
                ax.plot( x_sims, h1_sims, color='c',ls='--', label=r'$e_1$ Simulations')
                ax.plot( x_sims, h2_sims, color='m',ls=':', label=r'$e_2$ Simulations')
                if not self.reweight_s16a: ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
                ax.set_title("All fields")
                ax.set_ylim(0,17)
                ax.xaxis.set_major_locator(ticker.MultipleLocator(0.1))
            else:
                dat_weight = 1.+0*dat_weight
                sim_weight = 1.+0*sim_weight
                h1_sims, ed_sims = np.histogram(sim_e1, bins=self.hist_sims_nbins,
                                                range=self.e_range,density=True, weights=sim_weight)
                h2_sims, ed_sims = np.histogram(sim_e2, bins=self.hist_sims_nbins,
                                                range=self.e_range,density=True, weights=sim_weight)
                h1_data, ed_data = np.histogram(dat_e1, bins=self.hist_data_nbins,
                                                range=self.e_range,density=True, weights=dat_weight)
                h2_data, ed_data = np.histogram(dat_e2, bins=self.hist_data_nbins,
                                                range=self.e_range,density=True, weights=dat_weight)
                ax.plot( x_data, h1_data, color='b', label=r'$e_1$ S16A')
                ax.plot( x_data, h2_data, color='r', label=r'$e_2$ S16A')
                ax.plot( x_sims, h1_sims, color='c',ls='--', label=r'$e_1$ Simulations')
                ax.plot( x_sims, h2_sims, color='m',ls=':', label=r'$e_2$ Simulations')
                ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
                ax.set_title("All fields, unweighted")
                ax.set_ylim(0,17)
                ax.xaxis.set_major_locator(ticker.MultipleLocator(0.1))
        #fig.suptitle("Distribution of PSF ellipticity for sims and data", fontsize=14)
        fig.tight_layout()
        fig.subplots_adjust(top=0.9)
        return fig

    def finalise(self):
        """
        Make plots
        """
        fig1 = self.plot_e()
        fig2 = self.plot_size()

        fig1.savefig(self.plot_pref+'psf_e.'+self.plot_type)
        fig2.savefig(self.plot_pref+'psf_size.'+self.plot_type)
        # save seeing ratio for later use
        np.save('per_field_seeing_ratio_%s.npy'%self.sims_version, self.ratio)
