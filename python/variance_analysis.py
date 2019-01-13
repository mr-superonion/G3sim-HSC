import numpy as np
from base_analysis import base_analysis
import pyfits
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['ps.useafm'] = True
matplotlib.rcParams['pdf.use14corefonts'] = True
matplotlib.rcParams['text.usetex'] = True
matplotlib.rcParams.update({'font.size': 16})
import os
import utilities

class variance_analysis(base_analysis):
    """
    Compares noise variance properties in two files.
    """

    def __init__(self,
                 file1='/home/rmandelb.proj/rmandelb/hsc-sims/noise/variance_v2_rand.fits',
                 file2='/home/rmandelb.proj/rmandelb/hsc-sims/noise/variance_s16a.fits',
                 **kwargs):
        """
        Initialise the analysis.

        Parameters
        ----------

        file1: string
            Name of first file to compare

        file2: string
            Name of second file to compare
        """
        super(variance_analysis, self).__init__(**kwargs)

        # Read in the data - only for the 800 subfields used
        self.data1 = pyfits.getdata(file1)[:800]
        # Second file is the one with FDFC cut, all HSC-wide
        # Using the whole thing will make the sims be more like the whole survey
        # without fluctuations due to limited subfield selection.
        self.data2 = pyfits.getdata(file2)
        self.data2 = self.data2[self.data2['fdfc']==1]

        # Range of histograms
        self.hist_nbins = 20
        self.var_max = 0.01

    def select_data_sims(self):
        """
        Decide whether the analysis is going to use sims, data, or both.
        """
        self.use_data = False
        if self.reweight_s16a:
            self.use_sims = True
        else:
            self.use_sims = False

    def process_data(self, data, fieldname):
        """
        Process the data
        """
        pass

    def process_sims(self, data, fieldname):
        """
        Process sims
        """
        if not hasattr(self, 'weight'):
            self.weight = np.zeros(len(self.data1)).astype(float)
        if self.reweight_s16a:
            mask = utilities.galaxy_selector(data, cancel_shape_noise=True,
                                             bright_star_mask=False,
                                             use_wl_flag=self.use_wl_flag)
            # find which subfields are covered by this data
            subfield_id = data['num'] / 10000
            sids = np.unique(subfield_id)
            for sid in sids:
                use_indx = np.where(subfield_id == sid)[0][0]
                self.weight[sid] = data['reweight_s16a'][use_indx]

    def plot_var(self):
        """
        Plot distribution of variance and ratio
        """

        # And do special check for 'v4' only.  The problem with v4 is
        # that not quite all subfields were produced, so we have to
        # allow for their non-existence.  This was automatically taken
        # care of by process_sims for the case of
        # self.reweight_s16a==True, so just do it for the other case.
        if not hasattr(self, 'weight') or not self.reweight_s16a:
            self.weight = np.ones(len(self.data1)).astype(float)
        if not self.reweight_s16a and 'v4' in self.sims_version:
            for subfield in range(len(self.data1)):
                filename = os.path.join(utilities.default_sim_dir, self.sims_version, 'src-%04d.fits' % subfield)
                if not os.path.exists(filename): self.weight[subfield] = 0
        use_subfield = np.where(self.weight > 0)[0]
        print 'Using %d subfields'%len(use_subfield)

        # Extract list of fields
        fields = list(set(self.data1['field_name']) | set(self.data2['field_name']))
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

                # Both files have all fields, so don't worry about checking.
                data1_i = self.data1[self.data1['field_name'] == fields[i]]
                wt1_i = self.weight[self.data1['field_name'] == fields[i]]
                data2_i = self.data2[self.data2['field_name'] == fields[i]]
                min_var_2 = np.min(data2_i['var'])
                h1, ed = np.histogram(data1_i['var'],
                                      bins=self.hist_nbins,
                                      weights=wt1_i,
                                      range=[min_var_2, self.var_max], density=True)
                h2, ed = np.histogram(data2_i['var'], bins=self.hist_nbins,
                                      range=[min_var_2, self.var_max], density=True)
                x = 0.5*(ed[:-1] + ed[1:])
                ax.plot(100*x, h2,  color='b', label='S16A')
                ax.plot(100*x, h1, color='m', ls='--', label='Simulations')
                ax.set_ylim(0,800)
                ax.set_xlim(0,1)
                ax.set_xticks([0,0.2,0.4,0.6,0.8,1])
                ax.set_title(utilities.field_keys[fieldname])
                ax.set_xlabel('100 $\sigma_n^2$')
                self.ratio[fieldname] = (ed,h2.astype(float)/h1)
            elif i == nfields:
                var1 = self.data1['var']
                var2 = self.data2['var']
                h1, ed = np.histogram(var1, bins=self.hist_nbins,
                                      weights=self.weight,
                                      range=[min(var2), self.var_max], density=True)
                h2, ed = np.histogram(var2, bins=self.hist_nbins,
                                      range=[min(var2), self.var_max], density=True)
                x = 0.5*(ed[:-1] + ed[1:])

                ax.plot( 100*x, h2, color='b', label='S16A')
                ax.plot( 100*x, h1, color='m',ls='--', label='Simulations')
                if not self.reweight_s16a: ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
                ax.set_title("All fields")
                ax.set_ylim(0,600)
                ax.set_xlim(0,1)
                ax.set_xticks([0,0.2,0.4,0.6,0.8,1])
                ax.set_xlabel('100 $\sigma_n^2$')
                self.ratio['all'] = (ed,h2.astype(float)/h1)
            else:
                # can rely on var1, var2 from last time
                use_weights = 0*self.weight + 1.0
                h1, ed = np.histogram(var1, bins=self.hist_nbins,
                                      weights=use_weights,
                                      range=[min(var2), self.var_max], density=True)
                ax.plot( 100*x, h2, color='b', label='S16A')
                ax.plot( 100*x, h1, color='m',ls='--', label='Simulations')
                ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
                ax.set_title("All fields, unweighted")
                ax.set_xlabel('100 $\sigma_n^2$')
                ax.set_ylim(0,600)
                ax.set_xlim(0,1)
                ax.set_xticks([0,0.2,0.4,0.6,0.8,1])
        fig.tight_layout()
        fig.subplots_adjust(top=0.9)

        # set infinities / nans to zero before plotting / saving
        for key in self.ratio.keys():
            ed, ratio = self.ratio[key]
            ratio[np.isnan(ratio)] = 0.0
            ratio[np.isinf(ratio)] = 0.0
            self.ratio[key] = (ed,ratio)

        fig2 = plt.figure(figsize=(10,10))
        for i in range(nfields+1):
            ax = fig2.add_subplot(3,3,i+1)
            if i < nfields:
                fieldname = fields[i]

                if fieldname in self.ratio.keys():
                    ed,y = self.ratio[fieldname]
                    x = 0.5*(ed[:-1] + ed[1:])
                    ax.plot(x, y)
                    ax.set_title(fieldname)
            else:
                ed,y = self.ratio['all']
                x = 0.5*(ed[:-1] + ed[1:])
                ax.plot(x, y)
                ax.set_title("All fields")
        fig2.suptitle("Ratio of histograms", fontsize=14)
        fig2.tight_layout()
        fig2.subplots_adjust(top=0.9)

        return fig, fig2

    def finalise(self):
        """
        Make plots
        """
        fig1, fig2 = self.plot_var()

        fig1.savefig(self.plot_pref+'noise_var_comparison.'+self.plot_type)
        fig2.savefig(self.plot_pref+'noise_var_comparison_ratio.'+self.plot_type)
        if not self.reweight_s16a:
            # Only save this for reweighting purposes if the plot was
            # not reweighted already.
            np.save('per_field_ratio_%s.npy'%self.sims_version, self.ratio)
