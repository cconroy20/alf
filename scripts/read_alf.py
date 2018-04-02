import sys
import warnings
from copy import deepcopy
import numpy as np
from numpy.polynomial.chebyshev import chebfit, chebval
from scipy import constants, interpolate
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from astropy.io import ascii
from astropy.table import Table, Column, hstack

class Alf(object):
    def __init__(self, outfiles, read_mcmc=True, info=None):
        self.outfiles = outfiles
        #self.legend = info['label']
        #self.imf_type = info['imf_type']
        self.nsample = None
        self.spectra = None

        if read_mcmc:
            self.mcmc = np.loadtxt('{0}.mcmc'.format(self.outfiles))
        results = ascii.read('{0}.sum'.format(self.outfiles))

        with open('{0}.sum'.format(self.outfiles)) as f:
            for line in f:
                if line[0] == '#':
                    if 'Nwalkers' in line:
                        self.nwalkers = float(line.split('=')[1].strip())
                    elif 'Nchain' in line:
                        self.nchain = float(line.split('=')[1].strip())
                    elif 'Nsample' in line:
                        self.nsample = float(line.split('=')[1].strip())

        self.labels = np.array([
                  'chi2','velz','sigma','logage','zH',
                  'FeH', 'a', 'C', 'N', 'Na', 'Mg', 'Si',
                  'K', 'Ca', 'Ti','V', 'Cr', 'Mn', 'Co',
                  'Ni', 'Cu', 'Sr','Ba', 'Eu', 'Teff',
                  'IMF1', 'IMF2', 'logfy', 'sigma2', 'velz2',
                  'logm7g', 'hotteff', 'loghot', 'fy_logage',
                  'logtrans', 'logemline_h', 'logemline_oii',
                  'logemline_oiii', 'logemline_sii', 'logemline_ni',
                  'logemline_nii', 'jitter', 'IMF3', 'logsky', 'IMF4',
                  'h3', 'h4', 'ML_v','ML_i','ML_k','MW_v', 'MW_i','MW_k'
                  ])

        results = Table(results, names=self.labels)

        """
        0:   Mean of the posterior
        1:   Parameter at chi^2 minimum
        2:   1 sigma error
        3-7: 2.5%, 16%, 50%, 84%, 97.5% CLs
        8-9: lower and upper priors
        """

        types = Column(['mean', 'chi2', 'error',
                        'cl25', 'cl16', 'cl50',
                        'cl84', 'cl98', 'lo_prior',
                        'hi_prior'],
                        name='Type')
        results.add_column(types, index=0)

        """
        Create separate table for abundances
        """
        self.xH = results['Type','a', 'C', 'N', 'Na', 'Mg',
                          'Si', 'K', 'Ca', 'Ti','V', 'Cr',
                          'Mn', 'Co', 'Ni', 'Cu', 'Sr','Ba',
                          'Eu']

        # Creating an empty dict
        # is filled in abundance_correct()
        self.xFe = {}

        self.results = results['Type', 'chi2', 'velz', 'sigma',
                               'logage', 'zH', 'FeH', 'Teff',
                               'IMF1', 'IMF2', 'logfy', 'sigma2',
                               'velz2', 'logm7g', 'hotteff',
                'loghot', 'fy_logage',
                  'logtrans', 'logemline_h', 'logemline_oii',
                  'logemline_oiii', 'logemline_sii', 'logemline_ni',
                  'logemline_nii', 'jitter', 'IMF3', 'logsky', 'IMF4',
                  'h3', 'h4', 'ML_v','ML_i','ML_k','MW_v', 'MW_i','MW_k']

        """
        Read in input data and best fit model

        This isn't going to work correctly if the file
        doesn't exist
        """
        try:
            m = np.loadtxt('{0}.bestspec'.format(self.outfiles))
        except:
            warning = ('Do not have the *.bestspec file')
            warnings.warn(warning)
        data = {}
        data['wave'] = m[:,0]/(1.+self.results['velz'][5]*1e3/constants.c)
        data['m_flux'] = m[:,1] # Model spectrum, normalization applied
        data['d_flux'] = m[:,2] # Data spectrum
        data['snr'] = m[:,3]  # Including jitter and inflated errors
        data['unc'] = 1/m[:,3]
        data['poly'] = m[:,4] # Polynomial used to create m_flux
        data['residual'] = (m[:,1] - m[:,2])/m[:,1] * 1e2
        self.spectra = data

        try:
            m = np.loadtxt('{0}.bestspec2'.format(self.outfiles))
            model = {}
            model['wave'] = m[:,0]
            #model['wave'] = m[:,0]/(1.+self.results['velz'][5]*1e3/constants.c)
            model['flux'] = m[:,1]
            self.ext_model = model
        except:
            self.ext_model = None

        """
        Check the values of the nuisance parameters
        and raise a warning if they are too large.
        """
        #warning = ('\n For {0} {1}={2}, which is '
        #           'larger than acceptable. \n')
        #if self.results['loghot'][0] > -1.0:
        #    warnings.warn(warning.format(self.path, 'loghot',
        #                  self.results['loghot'][0]))

    def get_total_met(self):

        zh = np.where(self.labels == 'zH')
        feh = np.where(self.labels == 'FeH')
        total_met = self.mcmc[:,zh] + self.mcmc[:,feh]

        #Computing errors directly from the chains.
        self.tmet = self.get_cls(total_met)

    def normalize_spectra(self):
        """
        Normalize the data and model spectra
        """
        self.spectra['m_flux_norm'] = deepcopy(self.spectra['m_flux'])
        self.spectra['d_flux_norm'] = deepcopy(self.spectra['d_flux'])
        self.spectra['unc_norm']    = deepcopy(self.spectra['unc'])

        chunks = 1000
        min_ = min(self.spectra['wave'])
        max_ = max(self.spectra['wave'])
        num  = int((max_ - min_)/chunks) + 1

        for i in range(num):
            k = ((self.spectra['wave'] >= min_ + chunks*i) &
                 (self.spectra['wave'] <= min_ + chunks*(i+1)))

            if len(self.spectra['d_flux_norm'][k]) < 10:
                continue

            coeffs = chebfit(self.spectra['wave'][k],
                             self.spectra['d_flux_norm'][k], 2)
            poly = chebval(self.spectra['wave'][k], coeffs)
            self.spectra['d_flux_norm'][k] = self.spectra['d_flux_norm'][k]/poly
            self.spectra['unc_norm'][k] = self.spectra['unc_norm'][k]/poly

            coeffs = chebfit(self.spectra['wave'][k],
                             self.spectra['m_flux_norm'][k], 2)
            poly = chebval(self.spectra['wave'][k], coeffs)
            self.spectra['m_flux_norm'][k] = self.spectra['m_flux_norm'][k]/poly

    def abundance_correct(self, s07=False, b14=False, m11=True):
        """
        Convert abundances from X/H to X/Fe.

        Correct the raw abundance values given
        by ALF.
        """

        # Correction factros from Schiavon 2007, Table 6
        # NOTE: Forcing factors to be 0 for [Fe/H]=0.0,0.2
        lib_feh = [-1.6, -1.4, -1.2, -1.0, -0.8,
                   -0.6, -0.4, -0.2, 0.0, 0.2]
        lib_ofe = [0.6, 0.5, 0.5, 0.4, 0.3, 0.2,
                   0.2, 0.1, 0.0, 0.0]

        if s07:
            #Schiavon 2007
            lib_mgfe = [0.4, 0.4, 0.4, 0.4, 0.29,
                        0.20, 0.13, 0.08, 0.05, 0.04]
            lib_cafe = [0.32, 0.3, 0.28, 0.26, 0.20,
                        0.12, 0.06, 0.02, 0.0, 0.0]
        elif b14:
            # Fitted from Bensby+ 2014
            lib_mgfe = [0.4 , 0.4, 0.4, 0.38, 0.37,
                        0.27, 0.21, 0.12, 0.05, 0.0]
            lib_cafe = [0.32, 0.3, 0.28, 0.26, 0.26,
                        0.17, 0.12, 0.06, 0.0, 0.0]
        elif m11 or (b14 is False and s07 is False):
            # Fitted to Milone+ 2011 HR MILES stars
            lib_mgfe = [0.4, 0.4, 0.4, 0.4, 0.34, 0.22,
                        0.14, 0.11, 0.05, 0.04]
            # from B14
            lib_cafe = [0.32, 0.3, 0.28, 0.26, 0.26,
                        0.17, 0.12, 0.06, 0.0, 0.0]

        # In ALF the oxygen abundance is used
        # a proxy for alpha abundance
        del_alfe = interpolate.interp1d(lib_feh, lib_ofe,
                                        kind='linear',
                                        bounds_error=False,
                                        fill_value='extrapolate')
        del_mgfe = interpolate.interp1d(lib_feh, lib_mgfe,
                                        kind='linear',
                                        bounds_error=False,
                                        fill_value='extrapolate')
        del_cafe = interpolate.interp1d(lib_feh, lib_cafe,
                                        kind='linear',
                                        bounds_error=False,
                                        fill_value='extrapolate')

        zh = np.where(self.labels == 'zH')
        al_corr = del_alfe(self.mcmc[:,zh])
        mg_corr = del_mgfe(self.mcmc[:,zh])
        ca_corr = del_cafe(self.mcmc[:,zh])

        # Assuming Ca~Ti~Si
        group1 = {'Ca', 'Ti', 'Si'}

        # These elements seem to show no net enhancemnt
        # at low metallicity
        group2 = {'C', 'Ca', 'N', 'Cr', 'Ni', 'Na'}

        # These elements we haven't yet quantified
        group3 = {'Ba', 'Eu', 'Sr', 'Cu', 'Co',
                  'K', 'V', 'Mn'}

        for i, col in enumerate(self.xH.colnames):
            feh = np.where(self.labels == 'FeH')
            xh = np.where(self.labels == col)
            xfe = (self.mcmc[:,xh] - self.mcmc[:,feh])
            if col=='Type':
                continue
            elif col=='a':
                xfe_vals = xfe + al_corr
            elif col=='Mg':
                xfe_vals = xfe + mg_corr
            elif col in group1:
                xfe_vals = xfe + ca_corr
            elif col in group2 or col in group3:
                xfe_vals = xfe

            self.xFe[col] = self.get_cls(xfe_vals)

    def plot_model(self, fname):

        chunks = 1000
        min_ = min(self.spectra['wave'])
        max_ = max(self.spectra['wave'])
        num = int((max_ - min_)/chunks) + 1

        with PdfPages(fname) as pdf:
            for i in range(num):
                fig = plt.figure(figsize=(14,9), facecolor='white')
                ax1 = plt.subplot2grid((3,2), (0,0), rowspan=2, colspan=2)
                ax2 = plt.subplot2grid((3,2), (2,0), rowspan=1, colspan=2)

                j = ((self.spectra['wave'] >= min_ + chunks*i) &
                     (self.spectra['wave'] <= min_ + chunks*(i+1)))
                ax1.plot(self.spectra['wave'][j],
                         self.spectra['d_flux_norm'][j],
                         'k-', lw=2, label='Data')

                ax1.plot(self.spectra['wave'][j],
                         self.spectra['m_flux_norm'][j],
                         color='#E32017', lw=2, label='Model')
                ax1.legend(frameon=False)

                ax2.plot(self.spectra['wave'][j], self.spectra['residual'][j],
                            color='#7156A5', lw=2, alpha=0.7)
                ax2.fill_between(self.spectra['wave'][j],
                        -(self.spectra['unc'][j])*1e2,
                        +(self.spectra['unc'][j])*1e2,
                        color='#CCCCCC')
                ax2.set_ylim(-4.9, 4.9)

                ax1.set_ylabel(r'Flux (arbitrary units)',
                               fontsize=22)
                ax2.set_ylabel(r'Residual $\rm \%$',
                               fontsize=22)

                ax2.set_xlabel(r'Wavelength $(\AA)$',
                               fontsize=22, labelpad=10)

                pdf.savefig()

    def plot_corner(self, outname, params, color='k', save=True):
        import corner

        labels = np.array(self.labels)

        use = np.in1d(labels, params)

        figure = corner.corner(self.mcmc[:,use],
                               labels=labels[use],
                               color=color,
                               plot_contours=True)

        plt.tight_layout()
        if save:
            plt.savefig(outname)

    def plot_traces(self, outname):
        plt.cla()
        plt.clf()

        self.nchain = 100
        self.nwalks = 512

        num = len(self.labels)
        data = np.zeros((self.nchain, self.nwalks, num))
        for i in range(0, self.nchain):
            for j in range(0,self.nwalks):
                data[i,j] = self.mcmc[i*510+j]

        full = hstack((self.xH, self.results))
        val = (full['Type_1'] == 'chi2')
        with PdfPages(outname) as pdf:
            for i, (label, trace) in enumerate(zip(self.labels, data.T)):
                fig = plt.figure(figsize=(8,6), facecolor='white')
                #if i == 0: # Don't care to see the chi^2 value
                #    continue
                plt.plot(np.arange(0, self.nchain),
                         data[:,:,i], color='k', alpha=0.1)
                plt.axhline(full[label][val], color='#3399ff')
                plt.xlabel('Step')
                plt.ylabel(label)
                pdf.savefig()
                plt.close()
                plt.cla()

    def plot_posterior(self, fname):
        plt.cla()
        plt.clf()

        fig, axarr = plt.subplots(7, 8, figsize=(40,40),facecolor='white')
        axarr = axarr.reshape(axarr.size,1).copy()
        plt.tick_params(axis='both', which='major', labelsize=15)
        plt.tick_params(axis='both', which='minor', labelsize=10)

        full = hstack((self.xH, self.results))
        val = (full['Type_1'] == 'chi2')
        for i, label in enumerate(self.labels):
            if (label=='ML_k' or label == 'MW_k' or
                np.isnan(full[label][val])==True):
                continue
            axarr[i-1][0].set_ylabel(label, fontsize=16, labelpad=30)

            axarr[i-1][0].hist(self.mcmc[:,i], bins=30,
                              histtype='step', color='k',
                              lw=2, alpha=0.9)
            axarr[i-1][0].axvline(full[label][val], color='#E32017',
                                   alpha=0.85)
            #axarr[i-1][0].autoscale(tight=True)

        plt.tight_layout()
        plt.savefig(fname)

    def get_cls(self, distribution):
        distribution = np.sort(np.squeeze(distribution))

        num = self.nwalkers*self.nchain/self.nsample
        lower = distribution[int(0.160*num)]
        median = distribution[int(0.500*num)]
        upper = distribution[int(0.840*num)]

        return {'cl50': median, 'cl84':  upper, 'cl16': lower}

    def write_params(self):
        pass

if __name__=='__main__':
    pass
