"""
Routines to read the output files of the
Absorption Line Fitter (ALF) code.
"""

import sys
import warnings
import numpy as np
from scipy import constants, interpolate
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from astropy.io import ascii

class Alf(object):
    def __init__(self, path, fname, legend):
        self.name = fname
        self.path = '{0}/{1}'.format(path, fname)
        self.legend = legend
        try:
            self.indata = np.loadtxt('{0}.dat'.format(self.path))
        except:
            warning = ('Do not have the input data file')
            warnings.warn(warning)
            self.indata = None
        self.spec   = np.loadtxt('{0}.bestspec'.format(self.path))
        self.mcmc   = None

        self.labels = ['chi2','velz','sigma','logage','zH',
                  'FeH', 'aH', 'CH', 'NH', 'NaH', 'MgH',
                  'SiH', 'KH', 'CaH', 'TiH','VH', 'CrH',
                  'MnH', 'CoH', 'NiH', 'CuH', 'SrH','BaH',
                  'EuH', 'Teff', 'IMF1', 'IMF2', 'logfy',
                  'sigma2', 'velz2', 'logm7g', 'hotteff',
                  'loghot','fy_logage','logtrans', 'logemline_H',
                  'logemline_Oiii','logemline_Sii', 'logemline_Ni',
                  'logemline_Nii','jitter','IMF3', 'logsky', 'IMF4',
                  'ML_r','ML_i','ML_k','MW_r', 'MW_i','MW_k']


        results = ascii.read('{0}.sum'.format(self.path), names=self.labels)

        if len(self.labels) != len(results.colnames):
            error = ('Label array and parameter array '
                     'have different lengths.')
            raise ValueError(error)

        """
        0:   Mean of the posterior
        1:   Parameter at chi^2 minimum
        2:   1 sigma error
        3-7: 2.5%, 16%, 50%, 84%, 97.5% CLs
        8-9: lower and upper priors
        """
        self.params = dict(zip(self.labels, results[0]))
        self.params_chi2 = dict(zip(self.labels, results[1]))
        self.errors = dict(zip(self.labels, results[2]))
        self.cl25 = dict(zip(self.labels, results[3]))
        self.cl16 = dict(zip(self.labels, results[4]))
        self.cl50 = dict(zip(self.labels, results[5]))
        self.cl84 = dict(zip(self.labels, results[6]))
        self.cl98 = dict(zip(self.labels, results[7]))
        self.lo_prior = dict(zip(self.labels, results[8]))
        self.lo_prior = dict(zip(self.labels, results[9]))

        """
        Check the values of the nuisance parameters
        and raise a warning if they are too large.
        """
        warning = ('\n For {0} {1}={2}, which is '
                   'larger than acceptable. \n')
        if self.params['logm7g'] > -1.0:
            warnings.warn(warning.format(self.name, 'logm7g',
                          self.params['logm7g']))
        elif self.params['Teff'] > -1.0:
            warnings.warn(warning.format(self.name, 'Teff',
                          self.params['Teff']))
        elif self.params['loghot'] > -1.0:
            warnings.warn(warning.format(self.name, 'loghot',
                          self.params['loghot']))

        ## Change to read in from *.bestp
        #self.nwalks = 1024
        #self.nchain = 100

    def abundance_correct(self):
        """
        Need to correct the raw abundance values given
        by ALF.

        Use the metallicity-dependent correction factors
        from the literature.

        To-Do:
            Only correcting the mean of the posterior values for now.
            Correct other parameters later.
        """

        # Schiavon 2008, Table 6
        lib_feh = [-1.6, -1.4, -1.2, -1.0, -0.8, -0.6, -0.4, -0.2, 0.0, 0.2]
        lib_ofe = [0.6, 0.5, 0.5, 0.4, 0.3, 0.2, 0.2, 0.1, 0.0, 0.0]
        # Bensby+ 2014
        lib_mgfe = [0.4, 0.4, 0.4, 0.38, 0.37, 0.27, 0.21, 0.12, 0.05, 0.0]
        lib_cafe = [0.32, 0.3, 0.28, 0.26, 0.26, 0.17, 0.12, 0.06, 0.0, 0.0]

        # In ALF the oxygen abundance is used a proxy for alpha abundance
        del_alfe = interpolate.UnivariateSpline(lib_feh, lib_ofe, s=1, k=1)
        del_mgfe = interpolate.UnivariateSpline(lib_feh, lib_mgfe, s=1, k=1)
        del_cafe = interpolate.UnivariateSpline(lib_feh, lib_cafe, s=1, k=1)

        alpha_correction = del_alfe(self.params['zH'] + self.params['FeH'])
        self.params['aH'] = self.params['aH'] + alpha_correction

        mg_correction = del_alfe(self.params['zH'] + self.params['FeH'])
        self.params['MgH'] = self.params['MgH'] + mg_correction

        ca_correction = del_alfe(self.params['zH'] + self.params['FeH'])
        self.params['CaH'] = self.params['CaH'] + ca_correction
        self.params['SiH'] = self.params['SiH'] + ca_correction
        self.params['TiH'] = self.params['TiH'] + ca_correction

    def plot_model(self, wrange, outpath, fname):
        fig = plt.figure(figsize=(14,9), facecolor='white')
        ax1 = plt.subplot2grid((3,2), (0,0), rowspan=2, colspan=2)
        ax2 = plt.subplot2grid((3,2), (2,0), rowspan=1, colspan=2)

        velz = self.params['velz']
        in_wave = self.indata[:,0]/(1.+velz*1e3/constants.c)
        mod_wave = self.spec[:,0]/(1.+velz*1e3/constants.c)

        model = np.interp(in_wave, mod_wave, self.spec[:,1])

        i = (in_wave > wrange[0]) & (in_wave < wrange[1])
        coeffs = np.polynomial.chebyshev.chebfit(in_wave[i],
                    self.indata[:,1][i], 2)
        poly = np.polynomial.chebyshev.chebval(in_wave[i],
                    coeffs)
        ax1.plot(in_wave[i], self.indata[:,1][i]/poly, 'k-', lw=2,
                    label='Data')

        i = (mod_wave > wrange[0]) & (mod_wave < wrange[1])
        coeffs = np.polynomial.chebyshev.chebfit(mod_wave[i],
                    self.spec[:,1][i], 2)
        poly = np.polynomial.chebyshev.chebval(mod_wave[i],
                    coeffs)
        ax1.plot(mod_wave[i], self.spec[:,1][i]/poly, color='#E32017',
                    lw=2, label='Model')
        ax1.legend(frameon=False)

        i = (in_wave > wrange[0]) & (in_wave < wrange[1])
        ax2.plot(in_wave[i], (model[i]-self.indata[:,1][i])/model[i]*1e2,
                    color='#7156A5', lw=2, alpha=0.7)
        ax2.fill_between(in_wave[i], -self.indata[:,2][i]/self.indata[:,1][i]*1e2,
                            self.indata[:,2][i]/self.indata[:,1][i]*1e2,
                            color='#CCCCCC')
        ax2.set_ylim(-4.9, 4.9)

        ax1.set_ylabel(r'Flux (arbitrary units)',fontsize=22)
        ax2.set_ylabel(r'Residual $\rm \%$',fontsize=22)

        ax2.set_xlabel(r'Wavelength $(\AA)$',fontsize=22, labelpad=10)
        ax2.set_xlabel(r'Wavelength $(\AA)$',fontsize=22, labelpad=10)

        plt.savefig('{0}/{1}'.format(outpath, fname))

    def plot_corner(self, outpath):
        import corner

        if self.mcmc is None:
            fname = '{0}.mcmc'.format(self.path)
            self.mcmc = np.loadtxt(fname)

        # This is a dumb way to remove the columns
        # that don't have any dynamic range, which
        # causes an error in corner. I'm sure if I
        # was using Pandas or something there'd be
        # a better way to do this.
        use = []
        for i in range(0,50):
            if not np.all(np.abs(self.mcmc[:,i] - 0.) < 1e-5):
                use.append(i)
        use = np.array(use)
        labels = np.array(self.labels)

        figure = corner.corner(self.mcmc[:,use], labels=labels[use])

        plt.tight_layout()
        plt.savefig('{0}/{1}_corner.pdf'.format(outpath, self.name))


    def plot_traces(self, path):
        outname = ('{0}/{1}_traces.pdf'.format(path, self.name))

        plt.cla()
        plt.clf()

        if self.mcmc is None:
            fname = '{0}.mcmc'.format(self.path)
            self.mcmc = np.loadtxt(fname)

        self.nchain = 100
        self.nwalks = 1020

        num = len(self.params)
        data = np.zeros((self.nchain, self.nwalks, num))
        for i in range(0, self.nchain):
            for j in range(0,self.nwalks):
                data[i,j] = self.mcmc[i*1020+j]

        with PdfPages(outname) as pdf:
            for i, (label, trace) in enumerate(zip(self.labels, data.T)):
                fig = plt.figure(figsize=(8,6), facecolor='white')
                if i == 0: # Don't care to see the chi^2 value
                    continue
                plt.plot(np.arange(0, self.nchain),
                         data[:,:,i], color='k', alpha=0.1)
                plt.axhline(self.params[label], color='#3399ff')
                plt.xlabel('Step')
                plt.ylabel(label)
                pdf.savefig()
                plt.close()
                plt.cla()

    def plot_posterior(self, path):
        plt.cla()
        plt.clf()

        if self.mcmc is None:
            fname = '{0}.mcmc'.format(self.path)
            self.mcmc = np.loadtxt(fname)

        fig, axarr = plt.subplots(7, 7, figsize=(40,40),facecolor='white')
        axarr = axarr.reshape(axarr.size,1).copy()
        plt.tick_params(axis='both', which='major', labelsize=15)
        plt.tick_params(axis='both', which='minor', labelsize=10)

        for i, label in enumerate(self.labels):
            if (label=='chi2' or label=='ML_k' or
                label == 'MW_k' or
                np.isnan(self.params[label])==True):
                continue

            axarr[i-1][0].set_ylabel(label, fontsize=16, labelpad=30)

            axarr[i-1][0].hist(self.mcmc[:,i], bins=30,
                              histtype='step', color='k',
                              lw=2, alpha=0.9)
            axarr[i-1][0].axvline(self.params[label], color='#E32017',
                                   alpha=0.85)
            axarr[i-1][0].autoscale(tight=True)

        plt.tight_layout()
        plt.savefig('{0}/{1}_posteriors.pdf'.format(path, self.name))

    def write_params(self):
        fname = '{0}_parameter_values.txt'.format(self.path)
        with open(fname, 'w') as f:
            for a in self.params.keys():
                f.write('{0:5}: {1:5.5} \n'.format(a, self.params[a]))

if __name__=='__main__':
    pass
