import sys
import warnings
import numpy as np
from scipy import constants, interpolate
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from astropy.io import ascii
from astropy.table import Table, Column, hstack

class Alf(object):
    def __init__(self, path, legend):
        self.path = path
        self.legend = legend
        try:
            self.indata = np.loadtxt('{0}.dat'.format(self.path))
        except:
            warning = ('Do not have the input data file')
            warnings.warn(warning)
            self.indata = None
        try:
            self.spec   = np.loadtxt('{0}.bestspec'.format(self.path))
        except:
            warning = ('Do not have the *.bestspec file')
            warnings.warn(warning)
            self.spec = None
        try:
            self.mcmc = np.loadtxt('{0}.mcmc'.format(self.path))
        except:
            warning = ('Do not have the *.mcmc file')
            warnings.warn(warning)
            self.mcmc = None

        results = ascii.read('{0}.sum'.format(self.path))

        with open('{0}.sum'.format(self.path)) as f:
                for line in f:
                    if line[0] == '#':
                        if 'Nwalkers' in line:
                            self.nwalkers = float(line.split('=')[1].strip())
                        elif 'Nchain' in line:
                            self.nchain = float(line.split('=')[1].strip())
                        elif 'Nsample' in line:
                            self.nsample = float(line.split('=')[1].strip())

        old = False
        if len(results.colnames) == 52:
           self.labels = np.array(['chi2','velz','sigma','logage','zH',
                      'FeH', 'a', 'C', 'N', 'Na', 'Mg',
                      'Si', 'K', 'Ca', 'Ti','V', 'Cr',
                      'Mn', 'Co', 'Ni', 'Cu', 'Sr','Ba',
                      'Eu', 'Teff', 'IMF1', 'IMF2', 'logfy',
                      'sigma2', 'velz2', 'logm7g', 'hotteff',
                      'loghot','fy_logage','logtrans', 'logemline_H',
                      'logemline_Oiii','logemline_Sii', 'logemline_Ni',
                      'logemline_Nii','jitter','IMF3', 'logsky', 'IMF4',
                      'h3', 'h4', 'ML_r','ML_i','ML_k','MW_r', 'MW_i','MW_k'])
        elif len(results.colnames) == 50:
            self.labels = np.array(['chi2','velz','sigma','logage','zH',
                      'FeH', 'a', 'C', 'N', 'Na', 'Mg',
                      'Si', 'K', 'Ca', 'Ti','V', 'Cr',
                      'Mn', 'Co', 'Ni', 'Cu', 'Sr','Ba',
                      'Eu', 'Teff', 'IMF1', 'IMF2', 'logfy',
                      'sigma2', 'velz2', 'logm7g', 'hotteff',
                      'loghot','fy_logage','logtrans', 'logemline_H',
                      'logemline_Oiii','logemline_Sii', 'logemline_Ni',
                      'logemline_Nii','jitter','IMF3', 'logsky', 'IMF4',
                      'ML_r','ML_i','ML_k','MW_r', 'MW_i','MW_k'])
            old = True

        results = Table(results, names=self.labels)
        if old:
            h3 = Column(np.zeros(len(results['chi2'])), name='h3')
            h4 = Column(np.zeros(len(results['chi2'])), name='h4')
            results.add_column(h3, index=43)
            results.add_column(h4, index=44)
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
        Split the big table up
        """

        self.basic = results['Type', 'chi2',
                             'velz', 'sigma',
                             'logage', 'zH',
                             'FeH']

        self.xH = results['Type','a', 'C', 'N', 'Na', 'Mg',
                          'Si', 'K', 'Ca', 'Ti','V', 'Cr',
                          'Mn', 'Co', 'Ni', 'Cu', 'Sr','Ba',
                          'Eu']

        # Creating an empty dict
        # is filled in abundance_correct()
        self.xFe = {}

        self.results = results['Type', 'Teff', 'IMF1',
                               'IMF2', 'logfy', 'sigma2',
                               'velz2', 'logm7g', 'hotteff',
                               'loghot','fy_logage', 'logtrans',
                               'logemline_H', 'logemline_Oiii',
                               'logemline_Sii', 'logemline_Ni',
                               'logemline_Nii','jitter','IMF3',
                               'logsky', 'IMF4', 'h3', 'h4',
                               'ML_r','ML_i', 'ML_k', 'MW_r',
                               'MW_i', 'MW_k']

        """
        Check the values of the nuisance parameters
        and raise a warning if they are too large.
        """
        #warning = ('\n For {0} {1}={2}, which is '
        #           'larger than acceptable. \n')
        #if self.results['loghot'][0] > -1.0:
        #    warnings.warn(warning.format(self.path, 'loghot',
        #                  self.results['loghot'][0]))

        ## Change to read in from *.bestp
        #self.nwalks = 1024
        #self.nchain = 100

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

    def plot_model(self, outpath, info, mock=False):
        val = (self.basic['Type'] == 'mean')
        velz = self.basic['velz'][val]
        in_wave = self.indata[:,0]/(1.+velz*1e3/constants.c)
        mod_wave = self.spec[:,0]/(1.+velz*1e3/constants.c)

        # Find overlapping wavelength range
        min_ = max([in_wave[0], mod_wave[0]])
        max_ = min([in_wave[-1], mod_wave[-1]])

        i = ((in_wave >= min_) & (in_wave <= max_))
        in_wave = in_wave[i]
        in_spec = self.indata[:,1][i]
        in_erro = self.indata[:,2][i]

        i = ((mod_wave >= min_) & (mod_wave <= max_))
        mod_wave = mod_wave[i]
        mod_spec = self.spec[:,1][i]

        model = np.interp(in_wave, mod_wave, mod_spec)

        chunk = 1000
        num = int(max_ - min_)/chunk + 1

        if not mock:
            fstring = (
                       '{0}/{1}_{2}_ssp{3}_fit{4}_imf{5}_'
                       'nad{6}_bh{7}_ns{8}_wd{9}_model_compare.pdf'
                       )
            fname = fstring.format(outpath,
                    self.legend.replace(' ', '_'),
                    info['instrument'], info['ssp_type'],
                    info['fit_type'], info['imf_type'],
                    info['nad'], info['bh_remnants'],
                    info['ns_remnants'], info['wd_remnants'])
        else:
            fname = '{0}/{1}_model_compare.pdf'.format(outpath, info['in_sigma'])
        print fname
        with PdfPages(fname) as pdf:
            for i in range(0, num):
                k = ((mod_wave >= min_+chunk*i) & (mod_wave <= min_+chunk*(i+1)))
                if not np.any(k) or len(mod_wave[k]) <= 10:
                    continue

                j = ((in_wave >= min(mod_wave[k])) & (in_wave <= max(mod_wave[k])))

                fig = plt.figure(figsize=(14,9), facecolor='white')
                ax1 = plt.subplot2grid((3,2), (0,0), rowspan=2, colspan=2)
                ax2 = plt.subplot2grid((3,2), (2,0), rowspan=1, colspan=2)

                coeffs = np.polynomial.chebyshev.chebfit(in_wave[j],
                            in_spec[j], 2)
                poly = np.polynomial.chebyshev.chebval(in_wave[j],
                            coeffs)
                ax1.plot(in_wave[j], in_spec[j]/poly, 'k-', lw=2,
                            label='Data')

                coeffs = np.polynomial.chebyshev.chebfit(mod_wave[k],
                            mod_spec[k], 2)
                poly = np.polynomial.chebyshev.chebval(mod_wave[k],
                            coeffs)
                ax1.plot(mod_wave[k], mod_spec[k]/poly, color='#E32017',
                            lw=2, label='Model')
                ax1.legend(frameon=False)

                ax2.plot(in_wave[j], (model[j]-in_spec[j])/model[j]*1e2,
                            color='#7156A5', lw=2, alpha=0.7)
                ax2.fill_between(in_wave[j],
                        -in_erro[j]/in_spec[j]*1e2,
                        in_erro[j]/in_spec[j]*1e2,
                        color='#CCCCCC')
                ax2.set_ylim(-4.9, 4.9)

                ax1.set_ylabel(r'Flux (arbitrary units)',fontsize=22)
                ax2.set_ylabel(r'Residual $\rm \%$',fontsize=22)

                ax2.set_xlabel(r'Wavelength $(\AA)$',fontsize=22, labelpad=10)

                pdf.savefig()

    def plot_corner(self, outpath):
        import corner

        if self.mcmc is None:
            fname = '{0}.mcmc'.format(self.path)
            self.mcmc = np.loadtxt(fname)

        labels = np.array(self.labels)
        use = np.where((labels=='ML_r') |
                        (labels=='ML_i') |
                        (labels=='IMF1') |
                        (labels=='IMF2')
                        )
        #print self.labels
        #sys.exit()

        figure = corner.corner(self.mcmc[:,use[0]], labels=labels[use[0]])

        plt.tight_layout()
        plt.savefig('{0}/{1}_corner.pdf'.format(outpath, self.legend))

    def plot_traces(self, outpath, info, mock=False):
        if not mock:
            fstring = (
                       '{0}/{1}_{2}_ssp{3}_fit{4}_imf{5}_'
                       'nad{6}_bh{7}_ns{8}_wd{9}_traces.pdf'
                       )
            outname = fstring.format(outpath,
                    self.legend.replace(' ', '_'),
                    info['instrument'], info['ssp_type'],
                    info['fit_type'], info['imf_type'],
                    info['nad'], info['bh_remnants'],
                    info['ns_remnants'], info['wd_remnants'])

        else:
            outname = '{0}/{1}_traces.pdf'.format(outpath, info['in_sigma'])
        plt.cla()
        plt.clf()

        if self.mcmc is None:
            fname = '{0}.mcmc'.format(self.path)
            self.mcmc = np.loadtxt(fname)

        self.nchain = 100
        self.nwalks = 510

        num = len(self.labels)
        data = np.zeros((self.nchain, self.nwalks, num))
        for i in range(0, self.nchain):
            for j in range(0,self.nwalks):
                data[i,j] = self.mcmc[i*510+j]

        full = hstack((self.basic, self.xH, self.results))
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

    def plot_posterior(self, path, info, mock=False):
        plt.cla()
        plt.clf()

        if self.mcmc is None:
            fname = '{0}.mcmc'.format(self.path)
            self.mcmc = np.loadtxt(fname)

        fig, axarr = plt.subplots(7, 8, figsize=(40,40),facecolor='white')
        axarr = axarr.reshape(axarr.size,1).copy()
        plt.tick_params(axis='both', which='major', labelsize=15)
        plt.tick_params(axis='both', which='minor', labelsize=10)

        full = hstack((self.basic, self.xH, self.results))
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
        if not mock:
            fstring = (
                       '{0}/{1}_{2}_ssp{3}_fit{4}_imf{5}_'
                       'nad{6}_bh{7}_ns{8}_wd{9}_posterior.pdf'
                       )
            fname = fstring.format(path,
                    self.legend.replace(' ', '_'),
                    info['instrument'], info['ssp_type'],
                    info['fit_type'], info['imf_type'],
                    info['nad'], info['bh_remnants'],
                    info['ns_remnants'], info['wd_remnants'])
        else:
            fname = '{0}/{1}_posterior.pdf'.format(path, info['in_sigma'])
        plt.savefig(fname)

    def get_cls(self, distribution):
        distribution = np.sort(np.squeeze(distribution))


        num = self.nwalkers*self.nchain/self.nsample
        lower = distribution[int(0.160*num)]
        median = distribution[int(0.500*num)]
        upper = distribution[int(0.840*num)]
        return {'cl50': median, 'cl84':  upper, 'cl16': lower}

    def write_params(self):
        fname = '{0}_parameter_values.txt'.format(self.path)
        with open(fname, 'w') as f:
            for a in self.params.keys():
                f.write('{0:5}: {1:5.5} \n'.format(a, self.params[a]))

if __name__=='__main__':
    pass
