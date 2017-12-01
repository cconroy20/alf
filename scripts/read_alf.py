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
import fsps

class Alf(object):
    def __init__(self, outfiles, read_mcmc=True, info=None):
        self.outfiles = outfiles
        #self.legend = info['label']
        #self.imf_type = info['imf_type']
        if read_mcmc:
            self.mcmc = np.loadtxt('{0}.mcmc'.format(self.outfiles))
        results = ascii.read('{0}.sum'.format(self.outfiles))

        self.nsample = None
        with open('{0}.sum'.format(self.outfiles)) as f:
            for line in f:
                if line[0] == '#':
                    if 'Nwalkers' in line:
                        self.nwalkers = float(line.split('=')[1].strip())
                    elif 'Nchain' in line:
                        self.nchain = float(line.split('=')[1].strip())
                    elif 'Nsample' in line:
                        self.nsample = float(line.split('=')[1].strip())
        if not self.nsample:
            # The old files don't have this
            # in the header. This is just a
            # guess at the default. Might
            # need to change.
            self.nsample = 1

        old = False
        if len(results.colnames) == 52:
           self.labels = np.array([
                      'chi2','velz','sigma','logage','zH',
                      'FeH', 'a', 'C', 'N', 'Na', 'Mg', 'Si',
                      'K', 'Ca', 'Ti','V', 'Cr', 'Mn', 'Co',
                      'Ni', 'Cu', 'Sr','Ba', 'Eu', 'Teff',
                      'IMF1', 'IMF2', 'logfy', 'sigma2', 'velz2',
                      'logm7g', 'hotteff', 'loghot','fy_logage',
                      'logtrans', 'logemline_H', 'logemline_Oiii',
                      'logemline_Sii', 'logemline_Ni', 'logemline_Nii',
                      'jitter','IMF3', 'logsky', 'IMF4', 'h3', 'h4',
                      'ML_v','ML_i','ML_k','MW_v', 'MW_i','MW_k'
                      ])
        elif len(results.colnames) == 50:
            self.labels = np.array([
                      'chi2','velz','sigma','logage','zH',
                      'FeH', 'a', 'C', 'N', 'Na', 'Mg', 'Si',
                      'K', 'Ca', 'Ti','V', 'Cr', 'Mn', 'Co',
                      'Ni', 'Cu', 'Sr','Ba', 'Eu', 'Teff',
                      'IMF1', 'IMF2', 'logfy', 'sigma2', 'velz2',
                      'logm7g', 'hotteff', 'loghot','fy_logage',
                      'logtrans', 'logemline_H', 'logemline_Oiii',
                      'logemline_Sii', 'logemline_Ni', 'logemline_Nii',
                      'jitter','IMF3', 'logsky', 'IMF4',
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
                               'loghot','fy_logage', 'logtrans',
                               'logemline_H', 'logemline_Oiii',
                               'logemline_Sii', 'logemline_Ni',
                               'logemline_Nii','jitter','IMF3',
                               'logsky', 'IMF4', 'h3', 'h4',
                               'ML_v','ML_i', 'ML_k', 'MW_v',
                               'MW_i', 'MW_k']

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
            self.data = None
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

        #self.mass = None

        """
        Check the values of the nuisance parameters
        and raise a warning if they are too large.
        """
        #warning = ('\n For {0} {1}={2}, which is '
        #           'larger than acceptable. \n')
        #if self.results['loghot'][0] > -1.0:
        #    warnings.warn(warning.format(self.path, 'loghot',
        #                  self.results['loghot'][0]))

    def get_total_met(obj):

        zh = np.where(obj.labels == 'zH')
        feh = np.where(obj.labels == 'FeH')
        total_met = obj.mcmc[:,zh] + obj.mcmc[:,feh]

        #Computing errors directly from the chains.
        return obj.get_cls(total_met)

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
        num  = (int(max_ - min_)/chunks) + 1

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

    def get_m2l(self, info, filters, in_=False, vega=False, remnants=True, mw=0):
        """
        Note: That this does not work correctly. It needs to work directly from
        the chains.
        """


        # Taken from alf_vars.f90
        imflo = 0.08
        imfhi = 100.0

        msto_t0 = +0.33250847
        msto_t1 = -0.29560944
        msto_z0 = +0.95402521
        msto_z1 = +0.21944863
        msto_z2 = +0.070565820

        krpa_imf1 = 1.3
        krpa_imf2 = 2.3
        krpa_imf3 = 2.3

        val = (self.results['Type'] == 'cl50')
        logage = self.results['logage'][val][0]
        zh = self.results['zH'][val][0]

        # line 546 in alf.f90
        msto = max(min(10**(msto_t0+msto_t1*logage)*\
                       (msto_z0+msto_z1*zh+msto_z2*zh**2), 3.0), 0.75)

        if mw == 1:
            mass = get_mass(imflo, msto, krpa_imf1, krpa_imf2,
                    krpa_imf3, remnants=remnants)
        else:
            if self.imf_type == 0:
                if in_ == False:
                    val = np.where(self.labels == 'IMF1')
                    imf1 = self.mcmc[:,val]
                else:
                    imf1 = info['in_imf1']
                mass = get_mass(imflo, msto, imf1, imf1, krpa_imf3,
                        remnants=remnants)

            elif self.imf_type == 1:
                # Double power-law IMF with a fixed low-mass cutoff
                if in_ == False:
                    val = np.where(self.labels == 'IMF1')
                    imf1 = self.mcmc[:,val]
                    val = np.where(self.labels == 'IMF2')
                    imf2 = self.mcmc[:,val]

                else:
                    imf1 = info['in_imf1']
                    imf2 = info['in_imf2']

                mass = get_mass(imflo, msto, imf1, imf2, krpa_imf3,
                        remnants=remnants)

            elif self.imf_type == 2:
                pass
            elif self.imf_type == 3:
                # Double power-law IMF with a variable low-mass cutoff
                if in_ == False:
                    val = np.where(self.labels == 'IMF1')
                    imf1 = self.mcmc[:,val]
                    val = np.where(self.labels == 'IMF2')
                    imf2 = self.mcmc[:,val]
                    val = np.where(self.labels == 'IMF3')
                    imf3 = self.mcmc[:,val]
                else:
                    imf1 = info['in_imf1']
                    imf2 = info['in_imf2']
                    imf3 = info['in_mcut']
                mass = get_mass(imf3, msto, imf1, imf2, krpa_imf3,
                        remnants=remnants)
            elif self.imf_type == 4:
                print "Not implemented yet"

        # Covert units of spectrum
        mypi   = 3.14159265
        lsun   = 3.839e33
        clight = 2.9979E10
        pc2cm  = 3.08568E18

        if self.ext_model is None:
            flux = self.spectra['m_flux']/self.spectra['poly']
            wave = self.spectra['wave']
        else:
            flux = self.ext_model['flux']
            wave = self.ext_model['wave']
        aspec = (flux*(wave**2/
                 (clight*1e8))*(lsun/(1e6*4*mypi*pc2cm**2)))

        #print filters
        band = fsps.find_filter(filters)
        if filters == 'v':
            band = band[21]
        else:
            band = band[0]
        band_info = fsps.get_filter(band) # need to fix if ever want more

        t_wave, trans = band_info.transmission
        interptrans = np.interp(wave, t_wave, trans, left=0, right=0)

        # put in check for when transmission curve is out of bounds
        tot_flux = np.trapz(aspec*interptrans,
                    np.log(wave))/np.trapz(interptrans, np.log(wave))

        # In AB
        mag = -2.5*np.log10(tot_flux) - 48.60

        if in_ == False and mw == 0:
            self.mass = self.get_cls(mass)

        msun = band_info.msun_ab
        lum = 10**(2./5 * (msun - mag))

        return mass/lum

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
        num = (int(max_ - min_)/chunks) + 1

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

    def plot_corner(self, params=None):
        """
        Note: still in progress. I'd like to make it so
        people can pass an argument of the parameters
        they'd like on this plot.
        """

        import corner

        labels = np.array(self.labels)
        params = ['chi2', 'velz', 'sigma',
                  'zH', 'logage', 'IMF1', 'IMF2',
                  'ML_r', 'MW_r']
        use = np.in1d(labels, params)

        figure = corner.corner(self.mcmc[:,use],
                               labels=labels[use])#,
                               #plot_contours=True)

        plt.tight_layout()

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

def get_mass(mlo, mto, imf1, imf2, imfup, remnants=True):#, imf3, imf4, timfnorm):

    # Taken from alf_vars.f90
    imflo = 0.08
    imfhi = 100.0

    # From getmass.f90
    # What are these?
    bhlim=40.0
    nslim=8.5
    m2 = 0.5
    m3 = 1.0

    # For IMF type NOT 4

    # Normalize the weights so that 1 Msun formed at t=0
    imfnorm = ((m2**(-imf1+2)-mlo**(-imf1+2))/(-imf1+2) +
               m2**(-imf1+imf2)*(m3**(-imf2+2)-m2**(-imf2+2))/(-imf2+2) +
               m2**(-imf1+imf2)*(imfhi**(-imfup+2)-m3**(-imfup+2))/(-imfup+2))

    # Stars still alive
    getmass = (m2**(-imf1+2)-mlo**(-imf1+2))/(-imf1+2)
    if mto < m3:
        getmass = getmass + m2**(-imf1+imf2)*(mto**(-imf2+2)-m2**(-imf2+2))/(-imf2+2)
    else:
        getmass = getmass + (m2**(-imf1+imf2)*(m3**(-imf2+2)-m2**(-imf2+2))/(-imf2+2) +
                             m2**(-imf1+imf2)*(mto**(-imfup+2)-m3**(-imfup+2))/(-imfup+2))

    getmass = getmass/imfnorm

    if remnants:
        # BH remnants
        # 40<M<imf_up leave behidn a 0.5*M BH
        getmass = getmass + 0.5*m2**(-imf1+imf2)*(imfhi**(-imfup+2)-bhlim**(-imfup+2))/(-imfup+2)/imfnorm

        # NS remnants
        # 8.5<M<40 leave behind a 1.4 Msun NS
        getmass = getmass +  1.4*m2**(-imf1+imf2)*(bhlim**(-imfup+1)-nslim**(-imfup+1))/(-imfup+1)/imfnorm

    # WD remnants
    # M<8.5 leave beind a 0.077*M+0.48 WD
    if mto < m3:
        getmass = getmass + 0.48*m2**(-imf1+imf2)*(nslim**(-imfup+1)-m3**(-imfup+1))/(-imfup+1)/imfnorm
        getmass = getmass + 0.48*m2**(-imf1+imf2)*(m3**(-imf2+1)-mto**(-imf2+1))/(-imf2+1)/imfnorm
        getmass = getmass + 0.077*m2**(-imf1+imf2)*(nslim**(-imfup+2)-m3**(-imfup+2))/(-imfup+2)/imfnorm
        getmass = getmass + 0.077*m2**(-imf1+imf2)*(m3**(-imf2+2)-mto**(-imf2+2))/(-imf2+2)/imfnorm
    else:
        getmass = getmass + 0.48*m2**(-imf1+imf2)*(nslim**(-imfup+1)-mto**(-imfup+1))/(-imfup+1)/imfnorm
        getmass = getmass + 0.077*m2**(-imf1+imf2)*(nslim**(-imfup+2)-mto**(-imfup+2))/(-imfup+2)/imfnorm

    # What's going on lines 221-223 of getmass.f90?

    return getmass

if __name__=='__main__':
    pass
