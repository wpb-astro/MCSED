""" SED fitting class using emcee for parameter estimation

.. moduleauthor:: Greg Zeimann <gregz@astro.as.utexas.edu>

"""
import logging
import sfh
import dust_abs
import dust_emission
import metallicity
import cosmology
import emcee
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import corner
import time
from scipy.integrate import simps
from scipy.interpolate import interp1d
from astropy.constants import c as clight
import numpy as np
from astropy.table import Table, vstack

plt.ioff() 

import seaborn as sns
sns.set_context("talk") # options include: talk, poster, paper
sns.set_style("ticks")
sns.set_style({"xtick.direction": "in","ytick.direction": "in",
               "xtick.top":True, "ytick.right":True,
               "xtick.major.size":12, "xtick.minor.size":4,
               "ytick.major.size":12, "ytick.minor.size":4,
               })


class Mcsed:
    def __init__(self, filter_matrix, wave, ssp_ages, ssp_met,
                 ssp_starspectra, ssp_nebspectra, emlinewave, ssp_emlineflux, 
                 sfh_class, dust_abs_class, dust_em_class, met_class=None,
                 nfreeparams=None, t_birth=None, 
                 starSSP=None, nebSSP=None, emlinefluxSSP=None, 
                 data_fnu=None, data_fnu_e=None, 
                 data_emline=None, data_emline_e=None, emline_dict=None,
                 use_emline_flux=None, linefluxCSPdict=None,
                 data_absindx=None, data_absindx_e=None, absindx_dict=None,
                 use_absorption_indx=None, absindxCSPdict=None,
                 fluxwv=None, fluxfn=None, medianspec=None, spectrum=None,
                 medianstarspec=None, starspectrum=None,
                 mediannebspec=None, nebspectrum=None, 
                 redshift=None, Dl=None, filter_flag=None, 
                 input_params=None, true_fnu=None, true_spectrum=None,
                 true_starspectrum=None, true_nebspectrum=None, 
                 sigma_m=0.1, nwalkers=40, nsteps=1000, 
                 chi2=None, tauISM_lam=None, tauIGM_lam=None):
        ''' Initialize the Mcsed class.

        Init
        ----
        filter_matrix : numpy array (2 dim)
            The filter_matrix has rows of wavelength and columns for each
            filter (can be much larger than the filters used for fitting)
        wave : numpy array (1 dim)
            wavelength for SSP models and all model spectra
        ssp_ages : numpy array (1 dim)
            ages of the SSP models
        ssp_met : numpy array (1 dim)
            metallicities of the SSP models
            assume a grid of values Z, where Z_solar = 0.019
        ssp_starspectra : numpy array (3 dim)
            single stellar population spectrum (stellar component) for each 
            age in ssp_ages and each metallicity in ssp_met 
        ssp_nebspectra : numpy array (3 dim)
            single stellar population spectrum (nebular component) for each 
            age in ssp_ages and each metallicity in ssp_met
            (must have same shape as ssp_starspectra)
            if None, assume stars+gas are combined in ssp_starspectra 
        emlinewave : numpy array (1 dim)
            Rest-frame wavelengths of requested emission lines (emline_dict)
            Corresponds to ssp_emline
        ssp_emlineflux : numpy array (3 dim)
            Emission line SSP grid spanning emlinewave, age, metallicity
            Only includes requested emission lines (from emline_dict)
            Only used for calculating model emission line strengths
            Spectral units are ergs / s / cm2 at 10 pc
        sfh_class : str
            Converted from str to class in initialization
            This is the input class for sfh.  Each class has a common attribute
            which is "sfh_class.get_nparams()" for organizing the total model_params.
            Also, each class has a key function, sfh_class.evaluate(t), with
            the input of time in units of Gyrs
        dust_abs_class : str 
            Converted from str to class in initialization
            This is the input class for dust absorption.
        dust_em_class : str
            Converted from str to class in initialization
            This is the input class for dust emission.
        met_class : str
            Converted from str to class in initialization
            This is the input class for stellar metallicity
        nfreeparams : int
            number of free model parameters
        t_birth : float
            Age of the birth cloud in Gyr
            set from the value provided in config.py
        starSSP : numpy array (2 dim)
            Grid of stellar SSP spectra at current guess of stellar metallicity
            (set from ssp_starspectra)
        nebSSP : numpy array (2 dim)
            Grid of nebular SSP spectra at current guess of stellar metallicity
            (set from ssp_nebspectra)
        emlinefluxSSP : numpy array (2 dim)
            Grid of emission line fluxes at each age in the SSP grid
            (set from ssp_emlineflux)
        data_fnu : numpy array (1 dim)
            Photometry for data.  Length = (filter_flag == True).sum()
        data_fnu_e : numpy array (1 dim)
            Photometric errors for data
        data_emline : Astropy Table (1 dim)
            Emission line fluxes in units ergs / cm2 / s
        data_emline_e : Astropy Table (1 dim)
            Emission line errors in units ergs / cm2 / s
        emline_dict : dictionary
            Keys are emission line names (str)
            Values are a two-element tuple:
                (rest-frame wavelength in Angstroms (float), weight (float))
            emline_list_dict defined in config.py, containing only the 
            emission lines that were also provided in the input file
            (i.e., only the measurements that will be used to constrain the model)
        use_emline_flux : bool
            If emline_dict contains emission lines, set to True. Else, False
        linefluxCSPdict : dict
            Emission-line fluxes for current SED model
        data_absindx : Astropy Table (1 dim)
            Absorption line indices
        data_absindx_e : Astropy Table (1 dim)
            Absorption line index errors
        absindx_dict : dict
            absorption_index_dict defined in config.py, containing only 
            measurements that were also provided in the input file
            (i.e., only the measurements that will be used to constrain the model)
        use_absorption_indx : bool
            True, if index measurements were included in the input file and
            should be used in the model selection
        absindxCSPdict : dict
            Absorption line index measurements for current SED model
        fluxwv : numpy array (1 dim)
            wavelengths of photometric filters
        fluxfn : numpy array (1 dim)
            flux densities of modeled photometry
        medianspec : numpy array (1 dim)
            best-fit SED model (same length as self.wave)
            set after fitting the model
        spectrum : numpy array (1 dim)
            current SED model (same length as self.wave) 
        medianstarspec : numpy array (1 dim)
            best-fit stellar SED model (same length as self.wave)
            set after fitting the model
        starspectrum : numpy array (1 dim)
            current stellar SED model (same length as self.wave) 
        mediannebspec : numpy array (1 dim)
            best-fit nebular SED model (same length as self.wave)
            set after fitting the model
        nebspectrum : numpy array (1 dim)
            current nebular SED model (same length as self.wave) 
        redshift : float
            Redshift of the source
        Dl : float
            Luminosity distance of the galaxy (in units of 10 pc)
        filter_flag : numpy array (1 dim)
            Length = filter_matrix.shape[1], True for filters matching data
        input_params : list
            input parameters for modeling.  Intended for testing fitting
            procedure.
        true_fnu : numpy array (1 dim)
            True photometry for test mode.  Length = (filter_flag == True).sum()
        true_spectrum : numpy array (1 dim)
            truth model spectrum in test model (realized from input_params)
        true_starspectrum : numpy array (1 dim)
            truth model stellar spectrum in test model (realized from input_params)
        true_nebspectrum : numpy array (1 dim)
            truth model nebular spectrum in test model (realized from input_params)
        sigma_m : float
            Fractional error expected from the models.  This is used in
            the log likelihood calculation.  No model is perfect, and this is
            more or less a fixed parameter to encapsulate that.
        nwalkers : int
            The number of walkers for emcee when fitting a model
        nsteps : int
            The number of steps each walker will make when fitting a model
        chi2 : dict
            keys: 'dof', 'chi2', 'rchi2'
            Track the degrees of freedom (accounting for data and model parameters)
            and the chi2 and reduced chi2 of the current fit
        tauISM_lam : numpy array (1 dim)
            Array of effective optical depths as function of wavelength 
            for MW dust correction
        tauIGM_lam : numpy array (1 dim)
            Array of effective optical depths as function of wavelength 
            for IGM gas correction
        '''
        # Initialize all argument inputs
        self.filter_matrix = filter_matrix
        self.wave = wave
        self.ssp_ages = ssp_ages
        self.ssp_met = ssp_met
        self.ssp_starspectra = ssp_starspectra
        self.ssp_nebspectra = ssp_nebspectra
        self.emlinewave = emlinewave
        self.ssp_emlineflux = ssp_emlineflux
        self.dnu = np.abs(np.hstack([0., np.diff(2.99792e18 / self.wave)]))
        self.sfh_class = getattr(sfh, sfh_class)()
        self.dust_abs_class = getattr(dust_abs, dust_abs_class)()
        self.dust_em_class = getattr(dust_emission, dust_em_class)()
        self.met_class = getattr(metallicity, 'stellar_metallicity')()
        self.param_classes = ['sfh_class', 'dust_abs_class', 'met_class',
                              'dust_em_class']
        self.nfreeparams = nfreeparams
        self.t_birth = t_birth
        self.starSSP = None
        self.nebSSP = None
        self.emlinefluxSSP = None
        self.data_fnu = data_fnu
        self.data_fnu_e = data_fnu_e
        self.data_emline = data_emline
        self.data_emline_e = data_emline_e
        self.emline_dict = emline_dict
        self.use_emline_flux = use_emline_flux
        self.linefluxCSPdict = None
        self.data_absindx = data_absindx
        self.data_absindx_e = data_absindx_e
        self.absindx_dict = absindx_dict
        self.use_absorption_indx = use_absorption_indx
        self.absindxCSPdict = None
        self.fluxwv = fluxwv
        self.fluxfn = fluxfn
        self.medianspec = medianspec
        self.spectrum = None
        self.medianstarspec = medianstarspec
        self.starspectrum = None
        self.mediannebspec = mediannebspec
        self.nebspectrum = None
        self.redshift = redshift
        if self.redshift is not None:
            self.set_new_redshift(self.redshift)
        self.Dl = Dl
        self.filter_flag = filter_flag
        self.input_params = input_params
        self.true_fnu = true_fnu
        self.true_spectrum = true_spectrum
        self.true_starspectrum = true_starspectrum
        self.true_nebspectrum = true_nebspectrum
        self.sigma_m = sigma_m
        self.nwalkers = nwalkers
        self.nsteps = nsteps
        self.chi2 = chi2
        self.tauISM_lam = tauISM_lam
        self.tauIGM_lam = tauIGM_lam

        # Set up logging
        self.setup_logging()

    def set_new_redshift(self, redshift):
        ''' Setting redshift

        Parameters
        ----------
        redshift : float
            Redshift of the source for fitting
        '''
        self.redshift = redshift
        # Need luminosity distance to adjust spectrum to distance of the source
        self.Dl = cosmology.Cosmology().luminosity_distance(self.redshift)
        self.sfh_class.set_agelim(self.redshift)

    def setup_logging(self):
        '''Setup Logging for MCSED

        Builds
        -------
        self.log : class
            self.log.info() is for general print and self.log.error() is
            for raise cases
        '''
        self.log = logging.getLogger('mcsed')
        if not len(self.log.handlers):
            # Set format for logger
            fmt = '[%(levelname)s - %(asctime)s] %(message)s'
            fmt = logging.Formatter(fmt)
            # Set level of logging
            level = logging.INFO
            # Set handler for logging
            handler = logging.StreamHandler()
            handler.setFormatter(fmt)
            handler.setLevel(level)
            # Build log with name, mcsed
            self.log = logging.getLogger('mcsed')
            self.log.setLevel(logging.DEBUG)
            self.log.addHandler(handler)

    def remove_waverange_filters(self, wave1, wave2, restframe=True):
        '''Remove filters in a given wavelength range

        Parameters
        ----------
        wave1 : float
            start wavelength of masked range (in Angstroms)
        wave2 : float
            end wavelength of masked range (in Angstroms)
        restframe : bool
            if True, wave1 and wave2 correspond to rest-frame wavelengths
        '''
        wave1, wave2 = np.sort([wave1, wave2])
        if restframe:
            wave_factor = 1. + self.redshift
        else:
            wave_factor = 1.
        loc1 = np.searchsorted(self.wave, wave1 * wave_factor)
        loc2 = np.searchsorted(self.wave, wave2 * wave_factor)
        # account for the case where indices are the same
        if (loc1 == loc2):
            loc2+=1
        maxima = np.max(self.filter_matrix, axis=0)
        try:
            newflag = np.max(self.filter_matrix[loc1:loc2, :], axis=0) < maxima * 0.1
        except ValueError:
            return
        maximas = np.max(self.filter_matrix[:, self.filter_flag], axis=0)
        newflags = np.max(self.filter_matrix[loc1:loc2, self.filter_flag], axis=0) < maximas * 0.1
        self.filter_flag = self.filter_flag * newflag
        if self.true_fnu is not None:
            self.true_fnu = self.true_fnu[newflags]
        self.data_fnu = self.data_fnu[newflags]
        self.data_fnu_e = self.data_fnu_e[newflags]


    def get_filter_wavelengths(self):
        '''Get central wavelengths of photometric filters 
        '''
        wave_avg = np.dot(self.wave, self.filter_matrix[:, self.filter_flag])
        return wave_avg

    def get_filter_fluxdensities(self, spectrum=None):
        '''Convert a spectrum to photometric fluxes for a given filter set.
        The photometric fluxes will be in the same units as the spectrum.
        The spectrum is in microjanskies(lambda) such that
        the photometric fluxes will be in microjanskies.

        Parameters
        ----------
        spectrum : None or 1d array
            if not None, measure the absorption indices using the input spectrum
            (must have same shape as self.wave)

        Returns
        -------
        f_nu : numpy array (1 dim)
            Photometric flux densities for an input spectrum
        '''
        if type(spectrum)==type(None):
            spectrum = self.spectrum.copy()
        f_nu = np.dot(spectrum, self.filter_matrix[:, self.filter_flag])
        return f_nu


    def measure_absorption_index(self, spectrum=None):
        '''
        measure absorption indices using current spectrum

        Parameters
        ----------
        spectrum : None or 1d array
            if not None, measure the absorption indices using the input spectrum
            (must have same shape as self.wave)

        Returns
        -------
        update self.absindxCSPdict, the dictionary of absorption line indices
        '''
        self.absindxCSPdict = {}
        if self.use_absorption_indx:
            # convert the spectrum from units of specific frequency to specific wavelength
            wave = self.wave.copy()
            factor = clight.to('Angstrom/s').value / wave**2.
            if type(spectrum)==type(None):
                spec = self.spectrum * factor
            else:
                spec = spectrum * factor

            for indx in self.absindx_dict.keys():
                wht, wave_indx, wave_blue, wave_red, unit = self.absindx_dict[indx]

                wave_indx = np.array(wave_indx) * (1. + self.redshift)
                wave_blue = np.array(wave_blue) * (1. + self.redshift)
                wave_red  = np.array(wave_blue) * (1. + self.redshift)

                # select appropriate data ranges for blue/red continuum and index
                sel_index = np.array([False]*len(wave))
                sel_index[np.argmin(abs(wave-wave_indx[0])):np.argmin(abs(wave-wave_indx[1]))] = True
                if abs(np.argmin(abs(wave-wave_indx[0]))-np.argmin(abs(wave-wave_indx[1])))<2:
                    sel_index[np.argmin(abs(wave-wave_indx[0])):np.argmin(abs(wave-wave_indx[0]))+2] = True
                sel_blue = np.array([False]*len(wave))
                sel_blue[np.argmin(abs(wave-wave_blue[0])):np.argmin(abs(wave-wave_blue[1]))] = True
                if abs(np.argmin(abs(wave-wave_blue[0]))-np.argmin(abs(wave-wave_blue[1])))<2:
                    sel_blue[np.argmin(abs(wave-wave_blue[0])):np.argmin(abs(wave-wave_blue[0]))+2] = True
                sel_red = np.array([False]*len(wave))
                sel_red[np.argmin(abs(wave-wave_red[0])):np.argmin(abs(wave-wave_red[1]))] = True
                if abs(np.argmin(abs(wave-wave_red[0]))-np.argmin(abs(wave-wave_red[1])))<2:
                    sel_red[np.argmin(abs(wave-wave_red[0])):np.argmin(abs(wave-wave_red[0]))+2] = True

                # estimate continuum in the index:
                fw_blue  = np.dot(spec[sel_blue][0:-1], np.diff(wave[sel_blue])) 
                fw_blue /= np.diff(wave[sel_blue][[0,-1]])
                fw_red   = np.dot(spec[sel_red][0:-1],  np.diff(wave[sel_red]))  
                fw_red  /= np.diff(wave[sel_red][[0,-1]])
                cont_waves = [np.median(wave_blue), np.median(wave_red)]
                cont_fw    = [fw_blue, fw_red]
                coeff = np.polyfit( cont_waves, cont_fw, 1)
                cont_index = coeff[0] * wave[sel_index] + coeff[1]

                # flux ratio of index and continuum
                spec_index = spec[sel_index] / cont_index

                if unit==0: # return measurement in equivalent width (Angstroms)
                    value = np.dot( 1. - spec_index[0:-1], np.diff(wave[sel_index]) )

                if unit==1: # return measurement in magnitudes
                    integral = np.dot( spec_index[0:-1], np.diff(wave[sel_index]) )
                    value = -2.5 * np.log10( integral / np.diff(wave[sel_index][[0,-1]]) ) 

                if unit==2: # return measurement as a flux density ratio (red / blue)
                    value = fw_red / fw_blue

                self.absindxCSPdict[indx] = float(value)


    def set_class_parameters(self, theta):
        ''' For a given set of model parameters, set the needed class variables
        related to SFH, dust attenuation, ect.

        Input
        -----
        theta : list
            list of input parameters for sfh, dust attenuation, 
            stellar metallicity, and dust emission
        '''
        start_value = 0
        ######################################################################
        # STAR FORMATION HISTORY
        self.sfh_class.set_parameters_from_list(theta, start_value)
        start_value += self.sfh_class.get_nparams()

        ######################################################################
        # DUST ATTENUATION
        self.dust_abs_class.set_parameters_from_list(theta, start_value)
        start_value += self.dust_abs_class.get_nparams()

        ######################################################################
        # STELLAR METALLICITY 
        self.met_class.set_parameters_from_list(theta, start_value)
        start_value += self.met_class.get_nparams()

        ######################################################################
        # DUST EMISSION
        self.dust_em_class.set_parameters_from_list(theta, start_value)
        start_value += self.dust_em_class.get_nparams()


    def get_ssp_spectrum(self):
        '''
        Calculate SSP for an arbitrary metallicity (self.met_class.met) given a
        model grid for a range of metallicities (self.ssp_met)

        if left as a free parameter, stellar metallicity (self.met_class.met)
        spans a range of log(Z / Z_solar)

        the SSP grid of metallicities (self.ssp_met) assumes values of Z
        (as opposed to log solar values)

        Returns
        -------
        starSSP : numpy array (2 dim)
            Grid of stellar SSP spectra at current guess of stellar metallicity
            (set from ssp_starspectra)
        nebSSP : numpy array (2 dim)
            Grid of nebular SSP spectra at current guess of stellar metallicity
            (set from ssp_nebspectra)
        emlinefluxSSP : 2-d array
            Single stellar population line fluxes for each age in self.ages

        '''
        if self.met_class.fix_met:
            if self.starSSP is not None:
                return self.starSSP, self.nebSSP, self.emlinefluxSSP
        Z = np.log10(self.ssp_met)
        Zsolar = 0.019
        z = self.met_class.met + np.log10(Zsolar)
        X = Z - z
        wei = np.exp(-(X)**2 / (2. * 0.15**2))
        wei /= wei.sum()
        self.starSSP = np.dot(self.ssp_starspectra, wei)
        if type(self.ssp_nebspectra)==type(None):
            self.nebSSP = None 
        else:
            self.nebSSP  = np.dot(self.ssp_nebspectra,  wei)
        if self.use_emline_flux:
            self.emlinefluxSSP = np.dot(self.ssp_emlineflux, wei)
        else:
            self.emlinefluxSSP = self.ssp_emlineflux[:,:,0]
        return self.starSSP, self.nebSSP, self.emlinefluxSSP

    def build_csp(self, sfr=None):
        '''Build a composite stellar population model for a given star
        formation history, dust attenuation law, and dust emission law.

        In addition to the returns it also modifies a lineflux dictionary

        Returns
        -------
        csp : numpy array (1 dim)
            Composite stellar population model (micro-Jy) at self.redshift
            (both stellar and nebular components)
        starcsp : numpy array (1 dim)
            Composite stellar population model (micro-Jy) at self.redshift
            (stellar component)
        nebcsp : numpy array (1 dim)
            Composite stellar population model (micro-Jy) at self.redshift
            (nebular component)
        mass : float
            Mass for csp given the SFH input
        mdust_eb : float
            Dust mass if dust emission is being fit AND assume energy balance
        '''
        # Collapse for metallicity
        starSSP, nebSSP, emlinefluxSSP = self.get_ssp_spectrum()

        # Need star formation rate from observation back to formation
        if sfr is None:
            sfr = self.sfh_class.evaluate(self.ssp_ages)
        ageval = 10**self.sfh_class.age # Gyr

        # Treat the birth cloud and diffuse component separately
        age_birth = self.t_birth 

        # Get dust-free CSPs, properly accounting for ages
        # ageval sets limit on ssp_ages that are useable in model calculation
        # age_birth separates birth cloud and diffuse components
        sel = (self.ssp_ages > age_birth) & (self.ssp_ages <= ageval)
        sel_birth = (self.ssp_ages <= age_birth) & (self.ssp_ages <= ageval)
        sel_age = self.ssp_ages <= ageval

        # The weight is the linear time between ages of each SSP
        weight = np.diff(np.hstack([0, self.ssp_ages])) * 1e9 * sfr
        weight_orig = weight.copy()
        weight_birth = weight.copy()
        weight_age = weight.copy()
        weight[~sel] = 0
        weight_birth[~sel_birth] = 0
        weight_age[~sel_age] = 0

        if self.ssp_ages.shape[0] != self.ssp_starspectra.shape[1]:
            print(( self.ssp_ages.shape , self.ssp_starspectra.shape) )

        # Cover the two cases where ssp_ages contains ageval and when not
        # A: index of last acceptable SSP age
        A = np.nonzero(self.ssp_ages <= ageval)[0][-1]
        # indices of SSP ages that are too old
        select_too_old = np.nonzero(self.ssp_ages >= ageval)[0]
        if len(select_too_old):
            # B: index of first SSP that is too old
            B = select_too_old[0]
            # only adjust weight if ageval falls between two SSP age gridpoints
            if A != B:
                lw = ageval - self.ssp_ages[A]
                wei = lw * 1e9 * np.interp(ageval, self.ssp_ages, sfr)
                if ageval > age_birth:
                    weight[B] = wei
                if ageval <= age_birth:
                    weight_birth[B] = wei
                weight_age[B] = wei

        # Cover two cases where ssp_ages contains age_birth and when not
        # A: index of last acceptable SSP age
        A = np.nonzero(self.ssp_ages <= age_birth)[0][-1]
        # indices of SSP ages that are too old
        select_too_old = np.nonzero(self.ssp_ages >= age_birth)[0]
        if (len(select_too_old)>0): 
            # B: index of first SSP that is too old
            B = select_too_old[0]
            if A != B:
                lw = age_birth - self.ssp_ages[A]
                wei = lw * 1e9 * np.interp(age_birth, self.ssp_ages, sfr)
                if ageval > age_birth:
                    weight[B] = weight_age[B] - wei
                if ageval >= age_birth:
                    weight_birth[B] = wei
                else:
                    weight_birth[B] = weight_age[B]

        # Finally, do the matrix multiplication using the weights
        starspec_dustfree = np.dot(self.starSSP, weight)
        starspec_birth_dustfree = np.dot(self.starSSP, weight_birth)
        if type(self.nebSSP) != type(None):
            nebspec_dustfree  = np.dot(self.nebSSP, weight)
            nebspec_birth_dustfree  = np.dot(self.nebSSP, weight_birth)
        emlineflux_dustfree = np.dot(self.emlinefluxSSP, weight)
        emlineflux_birth_dustfree = np.dot(self.emlinefluxSSP, weight_birth)
        mass = np.sum(weight_age)

        # Need to correct spectrum for dust attenuation
        Alam = self.dust_abs_class.evaluate(self.wave)
        starspec_dustobscured = starspec_dustfree * 10**(-0.4 * Alam)
        if type(self.nebSSP) != type(None):
            nebspec_dustobscured  = nebspec_dustfree * 10**(-0.4 * Alam)

        # Correct the corresponding birth cloud spectrum separately
        Alam_birth = Alam / self.dust_abs_class.EBV_old_young
        starspec_birth_dustobscured = starspec_birth_dustfree * 10**(-0.4 * Alam_birth)
        if type(self.nebSSP) != type(None):
            nebspec_birth_dustobscured  = nebspec_birth_dustfree * 10**(-0.4 * Alam_birth)

        # Compute attenuation for emission lines
        Alam_emline = self.dust_abs_class.evaluate(self.emlinewave,new_wave=True)
        Alam_emline_birth = Alam_emline / self.dust_abs_class.EBV_old_young
        emlineflux_dustobscured = emlineflux_dustfree * 10**(-0.4*Alam_emline)
        emlineflux_birth_dustobscured = emlineflux_birth_dustfree * 10**(-0.4*Alam_emline_birth)

        # Combine the young and old components
        starspec_dustfree     += starspec_birth_dustfree
        starspec_dustobscured += starspec_birth_dustobscured
        if type(self.nebSSP) != type(None):
            nebspec_dustfree      += nebspec_birth_dustfree
            nebspec_dustobscured  += nebspec_birth_dustobscured
        emlineflux_dustfree     += emlineflux_birth_dustfree
        emlineflux_dustobscured += emlineflux_birth_dustobscured

        # Combine the stellar and nebular components
        if type(self.nebSSP) != type(None):
            spec_dustfree     = starspec_dustfree + nebspec_dustfree
            spec_dustobscured = starspec_dustobscured + nebspec_dustobscured 
        else:
            spec_dustfree     = starspec_dustfree.copy() 
            spec_dustobscured = starspec_dustobscured.copy() 

        if self.dust_em_class.assume_energy_balance:
            # Bolometric luminosity of dust attenuation (for energy balance)
            L_bol = (np.dot(self.dnu, spec_dustfree) - np.dot(self.dnu, spec_dustobscured)) 
            dust_em = self.dust_em_class.evaluate(self.wave)
            L_dust = np.dot(self.dnu,dust_em)
            mdust_eb = L_bol/L_dust
            spec_dustobscured += mdust_eb * dust_em
            if type(self.nebSSP) != type(None):
                nebspec_dustobscured += mdust_eb * dust_em
        else:
            spec_dustobscured += self.dust_em_class.evaluate(self.wave)
            if type(self.nebSSP) != type(None):
                nebspec_dustobscured += self.dust_em_class.evaluate(self.wave)

        # Redshift the spectrum to the observed frame
        csp = np.interp(self.wave, self.wave * (1. + self.redshift),
                        spec_dustobscured * (1. + self.redshift))
        if type(self.nebSSP) != type(None):
            starcsp = np.interp(self.wave, self.wave * (1. + self.redshift),
                                starspec_dustobscured * (1. + self.redshift))
            nebcsp  = np.interp(self.wave, self.wave * (1. + self.redshift),
                                nebspec_dustobscured * (1. + self.redshift))
        else:
            starcsp = np.zeros(csp.shape)
            nebcsp  = np.zeros(csp.shape)

        # Correct for ISM and/or IGM (or neither)
        if self.tauIGM_lam is not None:
            csp     *= np.exp(-self.tauIGM_lam)
            starcsp *= np.exp(-self.tauIGM_lam)
            nebcsp  *= np.exp(-self.tauIGM_lam)
        if self.tauISM_lam is not None:
            csp     *= np.exp(-self.tauISM_lam)
            starcsp *= np.exp(-self.tauISM_lam)
            nebcsp  *= np.exp(-self.tauISM_lam)

        # Correct spectra from 10pc to redshift of the source
        csp     /= self.Dl**2
        starcsp /= self.Dl**2
        nebcsp  /= self.Dl**2

        # Update dictionary of modeled emission line fluxes
        linefluxCSPdict = {}
        if self.use_emline_flux:
            for emline in self.emline_dict.keys():
                indx = np.argmin(np.abs(self.emlinewave 
                                        - self.emline_dict[emline][0]))
                # flux is given in ergs / s / cm2 at 10 pc
                flux = emlineflux_dustobscured[indx]
                # Correct flux from 10pc to redshift of source
                linefluxCSPdict[emline] = flux / self.Dl**2
        self.linefluxCSPdict = linefluxCSPdict

        # Update dictionary of modeled absorption line indices
        self.measure_absorption_index(spectrum=csp)

        if self.dust_em_class.assume_energy_balance:
            return csp, starcsp, nebcsp, mass, mdust_eb
        else:
            return csp, starcsp, nebcsp, mass

    def measure_chi2(self, spectrum):
        '''
        Measure chi2 from the input spectrum. Used in measuring chi2 from 
        the median spectrum, emline fluxes, and absorption line indices

        Parameters
        ----------
        spectrum : 1d array
            same shape as self.wave

        Returns
        -------
        update the chi2 dictionary 
        '''
        # likelihood contribution from the photometry
        model_y = self.get_filter_fluxdensities(spectrum=spectrum)
        inv_sigma2 = 1.0 / (self.data_fnu_e**2 + (model_y * self.sigma_m)**2)
        chi2_term = np.sum((self.data_fnu - model_y)**2 * inv_sigma2)

        # calculate the degrees of freedom and store the current chi2 value
        if not self.chi2:
            dof_wht = list(np.ones(len(self.data_fnu)))

        # likelihood contribution from the absorption line indices
        if self.use_absorption_indx:
            self.measure_absorption_index(spectrum=spectrum)
            for indx in self.absindx_dict.keys():
                unit = self.absindx_dict[indx][-1]
                # if null value, ignore it (null = -99)
                if (self.data_absindx['%s_INDX' % indx]+99 > 1e-10):
                    indx_weight = self.absindx_dict[indx][0]
                    if indx_weight < 1e-10:
                        continue
                    model_indx = self.absindxCSPdict[indx]
                    if unit == 1: # magnitudes
                        model_err = 2.5*np.log10(1.+self.sigma_m)
                    else:
                        model_err = model_indx * self.sigma_m
                    obs_indx = self.data_absindx['%s_INDX' % indx]
                    obs_indx_e = self.data_absindx_e['%s_Err' % indx]
                    sigma2 = obs_indx_e**2. + model_err**2.
                    chi2_term += ( (model_indx - obs_indx)**2 /
                                  sigma2) * indx_weight
                    if not self.chi2:
                        dof_wht.append(indx_weight)

        if self.use_emline_flux:
            # if all lines have null line strengths, ignore 
            if not min(self.data_emline) == max(self.data_emline) == -99:
                for emline in self.emline_dict.keys():
                    if self.data_emline['%s_FLUX' % emline] > -99: # null value
                        emline_wave, emline_weight = self.emline_dict[emline]
                        if emline_weight < 1e-10:
                            continue
                        model_lineflux = self.linefluxCSPdict[emline]
                        model_err = model_lineflux * self.sigma_m
                        lineflux  = self.data_emline['%s_FLUX' % emline]
                        elineflux = self.data_emline_e['%s_ERR' % emline]
                        sigma2 = elineflux**2. + model_err**2.
                        chi2_term += ( (model_lineflux - lineflux)**2 /
                                      sigma2) * emline_weight
                        if not self.chi2:
                            dof_wht.append(emline_weight)

        # record current chi2 and degrees of freedom
        if not self.chi2:
            self.chi2 = {}
            dof_wht = np.array(dof_wht)
            npt = ( sum(dof_wht)**2. - sum(dof_wht**2.) ) / sum(dof_wht) + 1
            self.chi2['dof'] = npt - self.nfreeparams
        self.chi2['chi2']  = chi2_term
        self.chi2['rchi2'] = self.chi2['chi2'] / (self.chi2['dof'] - 1.)

    def lnprior(self):
        ''' Simple, uniform prior for input variables

        Returns
        -------
        0.0 if all parameters are in bounds, -np.inf if any are out of bounds
        '''
        flag = True
        for par_cl in self.param_classes:
            flag *= getattr(self, par_cl).prior()
        if not flag:
            return -np.inf
        else:
            return 0.0

    def lnlike(self):
        ''' Calculate the log likelihood and return the value and stellar mass
        of the model as well as other derived parameters

        Returns
        -------
        log likelihood, mass, sfr10, sfr100, fpdr, mdust_eb : (all float)
            The log likelihood includes a chi2_term and a parameters term.
            The mass comes from building of the composite stellar population
            The parameters sfr10, sfr100, fpdr, mdust_eb are derived in get_derived_params(self)
        '''
        if self.dust_em_class.assume_energy_balance:
            self.spectrum, self.starspectrum, self.nebspectrum, mass, mdust_eb = self.build_csp()
        else:
            self.spectrum, self.starspectrum, self.nebspectrum, mass = self.build_csp()
            mdust_eb = None

        sfr10,sfr100,fpdr = self.get_derived_params()

        # likelihood contribution from the photometry
        model_y = self.get_filter_fluxdensities()
        inv_sigma2 = 1.0 / (self.data_fnu_e**2 + (model_y * self.sigma_m)**2)
        chi2_term = -0.5 * np.sum((self.data_fnu - model_y)**2 * inv_sigma2)
        parm_term = -0.5 * np.sum(np.log(1 / inv_sigma2))

        # calculate the degrees of freedom and store the current chi2 value
        if not self.chi2:
            dof_wht = list(np.ones(len(self.data_fnu)))

        # likelihood contribution from the absorption line indices
        if self.use_absorption_indx:
            for indx in self.absindx_dict.keys():
                unit = self.absindx_dict[indx][-1]
                # if null value, ignore it (null = -99)
                if (self.data_absindx['%s_INDX' % indx]+99 > 1e-10):
                    indx_weight = self.absindx_dict[indx][0]
                    if indx_weight < 1e-10:
                        continue
                    model_indx = self.absindxCSPdict[indx]
                    if unit == 1: # magnitudes
                        model_err = 2.5*np.log10(1.+self.sigma_m)
                    else:
                        model_err = model_indx * self.sigma_m
                    obs_indx = self.data_absindx['%s_INDX' % indx]
                    obs_indx_e = self.data_absindx_e['%s_Err' % indx]
                    sigma2 = obs_indx_e**2. + model_err**2.
                    chi2_term += (-0.5 * (model_indx - obs_indx)**2 /
                                  sigma2) * indx_weight
                    parm_term += -0.5 * np.log(indx_weight * sigma2)
                    if not self.chi2:
                        dof_wht.append(indx_weight) 

        # likelihood contribution from the emission lines
        if self.use_emline_flux:
            # if all lines have null line strengths, ignore 
            if not min(self.data_emline) == max(self.data_emline) == -99:
                for emline in self.emline_dict.keys():
                    if self.data_emline['%s_FLUX' % emline] > -99: # null value
                        emline_wave, emline_weight = self.emline_dict[emline]
                        if emline_weight < 1e-10:
                            continue
                        model_lineflux = self.linefluxCSPdict[emline]
                        model_err = model_lineflux * self.sigma_m
                        lineflux  = self.data_emline['%s_FLUX' % emline]
                        elineflux = self.data_emline_e['%s_ERR' % emline]
                        sigma2 = elineflux**2. + model_err**2.
                        chi2_term += (-0.5 * (model_lineflux - lineflux)**2 /
                                      sigma2) * emline_weight
                        parm_term += -0.5 * np.log(emline_weight * sigma2)
                        if not self.chi2:
                            dof_wht.append(emline_weight)

        # record current chi2 and degrees of freedom
        if not self.chi2:
            self.chi2 = {}
            dof_wht = np.array(dof_wht)
            npt = ( sum(dof_wht)**2. - sum(dof_wht**2.) ) / sum(dof_wht) + 1
            self.chi2['dof'] = npt - self.nfreeparams
        self.chi2['chi2']  = -2. * chi2_term
        self.chi2['rchi2'] = self.chi2['chi2'] / (self.chi2['dof'] - 1.)

        return (chi2_term + parm_term, mass,sfr10,sfr100,fpdr,mdust_eb)

    def lnprob(self, theta):
        ''' Calculate the log probabilty and return the value and stellar mass 
        (as well as derived parameters) of the model

        Returns
        -------
        log prior + log likelihood, [mass,sfr10,sfr100,fpdr,mdust_eb]: (all floats) 
            The log probability is just the sum of the logs of the prior and
            likelihood.  The mass comes from the building of the composite
            stellar population. The other derived parameters are calculated in get_derived_params()
        '''
        self.set_class_parameters(theta)
        lp = self.lnprior()
        if np.isfinite(lp):
            lnl,mass,sfr10,sfr100,fpdr,mdust_eb = self.lnlike()
            if not self.dust_em_class.fixed:
                if self.dust_em_class.assume_energy_balance:
                    return lp + lnl, np.array([mass,sfr10,sfr100,fpdr,mdust_eb])
                else:
                    return lp + lnl, np.array([mass, sfr10, sfr100, fpdr])
            else:
                return lp + lnl, np.array([mass, sfr10, sfr100])
        else:
            if not self.dust_em_class.fixed:
                if self.dust_em_class.assume_energy_balance:
                    return -np.inf, np.array([-np.inf, -np.inf, -np.inf, -np.inf, -np.inf])
                else:
                    return -np.inf, np.array([-np.inf, -np.inf, -np.inf, -np.inf])
            else:
                return -np.inf, np.array([-np.inf, -np.inf, -np.inf])

    def get_init_walker_values(self, kind='ball', num=None):
        ''' Before running emcee, this function generates starting points
        for each walker in the MCMC process.

        Returns
        -------
        pos : np.array (2 dim)
            Two dimensional array with Nwalker x Ndim values
        '''
        # We need an initial guess for emcee so we take it from the model class
        # parameter values and deltas
        init_params, init_deltas, init_lims = [], [], []
        for par_cl in self.param_classes:
            init_params.append(getattr(self, par_cl).get_params())
            init_deltas.append(getattr(self, par_cl).get_param_deltas())
            if len(getattr(self, par_cl).get_param_lims()):
                init_lims.append(getattr(self, par_cl).get_param_lims())
        theta = list(np.hstack(init_params))
        thetae = list(np.hstack(init_deltas))
        theta_lims = np.vstack(init_lims)

        if num is None:
            num = self.nwalkers
        if kind == 'ball':
            pos = emcee.utils.sample_ball(theta, thetae, size=num)
        else:
### WPBWPB delete
#            pos = (np.random.rand(num)[:, np.newaxis] *
#                   (theta_lims[:, 1]-theta_lims[:, 0]) + theta_lims[:, 0])
            ran = (theta_lims[:, 1]-theta_lims[:, 0])[np.newaxis, :]
            pos = (np.random.rand(num, len(theta_lims))*
                   ran*0.8 + theta_lims[np.newaxis, :, 0]+0.1*ran)

        return pos


    def get_param_names(self):
        ''' Grab the names of the parameters for plotting

        Returns
        -------
        names : list
            list of all parameter names
        '''
        names = []
        for par_cl in self.param_classes:
            names.append(getattr(self, par_cl).get_names())
        names = list(np.hstack(names))
        return names

    def get_params(self):
        ''' Grab the the parameters in each class

        Returns
        -------
        vals : list
            list of all parameter values
        '''
        vals = []
        for par_cl in self.param_classes:
            vals.append(getattr(self, par_cl).get_params())
        vals = list(np.hstack(vals))
        self.nfreeparams = len(vals)
        return vals

    def get_param_lims(self):
        ''' Grab the limits of the parameters for making mock galaxies

        Returns
        -------
        limits : numpy array (2 dim)
            an array with parameters for rows and limits for columns
        '''
        limits = []
        for par_cl in self.param_classes:
            limits.append(getattr(self, par_cl).get_param_lims())
        limits = np.array(sum(limits, []))
        return limits

    def fit_model(self):
        ''' Using emcee to find parameter estimations for given set of
        data measurements and errors
        '''
        # Need to verify data parameters have been set since this is not
        # a necessity on initiation
        self.log.info('Fitting model using emcee')
        check_vars = ['data_fnu', 'data_fnu_e', 'redshift', 'filter_flag']
        for var in check_vars:
            if getattr(self, var) is None:
                self.error('The variable %s must be set first' % var)

### WPBWPB do I want ball or not?
        pos = self.get_init_walker_values(kind='not_ball')
        ndim = pos.shape[1]
        start = time.time()
        sampler = emcee.EnsembleSampler(self.nwalkers, ndim, self.lnprob,
                                        a=2.0)
        # Do real run
        sampler.run_mcmc(pos, self.nsteps, rstate0=np.random.get_state())
        end = time.time()
        elapsed = end - start
        self.log.info("Total time taken: %0.2f s" % elapsed)
        self.log.info("Time taken per step per walker: %0.2f ms" %
                      (elapsed / (self.nsteps) * 1000. /
                       self.nwalkers))
        # Calculate how long the run should last
        tau = np.max(sampler.acor)
        burnin_step = int(tau*3)
        self.log.info("Mean acceptance fraction: %0.2f" %
                      (np.mean(sampler.acceptance_fraction)))
        self.log.info("AutoCorrelation Steps: %i, Number of Burn-in Steps: %i"
                      % (np.round(tau), burnin_step))

        if self.dust_em_class.fixed: 
            numderpar = 3
        else: 
            if self.dust_em_class.assume_energy_balance:
                numderpar = 5
            else:
                numderpar = 4
        new_chain = np.zeros((self.nwalkers, self.nsteps, ndim+numderpar+1))
        new_chain[:, :, :-(numderpar+1)] = sampler.chain
        self.chain = sampler.chain
        for i in np.arange(len(sampler.blobs)):
            for j in np.arange(len(sampler.blobs[0])):
                for k in np.arange(len(sampler.blobs[0][0])):
                    x = sampler.blobs[i][j][k]
                    # stellar mass and dust mass
                    if k==0 or k==4: 
                        new_chain[j, i, -(numderpar+1)+k] = np.where((np.isfinite(x)) * (x > 10.),
                                               np.log10(x), -99.)
                    # other derived parameters 
                    else: 
                        new_chain[j, i, -(numderpar+1)+k] = np.where((np.isfinite(x)),np.log10(x), -99.) 
        new_chain[:, :, -1] = sampler.lnprobability
        self.samples = new_chain[:, burnin_step:, :].reshape((-1, ndim+numderpar+1))


    def get_derived_params(self):
        ''' These are not free parameters in the model, but are instead
        calculated from free parameters
        '''
        # Lookback times for past 10 and 100 Myr (avoid t=0 for log purposes)
        t_sfr100 = np.linspace(1.0e-9,0.1,num=251)
        t_sfr10 = np.linspace(1.0e-9,0.01,num=251)
        # Time-averaged SFR over the past 10 and 100 Myr 
        sfrarray = self.sfh_class.evaluate(t_sfr100)
        sfr100 = simps(sfrarray,x=t_sfr100)/(t_sfr100[-1]-t_sfr100[0])
        sfrarray = self.sfh_class.evaluate(t_sfr10)
        sfr10 = simps(sfrarray,x=t_sfr10)/(t_sfr10[-1]-t_sfr10[0])

        if self.dust_em_class.fixed:
            fpdr = None
        else:
            if self.dust_em_class.assume_energy_balance:
                umin,gamma,qpah = self.dust_em_class.get_params()
            else:
                umin,gamma,qpah,mdust = self.dust_em_class.get_params()
            umax = 1.0e6
            fpdr = gamma*np.log(umax/100.) / ((1.-gamma)*(1.-umin/umax) + gamma*np.log(umax/umin))

        return sfr10,sfr100,fpdr


    def set_median_fit(self,rndsamples=200,lnprobcut=7.5):
        '''
        set attributes
        median spectrum and filter flux densities for rndsamples random samples

        Input
        -----
        rndsamples : int
            number of random samples over which to compute medians
        lnprobcut : float
            Some of the emcee chains include outliers.  This value serves as
            a cut in log probability space with respect to the maximum
            probability.  For reference, a Gaussian 1-sigma is 2.5 in log prob
            space.

        Returns
        -------
        self.fluxwv : list (1d)
            wavelengths of filters
        self.fluxfn : list (1d)
            median flux densities of filters
        self.medianspec : list (1d)
            median spectrum (stellar and nebular)
        self.medianstarspec : list (1d)
            median stellar spectrum
        self.mediannebspec : list (1d)
            median nebular spectrum
        self.absindxCSPdict : dict
            update with median absorption line index measurements
        self.linefluxCSPdict : dict
            update with median emission line flux measurements
        self.chi2 : dict
            update dictionary (chi2, reduced chi2, and degrees of freedom)
            as measured from medianspec and median line fluxes / indexes
        '''
### WPBWPB delete
        Table(self.samples).write('samples.dat',format='ascii',overwrite=True)
        chi2sel = (self.samples[:, -1] >
                   (np.max(self.samples[:, -1], axis=0) - lnprobcut))
        nsamples = self.samples[chi2sel, :]

### WPBWPB delete
        print(len(nsamples))
        wv = self.get_filter_wavelengths()
        sp, starsp, nebsp, fn = ([], [], [], [])
        temline, tabsindx, tchi2 = ( Table(), Table(), Table() )
        for i in np.arange(rndsamples):
            ind = np.random.randint(0, nsamples.shape[0])
            self.set_class_parameters(nsamples[ind, :])
            if self.dust_em_class.assume_energy_balance:
                self.spectrum, self.starspectrum, self.nebspectrum, mass, mdust_eb = self.build_csp()
            else:
                self.spectrum, self.starspectrum, self.nebspectrum, mass = self.build_csp()
            fnu = self.get_filter_fluxdensities()
            sp.append(self.spectrum * 1.)
            starsp.append(self.starspectrum * 1.)
            nebsp.append(self.nebspectrum * 1.)
            fn.append(fnu * 1.)
            if self.use_emline_flux:
                tloc = Table()
                if not len(temline):
                    cols = list(self.linefluxCSPdict.keys())
                else:
                    cols = temline.colnames
                for emline in cols:
                    tloc[emline] = [self.linefluxCSPdict[emline]]
                temline = vstack([temline, tloc])
            if self.use_absorption_indx:
                tloc = Table()
                if not len(tabsindx):
                    cols = list(self.absindxCSPdict.keys())
                else:
                    cols = tabsindx.colnames
                for indx in cols:
                    tloc[indx] = [self.absindxCSPdict[indx]]
                tabsindx = vstack([tabsindx, tloc])
        self.medianspec = np.median(np.array(sp), axis=0)
        self.medianstarspec = np.median(np.array(starsp), axis=0)
        self.mediannebspec = np.median(np.array(nebsp), axis=0)
        self.fluxwv = wv
        self.fluxfn = np.median(np.array(fn), axis=0)
        if self.use_emline_flux:
            for emline in temline.colnames:
                self.linefluxCSPdict[emline] = np.median(temline[emline])
        if self.use_absorption_indx:
            for indx in tabsindx.colnames:
                self.absindxCSPdict[indx] = np.median(tabsindx[indx])
        self.measure_chi2(spectrum=self.medianspec)


    def spectrum_plot(self, ax, color=[0.996, 0.702, 0.031], alpha=0.1):
        ''' Make spectum plot for current model '''
        if self.dust_em_class.assume_energy_balance:
            self.spectrum, self.starspectrum, self.nebspectrum, mass, mdust_eb = self.build_csp()
        else:
            self.spectrum, self.starspectrum, self.nebspectrum, mass = self.build_csp()
        self.true_spectrum = self.spectrum.copy()
        self.true_starspectrum = self.starspectrum.copy()
        self.true_nebspectrum = self.nebspectrum.copy()
        ax.plot(self.wave, self.spectrum, color=color, alpha=alpha)

    def add_sfr_plot(self, ax1):
        ax1.set_xscale('log')
        ax1.set_yscale('log')
        ax1.set_ylabel(r'SFR [$M_{\odot} yr^{-1}$]')
        ax1.set_xlabel('Lookback Time') 
        ax1.set_xticks([1e-3, 1e-2, 1e-1, 1])
        ax1.set_xticklabels(['1 Myr', '10 Myr', '100 Myr', '1 Gyr'])
        ax1.set_yticks([1e-2, 1e-1, 1, 1e1, 1e2, 1e3])
        ax1.set_yticklabels(['0.01', '0.1', '1', '10', '100', '1000'])
        ax1.set_xlim([10**-3, max(10**self.sfh_class.age, 1.)])
        ax1.set_ylim([10**-2.3, 1e3])
        ax1.minorticks_on()

    def add_dust_plot(self, ax2):
        ax2.set_xscale('log')
        xtick_pos = [2000, 4000, 8000]
        xtick_lbl = ['2000', '4000', '8000']
        ax2.set_xticks(xtick_pos)
        ax2.set_xticklabels(xtick_lbl)
        ax2.set_xlim([1000, 10000])
        ax2.set_ylim([0.01, 4])
        ax2.set_ylabel(r'$A_\lambda$ [mag]')
        ax2.set_xlabel(r'Wavelength [$\AA$]')

    def add_spec_plot(self, ax3):
        ax3.set_xscale('log')
        if self.dust_em_class.fixed:
            xtick_pos = [1000, 3000, 5000, 10000, 20000, 40000]
            xtick_lbl = ['0.1', '0.3', '0.5', '1', '2', '4']
            xlims = (1. + self.redshift) * np.array([1150, 20000])
            xlims[0] = min( xlims[0], min(self.fluxwv) - 200)
            xlims[1] = max( xlims[1], max(self.fluxwv) + 5000) 
        else:
            xtick_pos = [3000, 5000, 10000, 40000, 100000, 400000, 1000000]
            xtick_lbl = ['0.3', '0.5', '1', '4', '10', '40', '100']
            xlims = (1. + self.redshift) * np.array([1150, 700000])
            xlims[0] = min( xlims[0], min(self.fluxwv) - 200)
            xlims[1] = max( xlims[1], max(self.fluxwv) + 50000) 
            ax3.set_yscale('log')
        ax3.set_xticks(xtick_pos)
        ax3.set_xticklabels(xtick_lbl)
        ax3.set_xlim(xlims)
        ax3.set_xlabel(r'Wavelength [$\mu$m]')
        ax3.set_ylabel(r'$F_{\nu}$ [$\mu$Jy]')

    def add_subplots(self, ax1, ax2, ax3, nsamples, rndsamples=200):
        ''' Add Subplots to Triangle plot below '''
        sp, fn = ([], []) 
        for i in np.arange(rndsamples):
            ind = np.random.randint(0, nsamples.shape[0])
            self.set_class_parameters(nsamples[ind, :])
            self.sfh_class.plot(ax1, alpha=0.1)
            self.dust_abs_class.plot(ax2, self.wave, alpha=0.1)
            self.spectrum_plot(ax3, alpha=0.1)

        ax3.plot(self.wave, self.medianspec, color='dimgray')
        ax3.scatter(self.fluxwv, self.fluxfn, marker='x', s=200,
                    color='dimgray', zorder=8)

        # plot "truths" if in test mode
        if self.input_params is not None:
            self.set_class_parameters(self.input_params)
            self.sfh_class.plot(ax1, color='k', alpha=1.0)
            self.dust_abs_class.plot(ax2, self.wave, color='k', alpha=1.0)
            self.spectrum_plot(ax3, color='k', alpha=0.5)
        if self.true_fnu is not None:
            p = ax3.scatter(self.fluxwv, self.true_fnu, marker='o', s=150,
                            color=[0.216, 0.471, 0.749], zorder=9)
            p.set_facecolor('none')

        ax3.errorbar(self.fluxwv, self.data_fnu, yerr=self.data_fnu_e, fmt='s',
                     fillstyle='none', markersize=150,
                     color=[0.510, 0.373, 0.529], zorder=10)
        ax3.scatter(self.fluxwv, self.data_fnu, marker='s', s=150,facecolors='none',
                    edgecolors=[0.510, 0.373, 0.529], linewidths=2, zorder=10)        
        sel = np.where((self.fluxwv > ax3.get_xlim()[0]) * (self.fluxwv < ax3.get_xlim()[1]))[0]
        ax3min = np.percentile(self.data_fnu[sel][self.data_fnu[sel]>0.0], 1)
        ax3max = np.percentile(self.data_fnu[sel][self.data_fnu[sel]>0.0], 99)
        ax3ran = ax3max - ax3min
        if not self.dust_em_class.fixed: 
            ax3max = max(max(self.data_fnu),max(self.medianspec))
            ax3.set_ylim([ax3min*0.5, ax3max + 0.4 * ax3ran])
            ax3.set_xlim(right=max(max(self.fluxwv),max(self.wave)))
        else:
            ax3.set_ylim([ax3min - 0.4 * ax3ran, ax3max + 0.6 * ax3ran])
        ax3.text((1.+self.redshift)*1400, ax3max,
                 r'${\chi}_{\nu}^2 = $%0.2f' % self.chi2['rchi2'])


    def triangle_plot(self, outname, lnprobcut=7.5, imgtype='png'):
        ''' Make a triangle corner plot for samples from fit

        Input
        -----
        outname : string
            The triangle plot will be saved as "triangle_{outname}.png"
        lnprobcut : float
            Some of the emcee chains include outliers.  This value serves as
            a cut in log probability space with respect to the maximum
            probability.  For reference, a Gaussian 1-sigma is 2.5 in log prob
            space.
        imgtype : string
            The file extension of the output plot
        '''
        # Make selection for three sigma sample
        chi2sel = (self.samples[:, -1] >
                   (np.max(self.samples[:, -1], axis=0) - lnprobcut))
        nsamples = self.samples[chi2sel, :]
        o = 0 
        names = self.get_param_names()[o:]
        names.append('Log Mass')
        if self.dust_em_class.assume_energy_balance:
            names.append("Log Dust Mass")
        if self.input_params is not None:
            truths = self.input_params[o:]
        else:
            truths = None
        percentilerange = [p for i, p in enumerate(self.get_param_lims())
                           if i >= o] + [[7, 11]]
        percentilerange = [.95] * len(names)
        if self.dust_em_class.fixed: 
            numderpar = 3
        else: 
            if self.dust_em_class.assume_energy_balance:
                numderpar = 5
            else:
                numderpar = 4
        if self.dust_em_class.assume_energy_balance:
            indarr = np.concatenate((np.arange(o,len(nsamples[0])-numderpar),np.array([-2]))) 
        else:
            indarr = np.arange(o,len(nsamples[0])-numderpar) 
        fsgrad = 11+int(round(0.75*len(indarr)))
        fig = corner.corner(nsamples[:, indarr], labels=names,
                            range=percentilerange, color='rebeccapurple',
                            truths=truths, truth_color='gainsboro',
                            label_kwargs={"fontsize": fsgrad}, show_titles=True,
                            title_kwargs={"fontsize": fsgrad-2},
                            quantiles=[0.16, 0.5, 0.84], bins=20,
                            **{'plot_density':False,
                               'plot_datapoints':False,
                               'fill_contours': True,
                               'plot_contours': True,
                               'no_fill_contours': False})
        w = fig.get_figwidth()
        fig.set_figwidth(w-(len(indarr)-13)*0.025*w)

        # Adding subplots
        w = fig.get_figwidth()
        fig.set_figwidth(w-(len(indarr)-13)*0.025*w)
        ax1 = fig.add_subplot(3, 1, 1)
        ax1.set_position([0.7-0.02*(len(indarr)-5), 0.60+0.001*(len(indarr)-5), 
                          0.28+0.02*(len(indarr)-5), 0.15+0.001*(len(indarr)-5)])
        ax2 = fig.add_subplot(3, 1, 2)
        ax2.set_position([0.7+0.008*(15-len(indarr)), 0.39, 
                          0.28-0.008*(15-len(indarr)), 0.15])
        ax3 = fig.add_subplot(3, 1, 3)
        ax3.set_position([0.38-0.008*(len(indarr)-4), 0.82-0.001*(len(indarr)-4), 
                          0.60+0.008*(len(indarr)-4), 0.15+0.001*(len(indarr)-4)])
        self.add_sfr_plot(ax1)
        self.add_dust_plot(ax2)
        self.add_spec_plot(ax3)
        self.add_subplots(ax1, ax2, ax3, nsamples)

        for ax_loc in fig.axes:
            ax_loc.minorticks_on() 
            ax_loc.set_axisbelow('False')

        fig.savefig("%s.%s" % (outname, imgtype), dpi=150)
        plt.close(fig)

    def sample_plot(self, outname, imgtype='png'):
        ''' Make a sample plot

        Input
        -----
        outname : string
            The sample plot will be saved as "sample_{outname}.png"
        imgtype : string
            The file extension of the output plot

        '''
        # Make selection for three sigma sample
        names = self.get_param_names()
        if self.input_params is not None:
            truths = self.input_params
        else:
            truths = None
        fig, ax = plt.subplots(self.chain.shape[2], 1, sharex=True,
                               figsize=(5, 2*self.chain.shape[2]))
        for i, a in enumerate(ax):
            for chain in self.chain[:, :, i]:
                a.plot(chain, 'k-', alpha=0.3)
            a.set_ylabel(names[i])
            if truths is not None:
                a.plot([0, self.chain.shape[1]], [truths[i], truths[i]], 'r--')
            if i == len(ax)-1:
                a.set_xlabel("Step")

        for ax_loc in fig.axes:
            ax_loc.minorticks_on()

        fig.savefig("%s.%s" % (outname, imgtype))
        plt.tight_layout()
        plt.close(fig)

    def add_fitinfo_to_table(self, percentiles, start_value=3, lnprobcut=7.5,
                             numsamples=1000):
        ''' Assumes that "Ln Prob" is the last column in self.samples'''
        chi2sel = (self.samples[:, -1] >
                   (np.max(self.samples[:, -1], axis=0) - lnprobcut))
        nsamples = self.samples[chi2sel, :-1]
        n = len(percentiles)
        for i, per in enumerate(percentiles):
            for j, v in enumerate(np.percentile(nsamples, per, axis=0)):
                self.table[-1][(i + start_value + j*n)] = v
        return (i + start_value + j*n)

    def add_truth_to_table(self, truth, start_value):
        for i, tr in enumerate(truth):
            self.table[-1][start_value + i + 1] = tr


