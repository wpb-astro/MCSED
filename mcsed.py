""" SED fitting class using emcee for parameter estimation

    CURRENT LIMITATIONS:
        A) Constant metallicity for input SSP
        B) Dust Emission is ad hoc from Draine and Li (2007)

    OPTIONAL FITTED PARAMETERS:
        A) SFH
            a) tau_sfh, age, a, b, c
        B) Dust law
            b) tau_dust, delta, Eb

    OUTPUT PRODUCTS:
        A) XXX Plot

.. moduleauthor:: Greg Zeimann <gregz@astro.as.utexas.edu>

"""
#import matplotlib
#matplotlib.use("Agg")
import logging
import sfh
import dust_abs
import dust_emission
import ssp
import cosmology
import emcee
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import corner
import time
# WPBWPB delete astrpy table
from astropy.table import Table
from scipy.integrate import simps
from scipy.interpolate import interp1d

import numpy as np

#WPBWPB re organize the arguments (aesthetic purposes)
class Mcsed:
    def __init__(self, filter_matrix, ssp_spectra,
                 emlinewave, ssp_emline, ssp_ages, ssp_masses,
                 ssp_met, wave, sfh_class, dust_abs_class, dust_em_class,
                 data_fnu=None, data_fnu_e=None, 
                 data_emline=None, data_emline_e=None, emline_dict=None,
                 redshift=None,
                 filter_flag=None, input_spectrum=None, input_params=None,
                 sigma_m=0.1, nwalkers=40, nsteps=1000, true_fnu=None):
        ''' Initialize the Mcsed class.

        Init
        ----
        filter_matrix : numpy array (2 dim)
            The filter_matrix has rows of wavelength and columns for each
            filter (can be much larger than the filters used for fitting)
        ssp_spectra : numpy array (3 dim)
            single stellar population spectrum for each age in ssp_ages
            and each metallicity in ssp_met 
        emlinewave : numpy array (1 dim)
            Rest-frame wavelengths of requested emission lines (emline_dict)
            Corresponds to ssp_emline
        ssp_emline : numpy array (3 dim)
            Emission line SSP grid spanning emlinewave, age, metallicity
            Only includes requested emission lines (from emline_dict)
            Only used for calculating model emission line strengths
            Spectral units are ergs / s / cm2 at 10 pc
        ssp_ages : numpy array (1 dim)
            ages of the SSP models
        ssp_masses : numpy array (1 dim)
            remnant masses of the SSP models
        ssp_met : numpy array (1 dim)
            metallicities of the SSP models
        wave : numpy array (1 dim)
            wavelength for SSP models and all model spectra
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
            This is the input class for dust absorption.
        data_fnu : numpy array (1 dim)
            Photometry for data.  Length = (filter_flag == True).sum()
WPBWPB units + are dimensions correct??
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
        use_emline_flux : bool
            If emline_dict contains emission lines, set to True. Else, False
        redshift : float
            Redshift of the source
        filter_flag : numpy array (1 dim)
            Length = filter_matrix.shape[1], True for filters matching data
        input_spectrum : numpy array (1 dim)
            F_nu(wave) for input
        input_params : list
            input parameters for modeling.  Intended for testing fitting
            procedure.
        sigma_m : float
            Fractional error expected from the models.  This is used in
            the log likelihood calculation.  No model is perfect, and this is
            more or less a fixed parameter to encapsulate that.
        nwalkers : int
            The number of walkers for emcee when fitting a model
        nsteps : int
            The number of steps each walker will make when fitting a model
        true_fnu : WPBWPB FILL IN
WPBWPB: describe self.t_birth, set using args and units of Gyr
        '''
        # Initialize all argument inputs
        self.filter_matrix = filter_matrix
        self.ssp_spectra = ssp_spectra
        self.emlinewave = emlinewave
        self.ssp_emline = ssp_emline
        self.ssp_ages = ssp_ages
        self.ssp_masses = ssp_masses
        self.ssp_met = ssp_met
        self.wave = wave
        self.dnu = np.abs(np.hstack([0., np.diff(2.99792e18 / self.wave)]))
        self.sfh_class = getattr(sfh, sfh_class)()
        self.dust_abs_class = getattr(dust_abs, dust_abs_class)()
# WPBWPB: is ssp_class used?
        self.ssp_class = getattr(ssp, 'fsps_freeparams')()
        self.dust_em_class = getattr(dust_emission, dust_em_class)()
# WPBWPB: describe SSP, lineSSP in comments... 
# ssp_spectra span many metallicities, SSP only span ages
        self.SSP = None
        self.lineSSP = None
# WPBWPB is ssp_class still a thing?
        self.param_classes = ['sfh_class', 'dust_abs_class', 'ssp_class',
                              'dust_em_class']
        self.data_fnu = data_fnu
        self.data_fnu_e = data_fnu_e
        self.data_emline = data_emline
        self.data_emline_e = data_emline_e
        self.emline_dict = emline_dict
        self.redshift = redshift
        self.filter_flag = filter_flag
        self.input_spectrum = input_spectrum
        self.input_params = input_params
        self.sigma_m = sigma_m
        self.nwalkers = nwalkers
        self.nsteps = nsteps
        self.true_fnu = true_fnu
        if self.redshift is not None:
            self.set_new_redshift(self.redshift)

        # Set up logging
        self.setup_logging()

        # Time array for sfh
        self.age_eval = np.logspace(-3, 1, 4000)

    def set_new_redshift(self, redshift):
        ''' Setting redshift

        Parameters
        ----------
        redshift : float
            Redshift of the source for fitting
        '''
        self.redshift = redshift
        # Need luminosity distance to adjust ssp_spectra from 10pc to Dl
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
        ''' FILL IN
        '''
        wave_avg = np.dot(self.wave, self.filter_matrix[:, self.filter_flag])
        return wave_avg

    def get_filter_fluxdensities(self):
        '''Convert a spectrum to photometric fluxes for a given filter set.
        The photometric fluxes will be in the same units as the spectrum.
        The spectrum is in microjanskies(lambda) such that
        the photometric fluxes will be in microjanskies.

        Returns
        -------
        f_nu : numpy array (1 dim)
            Photometric flux densities for an input spectrum
        '''
## WPBWPB delete
#        print('shape of spectrum, filter_matrix, filter_flag:')
#        print((self.spectrum.shape, self.filter_matrix.shape, self.filter_flag.shape))
        f_nu = np.dot(self.spectrum, self.filter_matrix[:, self.filter_flag])
        return f_nu

    def set_class_parameters(self, theta):
        ''' For a given set of model parameters, set the needed class variables
        related to SFH, dust attenuation, ect.

        Input
        -----
        theta : list
            list of input parameters for sfh, dust att., and dust em.
        '''
        start_value = 0
        ######################################################################
        # STAR FORMATION HISTORY
        self.sfh_class.set_parameters_from_list(theta, start_value)
        # Keeping track of theta index for age of model and other classes
        start_value += self.sfh_class.get_nparams()

        ######################################################################
        # DUST ATTENUATION
        self.dust_abs_class.set_parameters_from_list(theta, start_value)
        start_value += self.dust_abs_class.get_nparams()
# WPBWPB modify: pass a dust_abs_birthcloud keyword, see if its the same, blah

        ######################################################################
        # SSP Parameters
## WPBWPB delete
#        print(start_value)
#        print(self.ssp_class.fix_met)
        self.ssp_class.set_parameters_from_list(theta, start_value)
        start_value += self.ssp_class.get_nparams()
## WPBWPB delete
#        print(start_value)
        ######################################################################
        # DUST EMISSION
        self.dust_em_class.set_parameters_from_list(theta, start_value)
        start_value += self.dust_em_class.get_nparams()

    def get_ssp_spectrum(self):
        '''
        Calculate SSP for an arbitrary metallicity (self.ssp_class.met) given a
        model grid for a range of metallicities (self.ssp_met)

        Returns
        -------
        SSP : 2-d array
            Single stellar population models for each age in self.ages
        '''
        if self.ssp_class.fix_met:
            if self.SSP is not None:
## WPBWPB delete
#                print('self.SSP is not None!')
                return self.SSP, self.lineSSP
        Z = np.log10(self.ssp_met)
        z = self.ssp_class.met + np.log10(0.019)
        X = Z - z
        wei = np.exp(-(X)**2 / (2. * 0.15**2))
        wei /= wei.sum()
        self.SSP = np.dot(self.ssp_spectra, wei)
# WPBWPB: this is where I would relax logU criteria, same as metallicity
# careful: needs to be dealt with when measuring line fluxes originally

        # only treat the emission line grid if requested
## WPBWPB delete
#        print('shape of emline SSP before/after wei')
#        print(self.ssp_emline.shape)
        if self.use_emline_flux:
            self.lineSSP = np.dot(self.ssp_emline, wei)
        else:
            self.lineSSP = self.ssp_emline[:,:,0]
## WPBWPB delete
#        print(self.lineSSP.shape)
        return self.SSP, self.lineSSP

    def build_csp(self, sfr=None):
        '''Build a composite stellar population model for a given star
        formation history, dust attenuation law, and dust emission law.

        In addition to the returns it also modifies a lineflux dictionary

        Returns
        -------
        csp : numpy array (1 dim)
            Composite stellar population model at self.redshift
WPBWPB units??
        mass : float
            Mass for csp given the SFH input
        '''
        # Collapse for metallicity
        SSP, lineSSP = self.get_ssp_spectrum()

        # Need star formation rate from observation back to formation
        if sfr is None:
            sfr = self.sfh_class.evaluate(self.ssp_ages)
        ageval = 10**self.sfh_class.age # Gyr

## WPBWPB delete
#        print('these are the sfr in SFH age grid in Gyr:')
#        sfh_ages = 10.**(np.array(self.sfh_class.ages)-9.)
#        print(sfh_ages)
#        print(self.sfh_class.evaluate(sfh_ages))

        # Treat the birth cloud and diffuse component separately
# WPBWPB: may want to modify: have this as user-defined setting...
        age_birth = self.t_birth #10**-2 # Gyr 

### WPBWPB delete - both in Gyr
##        print('this is the ageval: %s' % ageval)
##        print('this is ssp ages: %s' % self.ssp_ages)
##        return

## WPBWPB delete
##        for ageval, age_birth in [ [0.008, 0.011 ], [0.005011872336272725, 0.011 ], [0.01, 0.011 ], [0.14, 0.011 ], [0.0116, 0.011 ], [0.1, 0.011 ], [0.0145, 0.01 ], [0.011, 0.011 ], [0.01, 0.01 ] ]:
#        for ageval, age_birth in [ [0.01116, 0.011] ]:
#            self.build_dustfree_CSP(sfr, ageval, age_birth)
#        return

        # Get dust-free CSPs, properly accounting for ages
        # ageval sets limit on ssp_ages that are useable in model calculation
        # age_birth separates birth cloud and diffuse components
# WPBWPB delete -- ageval, ssp_ages, age_birth are in units Gyr
        sel = (self.ssp_ages > age_birth) & (self.ssp_ages <= ageval)
        sel_birth = (self.ssp_ages <= age_birth) & (self.ssp_ages <= ageval)
        sel_age = self.ssp_ages <= ageval

        # The weight is the time between ages of each SSP
        weight = np.diff(np.hstack([0, self.ssp_ages])) * 1e9 * sfr
# WPBWPB delete
        weight_orig = weight.copy()
        weight_birth = weight.copy()
        weight_age = weight.copy()
        # Ages greater than ageval should have zero weight in CSP
        # weight should only include populations younger than ageval
        # and older than age_birth
        # weight_birth should only include populations younger than ageval
        # and no older than age_birth
        # weight_age only considers the age of the system (for mass)
        weight[~sel] = 0
        weight_birth[~sel_birth] = 0
        weight_age[~sel_age] = 0

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
        if (len(select_too_old)>0): # & (ageval>=age_birth):
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
        #print "Max(SSP) = %.3e"%(np.amax(self.SSP))
        spec_dustfree = np.dot(self.SSP, weight)
        spec_birth_dustfree = np.dot(self.SSP, weight_birth)
        linespec_dustfree = np.dot(self.lineSSP, weight_birth)
        mass = np.sum(weight_age * self.ssp_masses)

## WPBWPB delete
#        print('These are the weights, and various properties:')
#        t = Table()
#        t['ageval'] = [ageval]*len(weight)
#        t['age_birth'] = [age_birth]*len(weight)
#        t['ssp_ages'] = self.ssp_ages
#        t['sfr'] = sfr
#        t['weight_orig'] = weight_orig
#        t['weight_young'] = weight_birth
#        t['weight_old'] = weight
#        t['weight_age'] = weight_age
#        t.write('CSP_weights.dat',format='ascii')
#        return

        # Need to correct spectrum for dust attenuation
        Alam = self.dust_abs_class.evaluate(self.wave)
        spec_dustobscured = spec_dustfree * 10**(-0.4 * Alam)

        # Correct the corresponding birth cloud spectrum separately
# WPBWPB: check which law using for the birth cloud
# if attenuating it directly tied to overall dust law,
# modified by coefficient between EBV_stars ~ gas, get it here
        Alam_birth = Alam / self.dust_abs_class.EBV_stars_gas
        spec_birth_dustobscured = spec_birth_dustfree * 10**(-0.4 * Alam_birth)

        # Combine the young and old components
        spec_dustfree += spec_birth_dustfree
        spec_dustobscured += spec_birth_dustobscured

        # compute attenuation for emission lines
        Alam_emline = (self.dust_abs_class.evaluate(self.emlinewave,new_wave=True)
                       / self.dust_abs_class.EBV_stars_gas)
# WPBWPB: else, use a separate model

## WPBWPB delete
#        print(emwaves)
#        print(Alam_emline)
#        print('shape of linespec_dustfree, emwaves : (%s, %s)' % (linespec_dustfree.shape, emwaves.shape))
        linespec_dustobscured = linespec_dustfree * 10**(-0.4*Alam_emline)

## WPBWPB compare Alam in diffuse, birthcloud components
#        print('diffuse, birth, emline Alam:')
#        print(Alam[ np.searchsorted(self.wave, self.emlinewave) ])
#        print(Alam_birth[ np.searchsorted(self.wave, self.emlinewave) ])
#        print(Alam_emline)


# WPB: exclude dust emission component altogether? Does it make a difference?
        # Change in bolometric Luminosity
        # L_bol = (np.dot(self.dnu, spec_dustfree) -
        #          np.dot(self.dnu, spec_dustobscured))
        #print "Max(spec_dustfree) = %.3e"%(max(spec_dustfree))
        #print "Max(spec_dustobscured) = %.3e"%(max(spec_dustobscured))
        #print "L_bol = %.3e"%(L_bol)
        # if not self.dust_em_class.fixed: 
        #     umin,gamma,qpah = self.dust_em_class.get_params()
        # else:
        #     umin,gamma,qpah = 2.0, 0.05, 2.5 #Default values
        # umax=1.0e6
        # P0 = 135.0 #Power absorbed per unit dust mass in radiation field U=1; units L_sun/M_sun
        # Lbolfac = 2.488e-24 #Convert from uJy*Hz at 10 pc to L_sun
        # uavg = (1.-gamma)*umin + gamma*umin*np.log(umax/umin) / (1.-umin/umax)
        # mdust = L_bol*Lbolfac/uavg/P0
        # if L_bol<0 or uavg<0 or ~np.isfinite(L_bol) or ~np.isfinite(uavg):
        #     print "Lbol = %.3e; <U> = %.3f; M_dust = %.3e"%(L_bol*Lbolfac,uavg,mdust)

        # Add dust emission
        if min(spec_dustobscured[self.wave>5.0e4])<0.0: 
            print("Before adding dust: min(spec_dustobscured[wave>5.0 um]) =",
                  min(spec_dustobscured[self.wave>5.0e4]))
        spec_dustobscured += self.dust_em_class.evaluate(self.wave)
        #print("After adding dust: min(spec_dustobscured[wave>5.0 um]) =",min(spec_dustobscured[self.wave>5.0e4]))

        # Redshift to observed frame
        csp = np.interp(self.wave, self.wave * (1. + self.redshift),
                        spec_dustobscured * (1. + self.redshift))

        # Update dictionary of modeled emission line fluxes
        linefluxCSPdict = {}
        if self.use_emline_flux:
            for emline in self.emline_dict.keys():
                indx = np.argmin(np.abs(self.emlinewave 
                                        - self.emline_dict[emline][0]))
                # flux is given in ergs / s / cm2 at 10 pc
                flux = linespec_dustobscured[indx]
                # Correct flux from 10pc to redshift of source
                linefluxCSPdict[emline] = linespec_dustobscured[indx] / self.Dl**2
        self.linefluxCSPdict = linefluxCSPdict

## WPBWPB delete
#        print( linefluxCSPdict )

        # Correct spectra from 10pc to redshift of the source
        return csp / self.Dl**2, mass
        # if not DMreturn: 
        #     return csp / self.Dl**2, mass
        # else:
        #     return csp/self.Dl**2, mass, mdust, dustmass*dusttohratio

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
        log likelihood, mass, t10, t50, t90, sfr10, sfr100 : float, float, float, float, float, float, float
            The log likelihood includes a chi2_term and a parameters term.
            The mass comes from building of the composite stellar population
            The parameters t10, t50, t90, sfr10, and sfr100 are derived in get_derived_params(self)
        '''
        # if not self.dust_em_class.fixed: 
        #     self.spectrum, mass, mdust, mdust2 = self.build_csp(DMreturn=True)
        # else:
        #     self.spectrum, mass = self.build_csp(DMreturn=False)
        #     mdust = None
        #     mdust2 = None
        self.spectrum, mass = self.build_csp()

        # compare input and model emission line fluxes
        emline_term = 0.0
        if self.use_emline_flux:
            # if all lines have null line strengths, ignore 
            if not min(self.data_emline) == max(self.data_emline) == -99:
## WPBWPB delete
#                print('this is emline_dict:' + str(self.emline_dict))
#                print('this is emline data, error:')
#                print(self.data_emline)
#                print(self.data_emline_e)
                for emline in self.emline_dict.keys():
                    if self.data_emline['%s_FLUX' % emline] > -99: # null value
                        emline_wave, emline_weight = self.emline_dict[emline]
                        model_lineflux = self.linefluxCSPdict[emline] 
                        lineflux  = self.data_emline['%s_FLUX' % emline]
                        elineflux = self.data_emline_e['%s_ERR' % emline]
                        emline_term += (-0.5 * (model_lineflux - lineflux)**2 /
                                        elineflux**2.) * emline_weight
## WPBWPB: delete
#                print('this is emline and term:')
#                print(emline)
#                print(emline_term)

        model_y = self.get_filter_fluxdensities()
        inv_sigma2 = 1.0 / (self.data_fnu_e**2 + (model_y * self.sigma_m)**2)
        chi2_term = -0.5 * np.sum((self.data_fnu - model_y)**2 * inv_sigma2)
        parm_term = -0.5 * np.sum(np.log(1 / inv_sigma2))
        sfr10,sfr100,fpdr = self.get_derived_params()
        #return (chi2_term + parm_term + emline_term, mass,sfr10,sfr100,fpdr,mdust,mdust2)
        return (chi2_term + parm_term + emline_term, mass,sfr10,sfr100,fpdr)

    def lnprob(self, theta):
        ''' Calculate the log probabilty and return the value and stellar mass (as well as derived parameters)
        of the model

        Returns
        -------
        log prior + log likelihood, [mass,sfr10,sfr100,fpdr]: float,float,float,float,float
            The log probability is just the sum of the logs of the prior and
            likelihood.  The mass comes from the building of the composite
            stellar population. The other derived parameters are calculated in get_derived_params()
        '''
        self.set_class_parameters(theta)
        lp = self.lnprior()
        if np.isfinite(lp):
            #lnl,mass,sfr10,sfr100,fpdr,mdust,mdust2 = self.lnlike()
            lnl,mass,sfr10,sfr100,fpdr = self.lnlike()
            if not self.dust_em_class.fixed:
                #return lp + lnl, np.array([mass, sfr10, sfr100, fpdr, mdust, mdust2])
                return lp + lnl, np.array([mass, sfr10, sfr100, fpdr])
            else:
                return lp + lnl, np.array([mass, sfr10, sfr100])
        else:
            if not self.dust_em_class.fixed:
                #return -np.inf, np.array([-np.inf, -np.inf, -np.inf, -np.inf, -np.inf, -np.inf])
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
            pos = (np.random.rand(num)[:, np.newaxis] *
                   (theta_lims[:, 1]-theta_lims[:, 0]) + theta_lims[:, 0])
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
## WPBWPB delete
#            print(par_cl)
#            print(vals)
        vals = list(np.hstack(vals))
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
        data magnitudes and errors
        '''
        # Need to verify data parameters have been set since this is not
        # a necessity on initiation
        self.log.info('Fitting model using emcee')
        check_vars = ['data_fnu', 'data_fnu_e', 'redshift', 'filter_flag']
        for var in check_vars:
            if getattr(self, var) is None:
                self.error('The variable %s must be set first' % var)

        pos = self.get_init_walker_values(kind='ball')
        ndim = pos.shape[1]
        start = time.time()
        # Time to set up the sampler and run the mcmc
        #dtype = [("log_mass", float),("t10", float),("t50", float),("t90", float),("sfr10", float),("sfr100", float)] #For blobs
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
            numderpar = 4
        new_chain = np.zeros((self.nwalkers, self.nsteps, ndim+numderpar+1))
        new_chain[:, :, :-(numderpar+1)] = sampler.chain
        self.chain = sampler.chain
        for i in xrange(len(sampler.blobs)):
            for j in xrange(len(sampler.blobs[0])):
                for k in xrange(len(sampler.blobs[0][0])):
                    x = sampler.blobs[i][j][k]
                    #if k==0 or k==4 or k==5: #Stellar mass or Dust mass--can't take log of negative numbers
                    if k==0: #Stellar mass--can't take log of negative numbers; also, want reasonable values
                        new_chain[j, i, -(numderpar+1)+k] = np.where((np.isfinite(x)) * (x > 10.),
                                               np.log10(x), -99.) #Stellar mass
                    else: 
                        new_chain[j, i, -(numderpar+1)+k] = np.where((np.isfinite(x)),np.log10(x), -99.) #Other derived params
        new_chain[:, :, -1] = sampler.lnprobability
        self.samples = new_chain[:, burnin_step:, :].reshape((-1, ndim+numderpar+1))

    def calc_gaw(self,t,sfr,frac,mass):
        ''' Calculate lookback time at which the fraction "frac" of the stellar mass in the galaxy was created'''
        ind = 0
        #print "Length of t =", len(t)
        #print "Total SFR integral over mass =", simps(sfr,t)/mass
        #stellartot = quad(sfr_f,t[0],t[-1])[0]
        #while quad(sfr_f,t[0],t[ind])[0]/mass < frac: 
        #print "SFR stats: Max = %.3e, Min = %.3e"%(max(sfr),min(sfr))
        #indage = (np.abs(t-10**self.sfh_class.age)).argmin()
        #print "indage =", indage
        #totmass = simps(sfr[:indage+1],t[:indage+1])/mass
        #print "Age: %.3f, Mass: %.3f, SFR Integral over time / Mass: %.3f"%(self.sfh_class.age,np.log10(mass),totmass)
        while (simps(sfr[:ind+1],t[:ind+1])/mass<(1.0-frac) and ind<len(t)-1):
            ind+=1
        #forward = quad(sfr_f,t[0],t[ind])[0]/mass - frac
        forward = abs(simps(sfr[:ind+1],t[:ind+1])/mass - 1.0 + frac) #This SHOULD be positive even without abs()
        #backward = frac - quad(sfr_f,t[0],t[ind-1])[0]/mass
        backward = abs(1.0-frac - simps(sfr[:ind],t[:ind])/mass) #This SHOULD be positive even without abs()
        tot = forward+backward
        #print "Time given for txx: %.3f"%(np.log10(1.0e-9*(forward/tot *t[ind-1] + backward/tot * t[ind])))
        return 1.0e-9*(forward/tot *t[ind-1] + backward/tot * t[ind]) #Linear interpolation to get more accurate result--put result back into Gyr; Note--this may not be as accurate IF either forward or backward were negative before the absolute value--but even in those cases, the solution should be quite close. The only case where this should theoretically happen is if the integral of SFR over time doesn't go over (1-frac)*mass before the age of the universe at the time of observation is reached--this should really hopefully not happen.

    def get_derived_params1(self,params,mass):
        ''' These are not free parameters in the model, but are instead
        calculated from free parameters
        '''
        #if agenum is not None: age = params[agenum]
        #else: age = self.sfh_class.age
        #ageval = 10**age #Age in Gyr
        C = cosmology.Cosmology()
        ageval = C.lookback_time(20)-C.lookback_time(self.redshift) #Don't want to restrict txx parameters to estimated age
        t_sfr100 = np.linspace(1.0e-9,0.1,num=1001) #From 100 Mya to present (observed time); avoid t=0 for log purposes
        t_sfr10 = np.linspace(1.0e-9,0.01,num=1001) #From 10 Mya to present (observed time); avoid t=0 for log purposes
        sfrarray = self.sfh_class.evaluate(t_sfr100,force_params=params)
        sfr100 = simps(sfrarray,x=t_sfr100)/(t_sfr100[-1]-t_sfr100[0]) #Mean value over last 100 My
        sfrarray = self.sfh_class.evaluate(t_sfr10,force_params=params)
        sfr10 = simps(sfrarray,x=t_sfr10)/(t_sfr10[-1]-t_sfr10[0]) #Mean value over last 10 My
        t_gaw = np.geomspace(1.0e-9,ageval,num=250) #From present day (basically) back to birth of galaxy in Gyr
        sfrfull = self.sfh_class.evaluate(t_gaw,force_params=params)
        t_gaw*=1.0e9 #Need it in years for calculation
        #sfr_f = interp1d(t_gaw,sfrfull,kind='cubic',fill_value="extrapolate")

        t10 = self.calc_gaw(t_gaw,sfrfull,0.1,mass)
        t50 = self.calc_gaw(t_gaw,sfrfull,0.5,mass)
        t90 = self.calc_gaw(t_gaw,sfrfull,0.9,mass)

        return [sfr10,sfr100,t10,t50,t90]

    def get_t_params(self,params,mass):
        #if agenum is not None: age = params[agenum]
        #else: age = self.sfh_class.age
        #ageval = 10**age #Age in Gyr
        C = cosmology.Cosmology()
        ageval = C.lookback_time(20)-C.lookback_time(self.redshift) #Don't want to restrict txx parameters to estimated age
        t_gaw = np.geomspace(1.0e-9,ageval,num=250) #From present day (basically) back to birth of galaxy in Gyr
        sfrfull = self.sfh_class.evaluate(t_gaw,force_params=params)
        #sfrfull_avg = self.sfh_class.evaluate(t_gaw)
        #print "Fractional difference between average sfr over time and this particular set of sfh params:", np.linalg.norm(sfrfull-sfrfull_avg)/np.linalg.norm(sfrfull_avg)
        t_gaw*=1.0e9 #Need it in years for calculation
        t10 = self.calc_gaw(t_gaw,sfrfull,0.1,mass)
        t50 = self.calc_gaw(t_gaw,sfrfull,0.5,mass)
        t90 = self.calc_gaw(t_gaw,sfrfull,0.9,mass)
        return np.log10(t10),np.log10(t50),np.log10(t90)


    def get_derived_params(self):
        ''' These are not free parameters in the model, but are instead
        calculated from free parameters
        '''
        ageval = 10**self.sfh_class.age #Age in Gyr
        t_sfr100 = np.linspace(1.0e-9,0.1,num=251) #From 100 Mya to present (observed time); avoid t=0 for log purposes
        t_sfr10 = np.linspace(1.0e-9,0.01,num=251) #From 10 Mya to present (observed time); avoid t=0 for log purposes
        sfrarray = self.sfh_class.evaluate(t_sfr100)
        sfr100 = simps(sfrarray,x=t_sfr100)/(t_sfr100[-1]-t_sfr100[0]) #Mean value over last 100 My
        sfrarray = self.sfh_class.evaluate(t_sfr10)
        sfr10 = simps(sfrarray,x=t_sfr10)/(t_sfr10[-1]-t_sfr10[0]) #Mean value over last 10 My

        if self.dust_em_class.fixed:
            fpdr = None
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
            median spectrum
        '''
        chi2sel = (self.samples[:, -1] >
                   (np.max(self.samples[:, -1], axis=0) - lnprobcut))
        nsamples = self.samples[chi2sel, :]
        wv = self.get_filter_wavelengths()
## WPBWPB delete
#        rndsamples = 200
        sp, fn = ([], [])
        #start = time.time()
        for i in np.arange(rndsamples):
        #for i in range(len(nsamples)):
            ind = np.random.randint(0, nsamples.shape[0])
            self.set_class_parameters(nsamples[ind, :])
            #self.set_class_parameters(nsamples[i, :])
            self.spectrum, mass = self.build_csp()
            fnu = self.get_filter_fluxdensities()
            sp.append(self.spectrum * 1.)
            fn.append(fnu * 1.)
            #if i%(len(nsamples)/20)==0: print("Finished %d/20 of making spectra"%(20*i/len(nsamples)))
        self.medianspec = np.median(np.array(sp), axis=0)
        self.fluxwv = wv
        self.fluxfn = np.median(np.array(fn), axis=0)
        #end = time.time()
        #elapsed = end-start
        #self.log.info("Total time taken for creating medianspec: %0.2f s" % elapsed)


    def spectrum_plot(self, ax, color=[0.996, 0.702, 0.031], alpha=0.1):
        ''' Make spectum plot for current model '''
        self.spectrum, mass = self.build_csp()
        ax.plot(self.wave, self.spectrum, color=color, alpha=alpha)

    def add_sfr_plot(self, ax1):
        ax1.set_xscale('log')
        ax1.set_yscale('log')
        ax1.set_ylabel(r'SFR $M_{\odot} yr^{-1}$')
        ax1.set_xlabel('Lookback Time (Gyr)')
        ax1.set_xticks([1e-3, 1e-2, 1e-1, 1])
        ax1.set_xticklabels(['1 Myr', '10 Myr', '100 Myr', '1 Gyr'])
        ax1.set_yticks([1e-3, 1e-2, 1e-1, 1, 1e1, 1e2, 1e3])
        ax1.set_yticklabels(['0.001', '0.01', '0.1', '1', '10', '100', '1000'])
        ax1.set_xlim([10**-3, 10**self.sfh_class.age_lims[1]])
        ax1.set_ylim([1e-3, 1e3])

    def add_dust_plot(self, ax2):
        ax2.set_xscale('log')
        xtick_pos = [1000, 3000, 10000]
        xtick_lbl = ['1000', '3000', '10000']
        ax2.set_xticks(xtick_pos)
        ax2.set_xticklabels(xtick_lbl)
        ax2.set_xlim([1000, 20000])
        ax2.set_ylim([0, 8])
        ax2.set_ylabel(r'Dust Attenuation (mag)')
        ax2.set_xlabel(r'Wavelength $\AA$')

    def add_spec_plot(self, ax3):
# WPBWPB: adjust wavelength range, depending on whether dust emission is fit
        ax3.set_xscale('log')
        if self.dust_em_class.fixed:
            xtick_pos = [3000, 5000, 10000, 20000, 40000]
            xtick_lbl = ['0.3', '0.5', '1', '2', '4']
            xlims = [3000, 50000]
        else:
            xtick_pos = [3000, 5000, 10000, 40000, 100000, 1000000]
            xtick_lbl = ['0.3', '0.5', '1', '4', '10', '100']
            xlims = [3000, 2000000]
            ax3.set_yscale('log')
        ax3.set_xticks(xtick_pos)
        ax3.set_xticklabels(xtick_lbl)
        ax3.set_xlim(xlims)
        ax3.set_xlabel(r'Wavelength $\mu m$')
        ax3.set_ylabel(r'$F_{\nu}$ ($\mu$Jy)')

    def add_subplots(self, ax1, ax2, ax3, nsamples, rndsamples=200):
        ''' Add Subplots to Triangle plot below '''
### WPBWPB -- I think all of this can be deleted - migrated to a separate method that does not require plotting (double-commented lines within this method)
##        wv = self.get_filter_wavelengths()
##        rndsamples = 200
        sp, fn = ([], []) 
        for i in np.arange(rndsamples):
            ind = np.random.randint(0, nsamples.shape[0])
            self.set_class_parameters(nsamples[ind, :])
            self.sfh_class.plot(ax1, alpha=0.1)
            self.dust_abs_class.plot(ax2, self.wave, alpha=0.1)
            self.spectrum_plot(ax3, alpha=0.1)

##            fnu = self.get_filter_fluxdensities()
##            sp.append(self.spectrum * 1.)
##            fn.append(fnu * 1.)
## WPB edit: plotting HBeta line
## used to have self.hbflux = self.measure_hb() --> changed
##            hbm.append(self.hbflux * 1.)
##        # Plotting median value:
##        self.medianspec = np.median(np.array(sp), axis=0)
###        self.hbmedian = np.median(hbm)
        ax3.plot(self.wave, self.medianspec, color='dimgray')
##        self.fluxwv = wv
##        self.fluxfn = np.median(np.array(fn), axis=0)

        ax3.scatter(self.fluxwv, self.fluxfn, marker='x', s=200,
                    color='dimgray', zorder=8)
        chi2 = (1. / (len(self.data_fnu) - 1) *
                (((self.data_fnu - self.fluxfn) / self.data_fnu_e)**2).sum())
        # WPBWPB: reduced chi^2 or not? properly accounting for number of data points, including emlines?
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
                     fillstyle='none', markersize=15,
                     color=[0.510, 0.373, 0.529], zorder=10)
        
        sel = np.where((self.fluxwv > 3000.) * (self.fluxwv < 50000.))[0]
        ax3min = np.percentile(self.data_fnu[sel][self.data_fnu[sel]>0.0], 5)
        ax3max = np.percentile(self.data_fnu[sel][self.data_fnu[sel]>0.0], 95)
        ax3ran = ax3max - ax3min
        if not self.dust_em_class.fixed: 
            ax3max = max(max(self.data_fnu),max(self.medianspec))
            ax3.set_ylim([ax3min*0.5, ax3max + 0.4 * ax3ran])
            ax3.set_xlim(right=max(max(self.fluxwv),max(self.wave)))
        else:
            ax3.set_ylim([ax3min - 0.4 * ax3ran, ax3max + 0.4 * ax3ran])
        ax3.text(4200, ax3max + 0.2 * ax3ran, r'${\chi}_{\nu}^2 = $%0.2f' % chi2)

    def triangle_plot(self, outname, lnprobcut=7.5, imgtype='png'):
        ''' Make a triangle corner plot for samples from fit
        * Doesn't include the derived parameters t10, t50, t90, sfr10, and sfr100 as that would make the plot too crowded

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
# WPBWPB: since median fits have already been set, may be able to remove lnprobcut -- just use attributes that are already set.
        # Make selection for three sigma sample
        chi2sel = (self.samples[:, -1] >
                   (np.max(self.samples[:, -1], axis=0) - lnprobcut))
        nsamples = self.samples[chi2sel, :]
# WPBWPB: understand this line....
        o = 0  # self.sfh_class.nparams
        names = self.get_param_names()[o:]
        names.append('Log Mass')
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
            numderpar = 4
        print("I'm starting to construct the triangle plot")
        fig = corner.corner(nsamples[:, o:-numderpar], labels=names,
                            range=percentilerange,
                            truths=truths, truth_color='gainsboro',
                            label_kwargs={"fontsize": 18}, show_titles=True,
                            title_kwargs={"fontsize": 16},
                            quantiles=[0.16, 0.5, 0.84], bins=30)
        print('made the corner')
        # Adding subplots
        ax1 = fig.add_subplot(3, 1, 1)
        ax1.set_position([0.7, 0.60, 0.25, 0.15])
        ax2 = fig.add_subplot(3, 1, 2)
        ax2.set_position([0.7, 0.40, 0.25, 0.15])
        ax3 = fig.add_subplot(3, 1, 3)
        ax3.set_position([0.38, 0.80, 0.57, 0.15])
        self.add_sfr_plot(ax1)
# WPBWPB delete:
        print("I've added the sfr plot")
        self.add_dust_plot(ax2)
# WPBWPB delete:
        print("I've added the dust plot")
        self.add_spec_plot(ax3)
# WPBWPB delete:
        print("I've added the spec plot")
        self.add_subplots(ax1, ax2, ax3, nsamples)
# WPBWPB delete:
        print("I've added the subplots")
# WPB edit: printing HBeta line flux on the figure
# used to have self.hbflux = self.measure_hb() --> changed
#        if self.sfh_class.hblim is not None:
#            fig.text(.5, .75, r'H$\beta$ input: %0.2f' %
#                     (self.sfh_class.hblim * 1e17), fontsize=18)
#        fig.text(.5, .70, r'H$\beta$ model: %0.2f' % (self.hbmedian * 1e17),
#                 fontsize=18)
        # fig.set_size_inches(15.0, 15.0)
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
        fig.savefig("%s.%s" % (outname, imgtype))
        plt.close(fig)

    def add_fitinfo_to_table(self, percentiles, start_value=3, lnprobcut=7.5,
                             numsamples=1000,numder=3):
        ''' Assumes that "Ln Prob" is the last column in self.samples'''
        if self.dust_em_class.fixed: 
            numderpar = 3
        else: 
            numderpar = 4
        chi2sel = (self.samples[:, -1] >
                   (np.max(self.samples[:, -1], axis=0) - lnprobcut))
        nsamples = self.samples[chi2sel, :-1]
        sfhnum = self.sfh_class.get_nparams()
        sfhnames = self.sfh_class.get_names()
        # if "Log Age" in sfhnames: 
        #     agenum = sfhnames.index("Log Age")
        # else: 
        #     agenum=None
        params = np.zeros((sfhnum,numsamples))
        t10,t50,t90 = np.zeros(numsamples),np.zeros(numsamples),np.zeros(numsamples)
        derpar = np.zeros((numsamples,numder))
        ranarray = np.random.choice(len(nsamples),size=numsamples) #Make sure mass and params all chosen from same random set of nsamples
        mass = nsamples[:,-numderpar][ranarray]
        mass = 10**mass #It was in log units before--we want to feed the function get_t_params linear units
        for k in range(sfhnum): #Get random values of SFH parameters based on their distributions
            #params[k] = np.random.choice(nsamples[:,k],size=numsamples)
            params[k] = nsamples[:,k][ranarray]

        for k2 in range(numsamples):
            derpar[k2] = self.get_t_params(params[:,k2],mass[k2])
            if k2%(numsamples/10)==0: print k2,params[:,k2],mass[k2],derpar[k2]

        n = len(percentiles)
        for i, per in enumerate(percentiles):
            for j, v in enumerate(np.percentile(nsamples, per, axis=0)):
                self.table[-1][(i + start_value + j*n)] = v
        current = i+start_value+j*n
        for i,per in enumerate(percentiles):
            for j,v in enumerate(np.percentile(derpar,per,axis=0)):
                self.table[-1][(current+i+j*n+1)] = v
        return (i + current + j*n+1)

    def add_truth_to_table(self, truth, start_value):
## WPBWPB generalize: what if sfh parameters are not first? will that ever occur? ensure working for general case
        sfhnum = self.sfh_class.get_nparams()
        sfhtruth = truth[:sfhnum]
        mass = 10**truth[-1] #Stellar mass
        derpar = self.get_derived_params1(sfhtruth,mass)
        for par in derpar: 
            truth.append(np.log10(par))
        for i, tr in enumerate(truth):
            self.table[-1][start_value + i + 1] = tr
        #last = start_value+i+1
        #for j in range(len(derpar)):
         #   self.table[-1][last+j+1]=derpar[j]
