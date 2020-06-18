""" MCSED - ssp.py

Single Stellar Population module for loading models

.. moduleauthor:: Greg Zeimann <gregz@astro.as.utexas.edu>

"""
import sfh
import sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import os.path as op
import scipy.interpolate as scint
from astropy.convolution import Gaussian1DKernel, convolve

plt.ioff() # UNCOMMENT THIS LINE IF RUNNING ON LINUX

def bin_ssp_ages(ssp_ages, ssp_spec, ssp_linespec, sfh_ages, galaxy_age, t_birth):
    '''
    Collapse the SSP grid into fewer ages, if possible, to significantly
    increase the computational efficiency.

    If using a binned SFH, the SFR within each age bin is assumed to 
    be constant, allowing the SSP spectra within the bin to be collapsed.

    The SSP spectra within each bin are combined via weighted average,
    where the weights are determined by the amount of time between the 
    SSP ages that fall in the bin.

    The last age will extend to the oldest possible age of the galaxy.

    Parameters
    ----------
    ssp_ages : 1d array
        SSP age grid in Gyr

    ssp_spec : 3d array
        SSP spectra in units of micro-Jy at a distance of 10 pc
        dimensions: (wavelengths, ages, metallicities)

    ssp_linespec : 3d array
        SSP emission line fluxes in units of ergs/s/cm2 at 10 pc
        dimensions: (wavelengths, ages, metallicities)

    sfh_ages : 1d array
        SFH age bins in Gyr

    galaxy_age : float
        maximum age of the galaxy in Gyr
        (age of Universe at the redshift of the galaxy)

    t_birth : float
        age of the birth cloud in Gyr

    Returns
    -------
    binned_ages : 1d array
        collapsed SSP age grid in Gyr
        include age of the birth cloud
        and all SFH age bins, accounting for the age of the galaxy

    binned_spec : 3d array
        new ssp_spec, accounting for the new age weighting

    binned_linespec : 3d array
        new ssp_linespec, accounting for the new age weighting

    '''
    # build the list of age bin points
    agebin_list = [0., t_birth, sfh_ages]

    # add the maximum SSP age of the galaxy and exclude all older ages
    agebin_list.append( galaxy_age )
    agebin = np.sort( np.hstack(agebin_list) )
    # exclude duplicate age points
    agebin = agebin[ np.hstack([1., np.diff(agebin)]) > 1e-10 ]
    agebin = agebin[ agebin <= galaxy_age ]

    # if there are no SSP ages between oldest SFH bin and galaxy age,
    # use the next age in the SSP grid
    if not len(ssp_ages[ (ssp_ages > agebin[-2]) & (ssp_ages <= agebin[-1]) ]):
        agebin[-1] = ssp_ages[ np.where(ssp_ages>agebin[-1])[0][0] ]

    # set up the output array
    binned_ages   = agebin[1:]
    nwave_ssp     = ssp_spec.shape[0]
    nwave_linessp = ssp_linespec.shape[0]
    nmet          = ssp_spec.shape[2]

    binned_spec     = np.zeros((nwave_ssp,     len(binned_ages), nmet))
    binned_linespec = np.zeros((nwave_linessp, len(binned_ages), nmet))

    for i in np.arange(len(binned_ages)):
        sel = np.where( (ssp_ages > agebin[i]) * (ssp_ages <= agebin[i+1]) )[0]
        wht = np.diff(np.hstack([0., 1e9 * ssp_ages[sel]]))
        wht[0] = np.diff( 1e9 * np.array([agebin[i], ssp_ages[sel][0]]) )
        for imet in np.arange(nmet):
            binned_spec[:,i,imet] = np.dot(ssp_spec[:,sel,imet],wht) / wht.sum()
            binned_linespec[:,i,imet] = np.dot(ssp_linespec[:,sel,imet],wht) / wht.sum()


    return binned_ages, binned_spec, binned_linespec


def get_coarser_wavelength_fsps(wave, spec, redwave=1e5):
    '''
    smooth the spectrum with a gaussian kernel to improve 
    computational efficiency

    only affects the wavelength grid
    (the age and metallicity grids remain unchanged)

    Parameters
    ----------
    wave : numpy array (1d)
        initial wavelength grid
    spec : numpy array
        initial SSP grid over (wave, age, metallicity)
    redwave : float
        red wavelength cutoff (in Angstroms)
    '''
# WPBWPB: coarse_wavelength, how specified by user?
    sel = np.where((wave > 500) * (wave < redwave))[0]
    spec = spec[sel, :]
    wave = wave[sel]
    G = Gaussian1DKernel(25)
    nsel = np.where(np.abs(np.diff(wave)-0.9) < 0.5)[0]
    for i in np.arange(spec.shape[1]):
        spec[nsel, i] = convolve(spec[nsel, i], G)
    ndw = 12.
    nw = np.arange(wave[nsel[0]], wave[nsel[-1]+1], ndw)
    nwb = np.hstack([nw, wave[nsel[-1]+1]])
    nwave = np.hstack([wave[:nsel[0]], nw, wave[(nsel[-1]+1):]])
    nspec = np.zeros((len(nwave), spec.shape[1]))
    for i, sp in enumerate(spec.swapaxes(0, 1)):
        hsp = (np.histogram(wave[nsel], nwb, weights=sp[nsel])[0] /
               np.histogram(wave[nsel], nwb)[0])
        nspec[:, i] = np.hstack([sp[:nsel[0]], hsp, sp[(nsel[-1]+1):]])
    return nwave, nspec


#def bin_ages_fsps_separate(args, ages, spec):
#    ''' FILL IN
#
#    Parameters
#    ----------
#    args : FILL IN
#    ages :
#        SSP age grid in Gyr
#
#    returns age (Gyr), blah
#    '''
#    weight = np.diff(np.hstack([0., 1e9 * ages]))
#
### WPBWPB delete
##    print('this is weight:')
##    print(weight)
#
#    # convert max SSP age arguments to units Gyr
#    max_ssp_ages = 10.**(np.array(args.max_ssp_age)-9.)
#
#    # start to build the list of age bin points
#    agebin_list = [0., 10.**(args.t_birth-9.)]
#
### WPBWPB delete
##    print('these are SSP ages before re-gridding')
##    print(ages)
#
#    # add the SFH age bins:
#    sfh_class = getattr(sfh, args.sfh)()
#    sfh_ages_Gyr = 10.**(np.array(sfh_class.ages)-9.)
#    agebin_list.append( sfh_ages_Gyr )
#
#    # account for the SSP ages that should not be binned:
#    agebin = np.sort( np.unique(np.hstack(agebin_list) ))
#    agebin = agebin[ agebin <= max_ssp_ages[0] ]
#    sel_unbinned = ages > max_ssp_ages[0]
#    age_unbinned, spec_unbinned = (ages[sel_unbinned], spec[:, sel_unbinned])
#
#    final_age = np.hstack([agebin[1:], age_unbinned])
#    nspec = np.zeros( (spec.shape[0], len(final_age)) ) 
#
### WPBWPB delete
##    print('this is age_unbinned:')
##    print(age_unbinned)
#
#    agebin = np.hstack([agebin, age_unbinned])
#
#    for i in np.arange(nspec.shape[1]):
#        # bin the original SSP grid, if possible
#        if True: #i < len(agebin)-1:
#            sel = np.where((ages > agebin[i]) * (ages <= agebin[i+1]))[0]
#
#            wht = np.diff(np.hstack([0., 1e9 * ages[sel]]))
#            wht[0] = np.diff( 1e9 * np.array([agebin[i], ages[sel][0]]) )
#
#            whtnew = wht.copy()
##            print('this is new wht:')
##            print(wht)
#
#            wht =  weight[sel]
## IF I apply the following line, I can combine the binned and unbinned ages
## into a single array
## But the first element is different from the original binning. hmmmmm
#            wht[0] = np.diff( 1e9 * np.array([agebin[i], ages[sel][0]]) )
#
#            if max(abs(whtnew-wht)>0):
#                print(agebin[i:i+2])
#                print('this is new wht:')
#                print(whtnew)
#
#                print('this is old wht:')
#                print(wht)
#
#            nspec[:,i] = np.dot(spec[:, sel], wht) / wht.sum()
#
#        # add in the remaining spectra from original SSP grid
#        else:
#            k = i - len(agebin) + 1
#            nspec[:,i] = spec_unbinned[:,k]
## questions:
## weight is time between SSP age and previous, or time in bin?
## combining binned spectra and unbinned SSP spectra
#
#
#    return final_age, nspec
#
#
#def bin_ages_fsps(args, ages, spec):
#    ''' FILL IN
#
#    Parameters
#    ----------
#    args : FILL IN
#    ages :
#        SSP age grid in Gyr
#
#    returns age (Gyr), blah
#    '''
#    weight = np.diff(np.hstack([0., 1e9 * ages]))
#
### WPBWPB delete
##    print('this is weight:')
##    print(weight)
#
#    # convert max SSP age arguments to units Gyr
#    max_ssp_ages = 10.**(np.array(args.max_ssp_age)-9.)
#
#    # start to build the list of age bin points
#    agebin_list = [0., 10.**(args.t_birth-9.)]
#
### WPBWPB delete
##    print('these are SSP ages before re-gridding')
##    print(ages)
#
#    # add the SFH age bins:
#    sfh_class = getattr(sfh, args.sfh)()
#    sfh_ages_Gyr = 10.**(np.array(sfh_class.ages)-9.)
#    agebin_list.append( sfh_ages_Gyr )
#
#    # account for the SSP ages that should not be binned:
#    agebin = np.sort( np.unique(np.hstack(agebin_list) ))
#    agebin = agebin[ agebin <= max_ssp_ages[0] ]
#    sel_unbinned = ages > max_ssp_ages[0]
#    age_unbinned, spec_unbinned = (ages[sel_unbinned], spec[:, sel_unbinned])
#
#    final_age = np.hstack([agebin[1:], age_unbinned])
#    nspec = np.zeros( (spec.shape[0], len(final_age)) )
#
### WPBWPB delete
##    print('this is age_unbinned:')
##    print(age_unbinned)
#
#    agebin = np.hstack([agebin, age_unbinned])
#
#    for i in np.arange(nspec.shape[1]):
#        sel = np.where((ages > agebin[i]) * (ages <= agebin[i+1]))[0]
#
#        wht = np.diff(np.hstack([0., 1e9 * ages[sel]]))
#        wht[0] = np.diff( 1e9 * np.array([agebin[i], ages[sel][0]]) )
#
#        nspec[:,i] = np.dot(spec[:, sel], wht) / wht.sum()
#
## questions:
## weight is time between SSP age and previous, or time in bin?
#
#    return final_age, nspec
#


def read_fsps_neb(filename):
    '''
    Returns
    -------
    Z : list (1 dim)
        metallicity grid in log solar units
    Age : list (1 dim)
        age grid in years (up to 10 Myr)
    logU : list (1 dim)
        ionization parameter grid
    spec : list (1 dim)
        elements are spectra (numpy array, 1 dim)
        relative line fluxes, normalized by number of ionizing photons,
        i.e., units of the table values are inverse (photons / s)
    wave : numpy array (1 dim)
        wavelength for each spectrum in Angstroms
    '''
    cnt = 0
    Z, Age, logU, spec = [], [], [], []
    with open(filename) as f:
        for lines in f:
            if cnt == 1:
                wave = np.array(lines.split(), dtype=float)
            if cnt > 1:
                l = lines.split()
                if len(l) == 3:
                    Z.append(float(l[0]))
                    Age.append(float(l[1]))
                    logU.append(float(l[2]))
                else:
                    spec.append(np.array(l, dtype=float))
            cnt += 1
    return Z, Age, logU, spec, wave


def read_fsps_file(args, metallicity):
    '''Read in the stellar population models from fsps for a given isochrone
    and metallicity.

    Parameters
    ----------
    args : class
        The args class from mcsed.parse_args()

    Returns
    -------
    ages : numpy array (1 dim)
        ages of each SPS (Gyr)
    wave : numpy array (1 dim)
        wavelength grid for each spectrum in units of Angstroms
    spec : numpy array (2 dim)
        Spectra in f_nu (micro Janskies, i.e., 1e-29 ergs/s/cm^2/Hz) at 10pc
    '''
    pc10 = 10. * 3.08567758e18
    solar_microjansky = 3.826e33 * 1e29 / (4. * np.pi * pc10**2)
    filename = op.join('SSP', 'fsps_%s_%0.4f.spec' % (args.isochrone,
                                                           metallicity))
    if not op.exists(filename):
        print('Tried to open %s' % filename)
        print('Metallicity entered, %0.4f, does not match any of the %s '
              'isochrones of the %s models' % (args.metallicity,
                                               args.isochrone, args.ssp))
        print('Metallicity options [')
        for met in args.metallicity_dict[args.ssp][args.isochrone]:
            print('%0.4f ' % met)
        print(']')
        sys.exit(1)

    # Interpret the file
    cnt = 0
    ages, spec = [], []
    with open(filename) as f:
        for lines in f:
            if cnt == 9:
                wave = np.array(lines.split(), dtype=float)
            if cnt > 9:
                l = lines.split()
                if len(l) == 4:
                    ages.append(float(l[0]))
                else:
                    spec.append(np.array(l, dtype=float))
            cnt += 1

    # convert from solar bolometric luminosity per Hz to micro-Jy at 10 pc
    spec = np.array(spec).swapaxes(0, 1) * solar_microjansky
    ages = np.array(ages)
    # exclude non-physical ages
    sel = (ages >= 6.) & (ages <= args.max_ssp_age[1]) 
    # add one additional point (if max age falls between grid points)
    if len( ages[ages > args.max_ssp_age[1]] ):
        sel_next = np.where(sel)[0][-1]+1
        sel[ sel_next ] = True

#    sel = (ages >= 6.)

    return 10**(ages[sel]-9), wave, spec[:, sel]

def get_nebular_emission(ages, wave, spec, logU, metallicity,
                         filename='nebular/ZAU_ND_pdva',
                         sollum=3.826e33, kind='both'):
    ''' 
    ages : numpy array (1 dim)
        ages of the SSP models
    wave : numpy array (1 dim)
        wavelength for SSP models
    spec : numpy array (3 dim)
        SSP spectrum for each age and each metallicity
    kind : str
        {'both', 'line', 'cont'}
        'line', 'cont' return only the line and continuum nebular emission
        'both' returns line and continuum nebular emission

    Returns
    -------
    nspec : numpy array (3 dim)
        nebular emission spectrum in units micro-Jy at 10 pc
        dimensions: (wavelengths, ages, metallicities)
    '''
    while kind not in ['line', 'cont', 'both']:
        kind = input("Invalid entry. Please enter 'line', 'cont', or 'both'")
    cont_file = filename + '.cont'
    lines_file = filename + '.lines'
    cont_res = [np.array(x) for x in read_fsps_neb(cont_file)]
    if kind != 'cont':
        lines_res = [np.array(x) for x in read_fsps_neb(lines_file)]
    # Make array of Z, age, U
    V = np.array([10**cont_res[0]*0.019, cont_res[1]/1e6,
                  cont_res[2]]).swapaxes(0, 1)
    if kind != 'line':
        # 1e48 factor avoids interpolating small numbers - removed later
        C = scint.LinearNDInterpolator(V, cont_res[3]*1e48)
    if kind != 'cont':
        L = scint.LinearNDInterpolator(V, lines_res[3]*1e48)
        garray = make_gaussian_emission(wave, lines_res[4])
# WPBWPB: uncomment
    nspec = spec * 0.
## WPBWPB delete: only when saving an array without gaussian emission:
##    print(nspec.shape)
#    nspec = nspec[0:len(lines_res[4]), :] 
#    garray = 1.
    for i, age in enumerate(ages):
        if age <= 1e-2:
            if kind != 'line':
                cont = C(metallicity, age*1e3, logU)
            if kind != 'cont':
                lines = L(metallicity, age*1e3, logU)
            # qq has units 1e-29 photons / s / cm^2 at 10 pc
            # sollum has units ergs / s
            # table values (lines) have units s / photons
            # garray has units per Hz
            # --> product has units micro-Jy at 10 pc
            qq = number_ionizing_photons(wave, spec[:, i]) / 1e48 * sollum
            if kind == 'both': 
                nspec[:, i] = (nspec[:, i] 
                               + np.interp(wave, cont_res[4], cont*qq)
                               + (garray * lines * qq).sum(axis=1))
            if kind == 'line':
## WPBWPB delete: only when saving an array without gaussian emission:
#                nspec[:, i] = lines * qq
# WPBWPB: uncomment - original code
                nspec[:, i] = (nspec[:, i] 
                               + (garray * lines * qq).sum(axis=1))
            if kind == 'cont':
                nspec[:, i] = (nspec[:, i] 
                               + np.interp(wave, cont_res[4], cont*qq))

## WPBWPB: only for case when no gaussian emission is included
#    np.savez('SSP_nongaussian', linewave0=lines_res[4], linespec0=nspec)
#    return 

    return nspec

def add_nebular_emission(ages, wave, spec, logU, metallicity,
                         filename='nebular/ZAU_ND_pdva',
                         sollum=3.826e33):
    cont_file = filename + '.cont'
    lines_file = filename + '.lines'
    cont_res = [np.array(x) for x in read_fsps_neb(cont_file)]
    lines_res = [np.array(x) for x in read_fsps_neb(lines_file)]
    # Make array of Z, age, U
    V = np.array([10**cont_res[0]*0.019, cont_res[1]/1e6,
                  cont_res[2]]).swapaxes(0, 1)
    C = scint.LinearNDInterpolator(V, cont_res[3]*1e48)
    L = scint.LinearNDInterpolator(V, lines_res[3]*1e48)
    garray = make_gaussian_emission(wave, lines_res[4])
    nspec = spec * 1.
    for i, age in enumerate(ages):
        if age <= 1e-2:
            cont = C(metallicity, age*1e3, logU)
            lines = L(metallicity, age*1e3, logU)
            qq = number_ionizing_photons(wave, spec[:, i]) / 1e48 * sollum
            nspec[:, i] = (nspec[:, i] + np.interp(wave, cont_res[4], cont*qq)
                           + (garray * lines * qq).sum(axis=1))
    return nspec

def collapse_emline_SSP(args, linewave, linespec, clight=2.99792e18):
    '''Speed up construction of emission line fluxes from the CSP

    Parameters
    ----------
    args : dictionary
        user-passed arguments from config.py and command line
    linewave : numpy array (1 dim)
        wavelength for SSP models
    linespec : numpy array (3 dim)
        SSP spectrum for each age and each metallicity

    Returns
    -------
    linewave : numpy array (1 dim)
        collapsed wavelength grid of emission lines
    linespec : numpy array (3 dim)
        collapsed SSP spectrum for each age and each metallicity
        in units of ergs / s / cm^2 at 10 pc
    ''' 

    if (not args.use_emline_flux) | (args.emline_list_dict=={}):
        return np.array([1000.,2000.]), linespec[0:2,:,:]

# loop through all emission line spectra for all ages, metallicities
# WPBWPB: generalize such that does not assume only grid over ages and metallicities, but maybe ionization parameter (or arbitrary number of properties)
# maybe I want to raise an error?
    emlines = args.emline_list_dict.keys() 
    emwaves = np.array(args.emline_list_dict.values())[:,0]
## WPBWPB delete
#    print('here are emwaves, type: %s, %s' % (emwaves, type(emwaves)))
    ssp_emline_collapsed = linespec[0:len(emlines),:,:]
## WPBWPB delete
#    waves_collapsed = []
    # loop through emission line spectra for all ages, metallicities
    dims = linespec.shape
    for i in range(dims[1]): # age array
        # if no emission lines at this age, skip metallicity grid
        if np.max(linespec[:,i,:]) <= 0:
            empty = np.zeros( (len(emlines), dims[2]) )
            ssp_emline_collapsed[:,i,:] = empty
            continue
        for j in range(dims[2]): # metallicity array
            spec = linespec[:,i,j]
            lineflux = []
            for emline in emlines:
## WPBWPB delete
#                print('here"s the emline dict and your element:')
#                print(args.emline_list_dict)
#                print([type(val) for val in args.emline_list_dict.values()])
#                print(emline)
#                print(args.emline_list_dict[emline])
#                print(type(args.emline_list_dict[emline]))

                w = args.emline_list_dict[emline][0]
                indx = np.searchsorted(linewave, w)
                if spec[indx] <= 0:
                    lineflux.append(0.)
                    continue
                # find indices when target line goes to zero
                indx_zero = np.where(spec==0)[0]
                indx_lo = indx - np.min( indx - indx_zero[indx_zero<indx])
                indx_hi = indx + np.min( indx_zero[indx_zero>indx] - indx)
                dnu = np.diff( clight / linewave[indx_lo:indx_hi+2])
                y = spec[indx_lo:indx_hi+1]
                # convert flux density in micro-Jy at 10 pc 
                # to total flux (ergs/s/cm2) at 10pc
                lineflux.append(np.dot( y, np.abs(dnu) ) / 1e29)
            ssp_emline_collapsed[:,i,j] = lineflux

    # return modified emission line wavelength and SSP grid: 
    # fluxes at discrete wavelengths
    return np.array(emwaves), ssp_emline_collapsed


def make_gaussian_emission(wavebig, wave, stddev=1., clight=2.99792e18):
    ''' 

    Parameters
    ----------
    wavebig, wave both in units of Angstroms

    Returns
    -------
    gspec : FILL IN
        in units per Hz
    '''
    gspec = np.zeros((len(wavebig), len(wave)))
    G = Gaussian1DKernel(stddev).array
    mid = len(G) / 2
    dw = np.diff(wavebig)
    for i, w in enumerate(wave):
        xl = np.argmin(np.abs(wavebig - w))
        if (xl > mid) and ((xl+mid) < len(wavebig)):
            gspec[(xl-mid):(xl+mid+1), i] = G / clight * w**2 / dw[xl]
    return gspec


def number_ionizing_photons(wave, spectrum, clight=2.99792e18,
                            hplanck=6.626e-27):
    '''

    Parameters
    ----------
    wave : numpy array (1 dim)
        wavelength grid in units of Angstroms
    spectrum : numpy array (1 dim)
        Spectrum in f_nu (micro Janskies, i.e., 1e-29 ergs/s/cm^2/Hz) at 10pc
 
    Returns
    -------
    float
        number of photons capable of ionizing Hydrogen
        in units of 1e-29 photons / s / cm^2 at 10 pc
    '''
    nu = clight / wave 
    xlim = np.searchsorted(wave, 912., side='right')
    x1 = np.abs(np.diff(nu[:xlim])) 
    x2 = spectrum[:xlim] / nu[:xlim] 
    x3 = (x2[:-1] + x2[1:]) / 2. 
    return np.sum(x1 * x3) / hplanck


def read_ssp_fsps(args):
    ''' Read in SPS model and return ages, wavelength, and spectra

    Parameters
    ----------
    args : class
        The args class from mcsed.parse_args()

    Returns
    -------
    ages : numpy array (1d)
        SSP age grid (Gyr)
    wave : numpy array (1d)
        SSP wavelength grid (Angstroms)
    spec : numpy array (3d)
        SSP spectra in units micro-Jy at a distance of 10 pc
        dimensions: (wave, ages, metallicities)
    metallicities : numpy array (1d)
        SSP metallicities (in values of Z, where Zsolar = 0.019)
    linewave : numpy array (1d)
        rest-frame wavelengths of emission-line fluxes 
        (if used in model calculation)
    linespec : numpy array (3d)
        line fluxes in units ergs / cm2 / s at distance of 10 pc
        dimensions: (linewave, ages, metallicities)

    '''
    metallicities = np.array(args.metallicity_dict[args.ssp][args.isochrone])

## WPBWPB: m list unused
#    s, ls, m = ([], [], [])
    s, ls = [], []
    for met in metallicities:
        if args.ssp.lower() == 'fsps':
            ages, wave, spec = read_fsps_file(args, met)
        # add new SSP subroutines here:
# WPBWPB: only carry linespec if going to measure emission lines?
        linespec = get_nebular_emission(ages, wave, spec, args.logU,
                                        met, kind='line')
        spec = add_nebular_emission(ages, wave, spec, args.logU,
                                        met)

## WPBWPB add comment
## WPBWPB testing if age binning is issue in binned SFH for the too-old bins
##        if args.sfh in ['binned_lsfr', 'binned_fmass']: 
#        if False:
##            print('YOU NEED TO DEAL WITH THE SSP AGE BINNING!!!!')
#            ages0 = ages.copy()
### WPBWPB delete: useful if comparing to old "bin ages" function
##            ages0 = np.log10(ages0)+9.
#            ages, spec = bin_ages_fsps(args, ages0, spec)
#            ages9, linespec = bin_ages_fsps(args, ages0, linespec)

        # do not smooth the emission line grid
        wave0 = wave.copy()
        if args.fit_dust_em:
            wave, spec = get_coarser_wavelength_fsps(wave0, spec, redwave=350e4)
        else:
            wave, spec = get_coarser_wavelength_fsps(wave0, spec)
## WPBWPB delete -- appears unused
#        wei = (np.diff(np.hstack([0., ages])) *
#               getattr(sfh, args.sfh)().evaluate(ages))
        s.append(spec)
        ls.append(linespec)
## WPBWPB following line appears unused... 
#        m.append(np.dot(spec, wei)/spec.shape[1])

    # save plot of SSP spectra
    if args.output_dict['template spec']:
        if args.metallicity:
            imet = np.argmin(abs(args.metallicity - metallicities))
        else: # show SSP grid for 40% solar, if no metallicity is set
            imet = np.argmin(abs(0.0077 - metallicities))
        fig = plt.figure(figsize=(8, 8))
        import seaborn as sns
        colors = sns.color_palette("coolwarm", s[imet].shape[1])
        wei = np.diff(np.hstack([0., ages]))
#        wei = np.ones(s[imet].shape[1])
        for i in np.arange(s[imet].shape[1]):
            plt.plot(wave, s[imet][:, i] * wei[i] / 1e8, color=colors[i])
        plt.xlim([900., 40000.])
        plt.xscale('log')
        plt.yscale('log')
        plt.ylim([1e-5, 20])
#        plt.ylim([1e-5, 10.**3.5])
        plt.xlabel('Wavelength [$\\rm{\AA}$]')
        plt.ylabel('Relative $f_\\nu$')
        plt.savefig('template_spectra_plot.%s' % args.output_dict['image format'])
        plt.close(fig)

    spec = np.moveaxis(np.array(s), 0, 2)
    linespec = np.moveaxis(np.array(ls), 0, 2)

## WPBWPB delete
#    linespec0 = linespec.copy()

    # Collapse the emission line SSP grid
    linewave, linespec = collapse_emline_SSP(args, wave0, linespec) 

## WPBWPB delete
#    np.savez('SSP', wave=wave, spec=spec, linewave=linewave, linespec=linespec, linewave0=wave0, linespec0=linespec0)

## WPB delete
#    print('these are ages: %s' % ages)

    return ages, wave, spec, metallicities, linewave,linespec


