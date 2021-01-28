""" script for running MCSED

.. moduleauthor:: Greg Zeimann <gregz@astro.as.utexas.edu>

"""

from __future__ import absolute_import
import sys
import argparse as ap
import numpy as np
import os.path as op
import logging
import config
import ism_igm
from ssp import read_ssp_fsps, bin_ssp_ages
from astropy.io import fits
from astropy.table import Table, vstack
from mcsed import Mcsed
from distutils.dir_util import mkpath
from cosmology import Cosmology

sys.path.insert(0,'3dhst_catalogs')
import filter_info
sys.path.insert(0, 'SSP')
import ssp_metallicity_info

def setup_logging():
    '''Setup Logging for MCSED, which allows us to track status of calls and
    when errors/warnings occur.

    Returns
    -------
    log : class
        log.info() is for general print and log.error() is for raise cases
    '''
    log = logging.getLogger('mcsed')
    if not len(log.handlers):
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
        log = logging.getLogger('mcsed')
        log.setLevel(logging.DEBUG)
        log.addHandler(handler)
    return log


def str2bool(v, log):
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        log.warning('Could not interpret "metallicity" argument, by '
                    'default it will be set to False')
        return False


def parse_args(argv=None):
    '''Parse arguments from commandline or a manually passed list

    Parameters
    ----------
    argv : list
        list of strings such as ['-f', 'input_file.txt', '-s', 'default.ssp']

    Returns
    -------
    args : class
        args class has attributes of each input, i.e., args.filename
        as well as astributes from the config file
    '''

    # Avoid an infinite loop between parallel and series functions
    already_parallel = False
    if type(argv) is list:
        if '--already_parallel' in argv:
            argv.remove('--already_parallel')
            already_parallel = True

    # Pass "count" keyword (indexing objects in test mode) 
    count = 1
    if type(argv) is list:
        if '--count' in argv:
            indx = argv.index('--count')
            count = int( argv[indx+1] )
            del argv[indx: indx+2]

    parser = ap.ArgumentParser(description="MCSED",
                               formatter_class=ap.RawTextHelpFormatter)

    parser.add_argument("-f", "--filename",
                        help='''File to be read for galaxy data''',
                        type=str, default=None)

    parser.add_argument("-o", "--output_filename",
                        help='''Output filename for given run''',
                        type=str, default='test.dat')

    parser.add_argument("-p", "--parallel",
                        help='''Running in parallel?''',
                        action="count", default=0)

    parser.add_argument("-t", "--test",
                        help='''Test mode with mock galaxies''',
                        action="count", default=0)

    parser.add_argument("-tf", "--test_field",
                        help='''Test filters will match the given field''',
                        type=str, default='cosmos')

    parser.add_argument("-no", "--nobjects",
                        help='''Number of test objects''',
                        type=int, default=None)

    parser.add_argument("-s", "--ssp",
                        help='''SSP Models, default fsps''',
                        type=str, default=None)

    parser.add_argument("-i", "--isochrone",
                        help='''Isochrone for SSP model, e.g. padova''',
                        type=str, default=None)

    parser.add_argument("-sfh", "--sfh",
                        help='''Star formation history, e.g. constant''',
                        type=str, default=None)

    parser.add_argument("-dl", "--dust_law",
                        help='''Dust law, e.g. calzetti''',
                        type=str, default=None)

    parser.add_argument("-de", "--dust_em",
                        help='''Dust emission class, e.g., DL07 (Draine & Li 2007)\n''' 
                            +'''or False (if dust emission should be ignored)''',
                        type=str, default=None)

    parser.add_argument("-aeb", "--assume_energy_balance",
                        help='''If selected, normalization of dust IR emission based'''
                            +'''on attenuation amount''',
                        action="count", default=0)

    parser.add_argument("-z", "--metallicity",
                        help='''Fixed metallicity for SSP models (0.019 is solar),\n'''
                            +'''or False if stellar metallicity is a free parameter''',
                        type=str, default=None)

    parser.add_argument("-nw", "--nwalkers",
                        help='''Number of walkers for EMCEE''',
                        type=int, default=None)

    parser.add_argument("-ns", "--nsteps",
                        help='''Number of steps for EMCEE''',
                        type=int, default=None)

    parser.add_argument("-lu", "--logU",
                        help='''Ionization Parameter for nebular gas''',
                        type=float, default=None)

    parser.add_argument("-ism", "--ISM_correct_coords",
                        help='''If a coordinate system is given, MW dust correction will
                        be performed; default None''',
                        type=str, default=None)

    parser.add_argument("-igm", "--IGM_correct",
                        help='''If selected, Madau statistical IGM correction will be done
                        (affecting wavelengths up to rest-frame Ly-alpha)''',
                        action="count", default=0)

    # Initialize arguments and log
    args = parser.parse_args(args=argv)
    args.log = setup_logging()

    # Use config values if none are set in the input
    arg_inputs = ['ssp', 'metallicity', 'isochrone', 'sfh', 'dust_law',
                  't_birth', 'nwalkers', 'nsteps', 'logU', 
                  'phot_floor_error', 'emline_floor_error', 'absindx_floor_error',  
                  'model_floor_error', 'nobjects', 'test_zrange', 'blue_wave_cutoff', 
                  'dust_em', 'Rv', 'EBV_old_young', 'wave_dust_em',
                  'emline_list_dict', 'emline_factor', 'use_input_data',
                  'absorption_index_dict', 'separate_stars_gas',
                  'output_dict', 'param_percentiles', 'reserved_cores', 
                  'assume_energy_balance', 'ISM_correct_coords', 'IGM_correct']
    for arg_i in arg_inputs:
        try:
            if getattr(args, arg_i) in [None, 0]:
                setattr(args, arg_i, getattr(config, arg_i))
        except AttributeError:
            setattr(args, arg_i, getattr(config, arg_i))

    # Read the filter information
    filterarg_inputs = ['filt_dict', 'catalog_filter_dict', 'catalog_maglim_dict']
    for arg_i in filterarg_inputs:
        setattr(args, arg_i, getattr(filter_info, arg_i))

    # Read the SSP metallicity information
    setattr(args, 'metallicity_dict', getattr(ssp_metallicity_info, 'metallicity_dict'))

    # If a test field is specified on the command line, initialize test mode
    if '-tf' in argv:
        args.test = True

    # If coords is not None, ISM correction will be applied
    if args.ISM_correct_coords is not None: 
        args.ISM_correct = True
    else:
        args.ISM_correct = False

    # Ignore ISM/IGM corrections in test mode
    if args.test:
        args.ISM_correct = False
        args.IGM_correct = False

    # Set the maximum SSP age (speeds calculation)
    args.max_ssp_age = get_max_ssp_age(args)

    # Set metallicity as free or fixed parameter
    try:
        if args.metallicity not in ['0','1']:
            args.metallicity = float(args.metallicity)
        else:
            if args.metallicity=='0': 
                args.metallicity = False
            else:
                args.metallicity = True
    except ValueError:
        args.metallicity = str2bool(str(args.metallicity),args.log)
        if args.metallicity:
            args.log.info("Fixing metallicity at Z = 0.0077 (Zsolar = 0.019)")
            args.metallicity = 0.0077

    # Avoid an infinite loop between parallel and series functions
    if already_parallel:
        args.already_parallel = True
    else:
        args.already_parallel = False

    # Pass "count" keyword (indexing objects in test mode) 
    args.count = count

    # Determine whether emission lines / absorption line indices are used
    # to constrain the models (not included in test mode)
    if (not args.test) & (not args.use_input_data):
        args.use_emline_flux = False
        args.use_absorption_indx = False
    else:
        args.use_emline_flux = True
        args.use_absorption_indx = True


    if (type(args.emline_list_dict)!=dict) | (args.test):
        args.emline_list_dict={}
    if (type(args.absorption_index_dict)!=dict) | (args.test):
        args.absorption_index_dict={}

    # Set up dust emission arguments
    if isinstance(args.dust_em, str):
        args.fit_dust_em = True
    elif args.dust_em==True:
        args.dust_em = 'DL07'
        args.fit_dust_em = True
    else:
        args.dust_em = 'DL07'
        args.fit_dust_em = False

    return args


def build_filter_matrix(args, wave):
    '''Build a filter matrix with each row being an index of wave and
    each column being a unique filter.  This makes computation from spectra
    to magnitudes quick and easy.

    Parameters
    ----------
    args : class
        The args class is carried from function to function with information
        from command line input and config.py
    wave : numpy array
        The wave array corresponds to the wavelengths of the SSP models being
        used.

    Returns
    -------
    Fil_matrix : numpy array (2 dim)
        As mentioned above, the Fil_matrix has rows of wavelength and
        columns for each filter in args.filt_dict/config.filt_dict
    '''
    nfilters = len(args.filt_dict)
    Fil_matrix = np.zeros((len(wave), nfilters))
    for i in np.arange(nfilters):
        wv, through = np.loadtxt(op.join('FILTERS', args.filt_dict[i]),
                                 unpack=True)
        new_through = np.interp(wave, wv, through, 0.0, 0.0)
        S = np.sum(new_through)
        if S == 0.:
            S = 1.
        Fil_matrix[:, i] = new_through / S

    return Fil_matrix


def get_test_filters(args):
    '''Used in test mode, this function loops through args.filt_dict and sets
    a flag to true if the filter is in args.test_filter_dict or false if it
    is not.  This filter_flag is used later in the quick calculation of
    filter magnitudes.

    Parameters
    ----------
    args : class
        The args class is carried from function to function with information
        from command line input and config.py

    Returns
    -------
    filter_flag : numpy array (bool)
        Explained above.
    '''
    nfilters = len(args.filt_dict)
    filter_flag = np.zeros((nfilters,), dtype=bool)
    for i in args.filt_dict.keys():
        if i in args.catalog_filter_dict[args.test_field]:
            filter_flag[i] = True
    return filter_flag


def get_maglim_filters(args):
    '''Used in test mode, this function loops through args.filt_dict and
    retrieves the 5-sigma magnitude limit for each filter (AB), and returns
    the appropriate microjansky 1-sigma error.

    Parameters
    ----------
    args : class
        The args class is carried from function to function with information
        from command line input and config.py

    Returns
    -------
    photerror : numpy array (float)
        Explained above.
    '''
    nfilters = len(args.filt_dict)
    photerror = np.zeros((nfilters,), dtype=float)
    for i in args.filt_dict.keys():
        if i in args.catalog_filter_dict[args.test_field]:
            maglim = args.catalog_maglim_dict[args.test_field][i]
            photerror[i] = 10**(-0.4 * (maglim - 23.9)) / 5.
    return photerror


def get_max_ssp_age(args, z=None):
    '''
    Identify the maximum SSP age for the sample (i.e., at upper/lower redshifts)

    Max SSP age is the time between redshift z=20 and redshift of the galaxy
 
    Parameters
    ----------
    args : class
        The args class is carried from function to function with information
        from command line input and config.py

    z : float (optional)
        if specific redshift is passed, evaluate max age at that redshift
        otherwise, evaluate for the range of redshifts of the sample

    Returns
    -------
    maxage : tuple of (float, float)
        the youngest and oldest maximum age (in log years) of sample galaxies
    '''
    C = Cosmology()
    if z is not None:
        maxage = np.log10(C.lookback_time(20)-C.lookback_time(z)) + 9.
        return maxage

    if not args.test:
        F = Table.read(args.filename, format='ascii')
        z = F['z']
        zrange = (min(z), max(z))
    else:
        zrange = args.test_zrange

    # ages in log years:
    maxage_lo = np.log10(C.lookback_time(20)-C.lookback_time(zrange[1])) + 9.
    maxage_hi = np.log10(C.lookback_time(20)-C.lookback_time(zrange[0])) + 9.
    return (maxage_lo, maxage_hi)


def read_input_file(args):
    '''This function reads a very specific input file and joins it with
    archived 3dhst catalogs.  The input file should have the following columns:
    Field, ID, z

    Parameters
    ----------
    args : class
        The args class is carried from function to function with information
        from command line input and config.py

    Returns
    -------
    y : numpy array (2 dim)
        Photometric flux densities (in units of micro-Janskies)
    yerr : numpy array (2 dim)
        Photometric errors
    z : numpy array (1 dim)
        Redshift from the file returned as a numpy array
    flag : numpy array (2 dim)
        Flag set to True for filters in the catalog_filter_dict in config.py
    em : Astropy Table (2 dim)
        Emission line fluxes in ergs / cm2 / s
        Desired lines are read from dictionary in config.py
    emerr : Astropy Table (2 dim)
        Emission line errors in ergs / cm2 / s
    absindx : Astropy Table (2 dim)
        Absorption line indices read from the input file
    absindx_e : Astropy Table (2 dim)
        Errors on the absorption line indices
    '''
    F = Table.read(args.filename, format='ascii')
    # keep track of which columns from the input file are utilized
    Fcols = F.colnames
    nobj = len(F['Field'])

    # redshift array
    z = F['z']
    Fcols.remove('z')

    # Skelton catalogs
    fields = ['aegis', 'cosmos', 'goodsn', 'goodss', 'uds']
    name_base = '_3dhst.v4.1.cat.FITS'
    field_dict = {}
    for field in fields:
        field_dict[field] = fits.open(op.join('3dhst_catalogs',
                                              field+name_base))[1]

    # check whether any additional photometry is provided by the user
    input_filters = [col.strip('f_') for col in Fcols if (len(col)>1) & (col[0:2]=='f_')]
    infilt_dict = {}
    if args.use_input_data:
        for fname in input_filters:
            if op.exists('FILTERS/%s.res' % fname):
                if '%s.res' % fname not in args.filt_dict.values():
                    findex = max(args.filt_dict.keys())+1
                else:
                    findex = args.filt_dict.keys()[args.filt_dict.values().index('%s.res' % fname)]
                infilt_dict[ findex ] = '%s.res' % fname
                Fcols = [c for c in Fcols if c not in ['f_'+fname, 'e_'+fname]]
                args.filt_dict.update(infilt_dict)
            else:
                args.log.info('*CAUTION* %s.res filter curve does not exist:' % fname)

    nfilters = len(args.filt_dict)
    y = np.zeros((nobj, nfilters))
    yerr = np.zeros((nobj, nfilters))
    flag = np.zeros((nobj, nfilters), dtype=bool)

    # convert from mag_zp = 25 to microjanskies (mag_zp = 23.9)
    fac = 10**(-0.4*(25.0-23.9))
    phot_fill_value = -99 # null value, should not be changed

    # assemble photometry
    for i, datum in enumerate(F):
        loc = datum['Field'].lower()

        for j, ind in enumerate(args.filt_dict.keys()):
            if loc in args.catalog_filter_dict.keys():
                if ind in args.catalog_filter_dict[loc].keys():
                    colname  = "f_"+args.catalog_filter_dict[loc][ind]
                    ecolname = "e_"+args.catalog_filter_dict[loc][ind]
                elif ind in infilt_dict.keys():
                    colname  = "f_"+infilt_dict[ind].split('.res')[0]
                    ecolname = "e_"+infilt_dict[ind].split('.res')[0]
                else:
                    y[i, j] = 0.0
                    yerr[i, j] = 0.0
                    flag[i, j] = False
                    continue
            else:
                if ind in infilt_dict.keys():
                    colname  = "f_"+infilt_dict[ind].split('.res')[0]
                    ecolname = "e_"+infilt_dict[ind].split('.res')[0]
                else:
                    y[i, j] = 0.0
                    yerr[i, j] = 0.0
                    flag[i, j] = False
                    continue
            if loc in field_dict.keys():
                if colname in field_dict[loc].columns.names:
                    fi  = field_dict[loc].data[colname][int(datum['ID'])-1]
                    fie = field_dict[loc].data[ecolname][int(datum['ID'])-1]
                elif colname in F.colnames:
                    fi  = datum[colname]
                    fie = datum[ecolname]
                else:
                    y[i, j] = 0.0
                    yerr[i, j] = 0.0
                    flag[i, j] = False
                    continue
            else:
                if colname in F.colnames:
                    fi  = datum[colname]
                    fie = datum[ecolname]
                else:
                    y[i, j] = 0.0
                    yerr[i, j] = 0.0
                    flag[i, j] = False
                    continue
            if (fi > phot_fill_value):
                y[i, j] = fi*fac
                flag[i, j] = True
                # use a floor error if necessary
                if fi != 0:
                    yerr[i, j] = np.abs(np.max([args.phot_floor_error,
                                        np.abs(fie/fi)]) * fi * fac)
                else:
                    yerr[i, j] = 0.0
                    flag[i, j] = False
            else:
                y[i, j] = 0.0
                yerr[i, j] = 0.0
                flag[i, j] = False


    # read in emission line fluxes, if provided
    line_fill_value = -99 # null value, should not be changed
    if args.use_emline_flux:
        em, emerr = Table(), Table()

        for emline in list(args.emline_list_dict.keys()):
            colname, ecolname = '%s_FLUX' % emline, '%s_ERR' % emline
            if colname in Fcols:
                em_arr = np.array(F[colname]  * args.emline_factor)
                emerr_arr = np.max([abs(F[ecolname]),
                                    args.emline_floor_error*abs(np.array(F[colname]))],0)
                emerr_arr *= args.emline_factor

                # account for objects with null measurements
                null_idx = np.where(abs(np.array(F[colname])-line_fill_value)<1e-10)[0]
                em_arr[null_idx] = line_fill_value
                emerr_arr[null_idx] = line_fill_value

                em[colname]     = em_arr
                emerr[ecolname] = emerr_arr

                Fcols = [c for c in Fcols if c not in [colname, ecolname]]
            else:
                del args.emline_list_dict[emline]
        if not args.emline_list_dict:
            em    = np.full((len(F),2), line_fill_value)
            emerr = np.full((len(F),2), line_fill_value)
    else:
        em    = np.full((len(F),2), line_fill_value)
        emerr = np.full((len(F),2), line_fill_value)

    # read in absorption line indices, if provided
    if args.use_absorption_indx:
        absindx, absindx_e = Table(), Table()
        for indx in list(args.absorption_index_dict.keys()):
            colname, ecolname = '%s_INDX' % indx, '%s_Err' % indx
            # note the index units (for applying the floor error)
            unit = args.absorption_index_dict[indx][-1]
            if colname in Fcols:
                indx_arr = np.array(F[colname])
                if unit == 1: # magnitudes
                    efloor = 2.5*np.log10(1.+args.absindx_floor_error)
                    efloor_arr = np.array([efloor]*len(F))
                else:
                    efloor_arr = args.absindx_floor_error*abs(np.array(F[colname]))
                indxerr_arr = np.max([abs(F[ecolname]), efloor_arr],0)

                # account for objects with null measurements
                null_idx = np.where(abs(np.array(F[colname])-line_fill_value)<1e-10)[0]
                indx_arr[null_idx] = line_fill_value
                indxerr_arr[null_idx] = line_fill_value

                absindx[colname]    = indx_arr
                absindx_e[ecolname] = indxerr_arr

                Fcols = [c for c in Fcols if c not in [colname, ecolname]]
            else:
                del args.absorption_index_dict[indx]
        if not args.absorption_index_dict:
            absindx   = np.full((len(F),2), line_fill_value)
            absindx_e = np.full((len(F),2), line_fill_value)
    else:
        absindx   = np.full((len(F),2), line_fill_value)
        absindx_e = np.full((len(F),2), line_fill_value)


    # warn of any unused columns from the input file
    Fcols = [c for c in Fcols if c not in ['Field', 'ID']]
    if (Fcols!=[]) & (args.use_input_data):
        Fcols_str = '['
        for c in Fcols:
            if c!=Fcols[-1]:
                Fcols_str+= c+', '
            else:
                Fcols_str+= c+']'
        args.log.info('*CAUTION* unread columns in the input file: '+Fcols_str)

    return y, yerr, z, flag, F['ID'], F['Field'], em, emerr, absindx, absindx_e


def draw_uniform_dist(nsamples, start, end):
    ''' Draw random samples from a uniform distribution

    Parameters
    ----------
    nsamples : int
        Number of draws
    start : float
        lower bound
    end : float
        higher bound

    Returns
    -------
    uniform_sample : numpy array (1 dim)
        randomly drawn variables from a uniform distribution
    '''
    return np.random.rand(nsamples)*(end-start) + start


def draw_gaussian_dist(nsamples, means, sigmas):
    ''' Draw random samples from a normal distribution

    Parameters
    ----------
    nsamples : int
        Number of draws
    means : numpy array or list
        Average values to draw from
    sigmas : numpy array or list
        Standard deviation of the return distributions

    Returns
    -------
    normal_sample : numpy array (1 dim)
        randomly drawn variables from a normal distribution
    '''
    m = len(means)
    N = np.random.randn(nsamples * m).reshape(nsamples, m)
    return sigmas * N + means


def mock_data(args, mcsed_model, nsamples=5, phot_error=0.05):
    ''' Create mock data to test quality of MCSED fits

    Parameters
    ----------
    args : class
        The args class is carried from function to function with information
        from command line input and config.py
    mcsed_model : class
        Mcsed class for building fake galaxies given input thetas

    Returns
    -------
    y : numpy array (2 dim)
        Photometric magnitudes for mock galaxies
    yerr : numpy array (2 dim)
        Photometric errors in magnitudes
    z : numpy array (1 dim)
        Redshift for mock galaxies
    truth : numpy array (2 dim)
        Mock input parameters for each fake galaxy, e.g. dust, sfh, mass
    '''
    np.random.seed()
    thetas = mcsed_model.get_init_walker_values(num=nsamples, kind='ball')
    zmin, zmax = args.test_zrange
    zobs = draw_uniform_dist(nsamples, zmin, zmax)
    params, y, yerr, true_y = [], [], [], []

    for theta, z in zip(thetas, zobs):
        mcsed_model.set_class_parameters(theta)
        mcsed_model.set_new_redshift(z)
        if mcsed_model.dust_em_class.assume_energy_balance:
            mcsed_model.spectrum, mass, mdust_eb = mcsed_model.build_csp()
        else:
            mcsed_model.spectrum, mass = mcsed_model.build_csp()
            mdust_eb = None
        sfr10,sfr100,fpdr = mcsed_model.get_derived_params()

        f_nu = mcsed_model.get_filter_fluxdensities()
        if args.test_field in args.catalog_maglim_dict.keys():
            f_nu_e = get_maglim_filters(args)[mcsed_model.filter_flag]
            f_nu_e = np.max([f_nu_e, f_nu * phot_error], axis=0)
        else:
            f_nu_e = f_nu * phot_error
        y.append(f_nu_e*np.random.randn(len(f_nu)) + f_nu)
        yerr.append(f_nu_e)
        true_y.append(f_nu)
        derived_param_list = [np.log10(mass)]
        for par in [sfr10, sfr100, fpdr, mdust_eb]:
            if par is not None:
                derived_param_list.append( np.log10(par) )
        params.append(list(theta) + derived_param_list)

    return y, yerr, zobs, params, true_y


def main(argv=None, ssp_info=None):
    '''
    Execute the main functionality of MCSED

    Test mode: "python run_mcsed_fit.py -t"

    Live mode: "python run_mcsed_fit.py -f test_data.dat"

    For a "live" run, the key input ("-f") is a file with three columns:
    Field ID z

    If using the Skelton catalog:
        The field options are: cosmos, goodsn, goodss, aegis, uds
        The ID is the skelton photometric id for the given field
        The redshift, z, is fixed in the fitting
    '''

    # Make output folder if it doesn't exist
    mkpath('output')

    # Get Inputs
    if argv == None:
        argv = sys.argv
        argv.remove('run_mcsed_fit.py')

    args = parse_args(argv)

    # Catch to run in parallel
    if (args.parallel) & (not args.already_parallel):
        import run_mcsed_parallel
        run_mcsed_parallel.main_parallel(argv=argv)
        return   

    # Load Single Stellar Population model(s)
    if ssp_info is None:
        args.log.info('Reading in SSP model')
        ages, wave, starSSP, nebSSP, met, emlinewave, emlinefluxSSP = read_ssp_fsps(args)
    else:
        ages, wave, starSSP, nebSSP, met, emlinewave, emlinefluxSSP = read_ssp_fsps(args)

    # If not returning best-fit stellar/nebular spectra separately, combine them
    if not args.separate_stars_gas:
        starSSP += nebSSP
        nebSSP = None

    # Read in input data, if not in test mode 
    if not args.test: 
        input_file_data = read_input_file(args) 
    else:
        input_file_data = None

    # Get ISM and/or ISM correction
    if args.IGM_correct:
        tauIGMf = ism_igm.get_tauIGMf()
    if args.ISM_correct:
        tauISMf = ism_igm.get_tauISMf()

    # Build Filter Matrix
    filter_matrix = build_filter_matrix(args, wave)

    # Make one instance of Mcsed for speed on initialization
    # (relevant variables are reassigned for each galaxy)
    mcsed_model = Mcsed(filter_matrix, wave, ages, met, starSSP, nebSSP, 
                        emlinewave, emlinefluxSSP,
                        args.sfh, args.dust_law, args.dust_em, 
                        nwalkers=args.nwalkers, nsteps=args.nsteps,
                        sigma_m=args.model_floor_error)

    # Communicate emission line measurement preferences
    mcsed_model.use_emline_flux = args.use_emline_flux
    mcsed_model.emline_dict = args.emline_list_dict
    mcsed_model.use_absorption_indx = args.use_absorption_indx
    mcsed_model.absindx_dict = args.absorption_index_dict

    # Adjust Rv in the dust attenuation model, if specified in config file
    # (otherwise, use the default value for the requested dust law)
    if args.Rv >= 0:
        mcsed_model.dust_abs_class.Rv = args.Rv
    else:
        args.Rv = mcsed_model.dust_abs_class.Rv

    # Adjust the relative attenuation between young/old populations in the dust model
    # E(B-V)_diffuse = EBV_old_young * E(B-V)_birthcloud
    mcsed_model.dust_abs_class.EBV_old_young = args.EBV_old_young

    # Specify the age of the birth cloud (suffer different attenuation)
    mcsed_model.t_birth = 10**(args.t_birth-9.) # Gyr

    # Specify whether metallicity is fixed 
    if args.metallicity:
        mcsed_model.met_class.fix_met = True
        Zsolar = 0.019
        mcsed_model.met_class.met = np.log10(args.metallicity/Zsolar)
    else:
        mcsed_model.met_class.fix_met = False

    # Specify whether dust emission is fixed
    if (not args.fit_dust_em) | (args.test):
        mcsed_model.dust_em_class.fixed = True
    else:
        mcsed_model.dust_em_class.fixed = False

    # Specify whether energy balance is assumed
    if args.assume_energy_balance:
        if args.fit_dust_em:
            mcsed_model.dust_em_class.assume_energy_balance = True
        else:
            mcsed_model.dust_em_class.assume_energy_balance = False
    else:
        mcsed_model.dust_em_class.assume_energy_balance = False

    # Build names for parameters and labels for table
    names = mcsed_model.get_param_names()
    names.append('Log Mass')
    names.append('SFR10')
    names.append('SFR100')
    if not mcsed_model.dust_em_class.fixed:
        names.append('fPDR')
    if mcsed_model.dust_em_class.assume_energy_balance:
        names.append("Mdust_EB")

    percentiles = args.param_percentiles 
    labels = ['Field', 'ID', 'z']
    for name in names:
        labels = labels + [name + '_%02d' % per for per in percentiles]
    formats = {}
    for label in labels:
        formats[label] = '%0.3f'

    # If test mode, add truth values for table labels
    if args.test:
        for name in names:
            labels.append(name + '_truth')
            formats[labels[-1]] = '%0.3f'
    formats['Field'], formats['ID'] = ('%s', '%05d')

    mcsed_model.table = Table(names=labels, dtype=['S10', 'i4'] +
                              ['f8']*(len(labels)-2))

    # MAIN FUNCTIONALITY
    if args.test:
        fl = get_test_filters(args)
        mcsed_model.filter_flag = fl * True
        default = mcsed_model.get_params()
        y, yerr, z, truth, true_y = mock_data(args, mcsed_model,
                                              phot_error=args.phot_floor_error,
                                              nsamples=args.nobjects)

        cnts = np.arange(args.count, args.count + len(z))

        for yi, ye, zi, tr, ty, cnt in zip(y, yerr, z, truth, true_y, cnts):
            mcsed_model.input_params = tr
            mcsed_model.filter_flag = fl * True
            mcsed_model.set_class_parameters(default)
            mcsed_model.data_fnu = yi
            mcsed_model.data_fnu_e = ye
            mcsed_model.true_fnu = ty
            mcsed_model.set_new_redshift(zi)
            mcsed_model.data_emline = [-99]
            mcsed_model.data_emline_e = [-99]
            mcsed_model.data_absindx = [-99]
            mcsed_model.data_absindx_e = [-99]

            # Remove filters containing Lyman-alpha (and those blueward)
            mcsed_model.remove_waverange_filters(0., args.blue_wave_cutoff, restframe=True)
            # Remove filters dominated by dust emission, if applicable
            if not args.fit_dust_em:
                mcsed_model.remove_waverange_filters(args.wave_dust_em*1e4,1e10,
                                                     restframe=True)

            mcsed_model.fit_model()
            mcsed_model.set_median_fit()
            if args.output_dict['sample plot']:
                mcsed_model.sample_plot('output/sample_fake_%05d_%s_%s' % 
                                        (cnt, args.sfh, args.dust_law),
                                        imgtype = args.output_dict['image format'])

            if args.output_dict['triangle plot']:
                mcsed_model.triangle_plot('output/triangle_fake_%05d_%s_%s' % 
                                          (cnt, args.sfh, args.dust_law),
                                          imgtype = args.output_dict['image format'])
            mcsed_model.table.add_row(['Test', cnt, zi] + [0.]*(len(labels)-3))

            last = mcsed_model.add_fitinfo_to_table(percentiles)
            mcsed_model.add_truth_to_table(tr, last)
            print(mcsed_model.table)

            if names[-1] != 'Ln Prob':
                names.append('Ln Prob')
            if args.output_dict['fitposterior']:
                T = Table(mcsed_model.samples, names=names)
                T.write('output/fitposterior_fake_%05d_%s_%s.dat' % (cnt, args.sfh, args.dust_law),
                        overwrite=True, format='ascii.fixed_width_two_line')
            if args.output_dict['bestfitspec']:
                bestfitspec_data = [mcsed_model.wave, mcsed_model.medianspec, mcsed_model.true_spectrum]
                bestfitspec_name = ['wavelength', 'spectrum', 'true_spectrum']
                if args.separate_stars_gas:
                    bestfitspec_data += [mcsed_model.medianstarspec, mcsed_model.true_starspectrum,
                                         mcsed_model.mediannebspec, mcsed_model.true_nebspectrum]
                    bestfitspec_name += ['starspectrum', 'true_starspectrum',
                                         'nebspectrum', 'true_nebspectrum']
                T = Table(bestfitspec_data, names=bestfitspec_name)
                T.write('output/bestfitspec_fake_%05d_%s_%s.dat' % (cnt, args.sfh, args.dust_law),
                        overwrite=True, format='ascii.fixed_width_two_line')
            if args.output_dict['fluxdensity']:
                T = Table([mcsed_model.fluxwv, mcsed_model.fluxfn,
                           mcsed_model.data_fnu, mcsed_model.data_fnu_e, mcsed_model.true_fnu],
                           names=['wavelength','model_fluxdensity',
                                  'fluxdensity', 'fluxdensityerror','true_fluxdensity'])
                T.write('output/filterflux_fake_%05d_%s_%s.dat' % (cnt, args.sfh, args.dust_law),
                        overwrite=True, format='ascii.fixed_width_two_line')

    else:
        # get input data
        y, yerr, z, flag, objid, field, em, emerr, absindx, absindx_e = input_file_data

        if args.ISM_correct:
            ebv_MW = ism_igm.get_MW_EBV(args)
        else:
            ebv_MW = np.zeros(len(y))

        iv = mcsed_model.get_params()

        for yi, ye, zi, fl, oi, fd, emi, emie, indx, indxe, ebvi in zip(y, yerr, z, flag, 
                                                                        objid, field, em, emerr,
                                                                        absindx, absindx_e, ebv_MW):
            mcsed_model.filter_flag = fl
            mcsed_model.set_class_parameters(iv)
            mcsed_model.data_fnu = yi[fl]
            mcsed_model.data_fnu_e = ye[fl]
            mcsed_model.set_new_redshift(zi)
            mcsed_model.data_emline = emi
            mcsed_model.data_emline_e = emie
            mcsed_model.data_absindx = indx
            mcsed_model.data_absindx_e = indxe

            if args.sfh == 'binned_lsfr':
                # number of (useful) age grid points depends on galaxy age,
                # so the starSSP attribute should be reset each time
                mcsed_model.starSSP=None
                sfh_ages_Gyr = 10.**(np.array(mcsed_model.sfh_class.ages)-9.)
                max_ssp_age = get_max_ssp_age(args, z=zi)
                maxage_Gyr = 10.**(max_ssp_age-9.)
                binned_ssp = bin_ssp_ages(ages, starSSP, nebSSP, 
                                          emlinefluxSSP, sfh_ages_Gyr,
                                          maxage_Gyr, mcsed_model.t_birth)
                binned_ages, binned_starspec, binned_nebspec, binned_emlineflux = binned_ssp

                mcsed_model.ssp_ages = binned_ages
                mcsed_model.ssp_starspectra = binned_starspec
                mcsed_model.ssp_nebspectra = binned_nebspec
                mcsed_model.ssp_emlineflux = binned_emlineflux

            # Remove filters containing Lyman-alpha (and those blueward)
            mcsed_model.remove_waverange_filters(0., args.blue_wave_cutoff, restframe=True)
            # Remove filters dominated by dust emission, if applicable
            if not args.fit_dust_em:
                mcsed_model.remove_waverange_filters(args.wave_dust_em*1e4,1e10, 
                                                     restframe=True)

            # Only relevant if there is a nonzero E(B-V) Milky Way value to be fit
            if ebvi>1.0e-12: 
                tauISM_lam = ebvi*tauISMf(mcsed_model.wave)/1.086
                mcsed_model.tauISM_lam = tauISM_lam
            else:
                mcsed_model.tauISM_lam = None
            if args.IGM_correct:
                tauIGM_lam = tauIGMf(mcsed_model.wave,mcsed_model.redshift)
                tauIGM_lam.reshape(len(mcsed_model.wave))
                mcsed_model.tauIGM_lam = tauIGM_lam
            else:
                mcsed_model.tauIGM_lam = None

            mcsed_model.fit_model()
            mcsed_model.set_median_fit()

            if args.output_dict['sample plot']:
                mcsed_model.sample_plot('output/sample_%s_%05d_%s_%s' % 
                                        (fd, oi, args.sfh, args.dust_law),
                                        imgtype = args.output_dict['image format'])

            if args.output_dict['triangle plot']:
                mcsed_model.triangle_plot('output/triangle_%s_%05d_%s_%s' %
                                          (fd, oi, args.sfh, args.dust_law),
                                          imgtype = args.output_dict['image format'])

            mcsed_model.table.add_row([fd, oi, zi] + [0.]*(len(labels)-3))
            names = mcsed_model.get_param_names()
            names.append('Log Mass')
            names.append('SFR10')
            names.append('SFR100')
            if not mcsed_model.dust_em_class.fixed:
                names.append('fPDR')
            if mcsed_model.dust_em_class.assume_energy_balance:
                names.append('Mdust_EB')
            names.append('Ln Prob')
            if args.output_dict['fitposterior']: 
                T = Table(mcsed_model.samples, names=names)
                T.write('output/fitposterior_%s_%05d_%s_%s.dat' % (fd, oi, args.sfh, args.dust_law),
                        overwrite=True, format='ascii.fixed_width_two_line')
            if args.output_dict['bestfitspec']:
                bestfitspec_data = [mcsed_model.wave, mcsed_model.medianspec]
                bestfitspec_name = ['wavelength', 'spectrum']
                if args.separate_stars_gas:
                    bestfitspec_data += [mcsed_model.medianstarspec, mcsed_model.mediannebspec]
                    bestfitspec_name += ['starspectrum', 'nebspectrum']
                T = Table(bestfitspec_data, names=bestfitspec_name)
                T.write('output/bestfitspec_%s_%05d_%s_%s.dat' % (fd, oi, args.sfh, args.dust_law),
                        overwrite=True, format='ascii.fixed_width_two_line')
            if args.output_dict['fluxdensity']:
                T = Table([mcsed_model.fluxwv, mcsed_model.fluxfn,
                           mcsed_model.data_fnu, mcsed_model.data_fnu_e],
                           names=['wavelength','model_fluxdensity',
                                  'fluxdensity', 'fluxdensityerror'])
                T.write('output/filterflux_%s_%05d_%s_%s.dat' % (fd, oi, args.sfh, args.dust_law),
                        overwrite=True, format='ascii.fixed_width_two_line')
            if (args.output_dict['lineflux']) & (mcsed_model.use_emline_flux):
                emlines = list(mcsed_model.emline_dict.keys())
                emwaves, wht, model_fl, fl, fle = [], [], [], [], []
                for emline in emlines:
                    emwaves.append(mcsed_model.emline_dict[emline][0])
                    wht.append(mcsed_model.emline_dict[emline][1])
                    model_fl.append( mcsed_model.linefluxCSPdict[emline] )
                    fl.append( mcsed_model.data_emline['%s_FLUX' % emline] )
                    fle.append( mcsed_model.data_emline_e['%s_ERR' % emline] )
                T = Table([emwaves, wht, model_fl, fl, fle],
                          names=['rest_wavelength', 'weight', 'model_lineflux',
                                 'lineflux', 'linefluxerror'])
                T.sort('rest_wavelength')
                if len(T):
                    T.write('output/lineflux_%s_%05d_%s_%s.dat' % (fd, oi, args.sfh, args.dust_law),
                            overwrite=True, format='ascii.fixed_width_two_line')
            if (args.output_dict['absindx']) & (mcsed_model.use_absorption_indx):
                abs_names = list(mcsed_model.absindx_dict.keys())
                # each name: name, weight, modeled, measured, error
                wht, model, measure, error = [], [], [], []
                for indx in abs_names:
                    wht.append( mcsed_model.absindx_dict[indx][0] )
                    model.append( mcsed_model.absindxCSPdict[indx] )
                    measure.append( mcsed_model.data_absindx['%s_INDX' % indx] )
                    error.append( mcsed_model.data_absindx_e['%s_Err' % indx] )
                T = Table([abs_names, wht, model, measure, error],
                          names=['INDX', 'weight', 'model',
                                 'measure', 'measure_error'])
                if len(T):
                    T.write('output/absindx_%s_%05d_%s_%s.dat' % (fd, oi, args.sfh, args.dust_law),
                            overwrite=True, format='ascii.fixed_width_two_line')

            last = mcsed_model.add_fitinfo_to_table(percentiles)
            print(mcsed_model.table)
    if args.parallel:
        return [mcsed_model.table, formats]
    else:
        if args.output_dict['parameters']:
            print(mcsed_model.table)
            mcsed_model.table.write('output/%s' % args.output_filename,
                                    format='ascii.fixed_width_two_line',
                                    formats=formats, overwrite=True)
        if args.output_dict['settings']:
            filename = open('output/%s.args' % args.output_filename, 'w')
            del args.log
            filename.write( str( vars(args) ) )
            filename.close()
if __name__ == '__main__':
    main()


