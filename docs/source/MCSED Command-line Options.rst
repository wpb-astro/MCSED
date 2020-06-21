.. _section-cmd-line:

MCSED Command-Line Options
==========================

The main script is ``run_mcsed_fit.py``. Control options for
``MCSED`` can either be entered on the command line, or stated in the
configuration file ``config.py``. To view all the command line options,
open a terminal in the ``MCSED`` directory and issue:

``python run_mcsed_fit.py -h``

You will see a drop down menu of all the available input arguments. This
menu is reproduced below.

| ``python run_mcsed_fit.py -h usage: -h [-h] [-f FILENAME] [-o OUTPUT_FILENAME] [-p] [-s SSP]`` 
| ``[-z METALLICITY] [-i ISOCHRONE] [-sfh SFH] [-dl DUST_LAW] [-nw NWALKERS] [-ns NSTEPS] [-lu LOGU]``
| ``[-de DUST_EM_TYPE] [-aeb] [-t] [-tf TEST_FIELD] [-no NOBJECTS]``

**Optional Arguments:** 

``-h, --help``: show this help message and exit 

``-f FILENAME, --filename FILENAME``: File to be read for galaxy data 

``-o OUTPUT_FILENAME, --output_filename OUTPUT_FILENAME``: Output filename for given run

``-p, --parallel``: Running in parallel (multi-core only)

``-s SSP, --ssp SSP``: SSP (single stellar population) Models, default fsps 

``-z METALLICITY, --metallicity METALLICITY``: Fixed metallicity for SSP models (0.02 is solar), False if free parameter

``-i ISOCHRONE, --isochrone ISOCHRONE``: Isochrone for SSP model, e.g. padova

``-sfh SFH, --sfh SFH``: Star formation history, e.g. constant 

``-dl DUST_LAW, --dust_law DUST_LAW``: Dust law, e.g. calzetti 

``-nw NWALKERS, --nwalkers NWALKERS``: Number of walkers for EMCEE 

``-ns NSTEPS, --nsteps NSTEPS``: Number of steps for EMCEE 

``-lu LOGU, --logU LOGU``: Ionization Parameter for nebular gas 

``-de, --dust_em``: If ``True`` or ``'DL07'``, (Draine & Li 2007) dust emission parameters are fitted

``-aeb, --assume_energy_balance``: If true, normalization of dust IR emission based on attenuation amount 

``-t, --test``: Test script with fake data 

``-tf TEST_FIELD, --test_field TEST_FIELD``: Test filters will match the given field 

``-no NOBJECTS, --nobjects NOBJECTS``: Number of test objects

``-ism, --ISM_correct_coords``: If a coordinate system is given, MW dust correction will be performed

``-igm, --IGM_correct``: If selected, Madau statistical IGM correction will be done (affecting wavelengths up to rest-frame Ly\ :math:`alpha`)