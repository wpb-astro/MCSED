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

| ``usage: -h [-h] [-f FILENAME] [-o OUTPUT_FILENAME] [-p] [-t] [-tf TEST_FIELD]``
| ``[-no NOBJECTS] [-s SSP] [-i ISOCHRONE] [-sfh SFH] [-dl DUST_LAW]``
| ``          [-de DUST_EM] [-aeb] [-z METALLICITY] [-nw NWALKERS] [-ns NSTEPS]``
| ``          [-lu LOGU] [-ism ISM_CORRECT_COORDS] [-igm]``


**Optional Arguments:** 

``-h, --help``: show this help message and exit 

``-f FILENAME, --filename FILENAME``: File to be read for galaxy data 

``-o OUTPUT_FILENAME, --output_filename OUTPUT_FILENAME``: Output filename for given run

``-p, --parallel``: Running in parallel?

``-t, --test``: Test mode with mock galaxies

``-tf TEST_FIELD, --test_field TEST_FIELD``: Test filters will match the given field

``-no NOBJECTS, --nobjects NOBJECTS``: Number of test objects

``-s SSP, --ssp SSP``: SSP (single stellar population) Models, default fsps 

``-i ISOCHRONE, --isochrone ISOCHRONE``: Isochrone for SSP model, e.g. padova

``-sfh SFH, --sfh SFH``: Star formation history, e.g. constant 

``-dl DUST_LAW, --dust_law DUST_LAW``: Dust law, e.g. calzetti 

``-de, --dust_em``: Dust emission class, e.g., DL07 

``-aeb, --assume_energy_balance``: If true, normalization of dust IR emission based on attenuation amount 

``-z METALLICITY, --metallicity METALLICITY``: Fixed metallicity for SSP models (0.019 is solar), False if free parameter

``-nw NWALKERS, --nwalkers NWALKERS``: Number of walkers for EMCEE

``-ns NSTEPS, --nsteps NSTEPS``: Number of steps for EMCEE

``-lu LOGU, --logU LOGU``: Ionization Parameter for nebular gas

``-ism, --ISM_correct_coords``: If a coordinate system is given, MW dust correction will be performed

``-igm, --IGM_correct``: If selected, Madau statistical IGM correction will be done (affecting wavelengths up to rest-frame Ly\ :math:`\alpha`)
