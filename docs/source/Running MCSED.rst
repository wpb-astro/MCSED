.. _sec:running-mcsed:

Running MCSED
=============

``MCSED`` calculates the likelihood of a solution in the usual fashion,
i.e.,

.. math:: \log L = -{1 \over 2} \sum_{i=1}^N \left[ \ln \left( w_i \sigma_i^2 \right) +  { w_i \left( x_i - \mu_i \right)^2 \over 2 \sigma_i^2} \right]

where :math:`\mu_i` are the modeled quantities, :math:`x_i` are the data points,
:math:`\sigma_i` are the uncertainties, and :math:`w_i` are the weights.
By default, :math:`w_i \equiv 1` for all photometric points, while emission 
lines and absorption line indices can have user-defined weights
(as described in :ref:`section:inputs`). 
The :math:`\sigma_i` terms includes contributions from the uncertainties associated 
with both the observations and the models. As pointed out in :ref:`subsec:columns`, 
the observed uncertainties may be based on the errors given in the input file 
or the minimum fractional uncertainties defined by ``phot_floor_error``,
``emline_floor_error``, and ``absorption_floor_error``. In addition, the
:math:`\sigma_i` values also include an additional term associated with
the uncertainties of the models themselves. This fractional error is defined in the file
``config.py`` via the ``model_floor_error`` keyword and is given the 
default value of :math:`\sigma_m = 0.1`. 
The final :math:`\sigma_i` term is calculated as

.. math:: \sigma_i^2 = \sigma_{obs}^2 + \left( \mu_i \sigma_m \right)^2
  

``MCSED`` can be run in 3 different modes: a live mode, where SED fits
are performed on a set of galaxies defined in an input file, a test
mode, where ``MCSED`` is run on mock galaxies with known parameters, and
a parallel mode, which makes use of a system that has a large number of
cores. We describe these below.

.. _subsec:livemode:

Live Mode
---------

Live mode is the most common mode of running ``MCSED``. In this mode,
``MCSED`` will take an input file as described in
:ref:`section:inputs` and fit each object, producing
the outputs described in :ref:`section:outputs`. To
run ``MCSED`` in live mode, type

``python run_mcsed_fit.py -f <input_file>``

with any of the additional parameters described above. (Use the ``-h``
option to see the complete list). While running ``MCSED``, a table will
appear onscreen as it is being built, as well as the mean acceptance
fraction, auto correlation steps, and number of burn-in steps for each
object.

.. _subsec:testmode:

Test Mode
---------

Test mode allows the user to explore how well ``MCSED`` can recover
model parameters under realistic observing conditions. In this mode,
``MCSED`` will construct a mock galaxy SED, simulate photometric
observations of the galaxy, and fit a model to the simulated photometry.

To run ``MCSED`` in test mode, use the ``-t`` option on the command
line:

``python run_mcsed_fit.py -t``

Since this mode involves simulating photometry for a model galaxy SED,
it is directly tied to the Skelton et al. (2014) photometric catalog
(and is currently a photometry-only mode). The test field is specified
by including the (optional) ``-tf FIELD`` option in the above command
line call, where ``FIELD`` can be one of ``aegis``, ``cosmos``,
``goodsn``, ``goodss``, or ``uds`` (the default is ``cosmos``). The
photometric filters and corresponding photometric errors are the filter
set associated with the selected CANDELS field.

The mock SED is realized from a set of “truth” model parameters
(corresponding to the parameter classes defined in ``config.py`` or
specified on the command line) and randomly drawn from the prior
distributions. The redshift of the mock galaxy is randomly drawn from a
uniform distribution spanning the range of redshifts specified by the
``test_zrange`` keyword in ``config.py``. Photometric data for the mock
galaxy are then simulated by measuring the filter flux densities from
the mock SED and perturbing those fluxes about their associated
photometric uncertainties. Finally, these mock observations are used to
fit an SED model and estimate the corresponding model parameters. Upon
completion, ``MCSED`` will return the calculated parameters from the
run, as well as the initially fixed parameters.

.. _subsec:parallelmode:

Parallel mode
-------------

Parallel mode will automatically utilize :math:`N-i` cores on a
multiprocessor system, where :math:`N` is the total number of cores on
the machine and :math:`i` is the number of cores that will not be used
in the calculation. The number of reserved cores can be specified by the
user in ``config.py`` using the ``reserved_cores`` keyword. This mode is
extremely useful for fitting large samples of galaxies.

To run ``MCSED`` in parallel mode, use the ``-p`` option on the command
line, e.g.,

``python run_mcsed_fit.py -f <input_file> -p``

along with any other desired command line options. Parallel mode can be
used with either live mode (as above) or test mode.
