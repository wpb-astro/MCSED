.. _section:outputs:

Outputs
=======

``MCSED`` can produce a number of output ascii files and plots, which
summarize the fitted parameters, show their uncertainties, and compare
the fits to the input data. All these data will be stored in a directory
entitled ``output`` with a common filename header ``FName``, where by
default, ``FName = ’test.dat’``. The user can change this default name
using the ``-o`` argument on the command line. The user can also select
which files to be return and which files to skip using the output
dictionary in ``config.py``.

.. _subsec:outputfile:

Summary File (``parameters = True``)
------------------------------------

If the ``parameters`` keyword in the output dictionary is set to
``True``, a summary output file is generated. In this file, one row will
be returned for each object. The columns of this file include at a
minimum, the field, the galaxy ID, :math:`z`, and the fitted parameters,
including the current log SFR (in :math:`M_\odot` yr\ :math:`^{-1}`), the
stellar mass (in :math:`M_\odot`), the amount of attenuation as
parameterized by :math:`E(B-V)`, and variables associated with the
description of the SFR history, attenuation law, dust emissivity, and
stellar metallicity. Each of these parameter will be printed with its 5,
16, 50, 84, and 95 percentile measurements. (These defaults can be
changed in ``config.py`` using the array ``param_percentiles``.) The
output file from test mode will also return columns of truth values for
each parameter. The default is ``parameters = True``.

.. _subsec:settingsfile:

Settings File (``settings = True``)
-----------------------------------

If the ``settings`` keyword in the output dictionary is set to ``True``,
a file will be generated that lists all the user-defined fitting
assumptions, including the choice of star formation rate history,
attenuation law, dust emission law, and all associated parameters. This
file will have the extension ``.args``. The default is
``settings = True``.

.. _subsec:posteriorfile:

Posterior Probability Distributions (``fitposterior = True``)
-------------------------------------------------------------

For each object, a file will be returned which contains the calculated
parameters and the corresponding log posterior probability distribution
for each MCMC chain. The number of chains depends on the number of
walkers and steps defined in ``config.py`` or on the command line using
the ``-nw`` and ``-ns`` arguments. The parameters listed are the same as
those in the main output file. The default for this option is
``fitposterior = False``.

.. _subsec:outputSEDs:

Best Fit SEDs (``bestfitspec = True``)
--------------------------------------

For each object, a file containing the best fit SED will be returned,
with columns of wavelength and flux density. The default wavelength
range is from 500 Å to 10 :math:`\mu`\ m in the observers frame. The
flux densities are given in units of :math:`\mu`\ Jy. The default is
``bestfitspec = True``.

.. _subsec:outputphotometry:

Observed and Modeled Flux Densities (``fluxdensity = True``)
------------------------------------------------------------

For each object, a file containing the input and best-fit photometric
flux densities is returned. The columns are the name of the emission
line, its wavelength (in Å), the measured flux density, the error on
this flux density, and the best-fit flux density, all in
:math:`\mu`\ Jy. The default for this option is ``fluxdensity = True``.

.. _subsec:outputlines:

Observed and Modeled Emission Line Fluxes (``lineflux = True``)
---------------------------------------------------------------

For each object, a file containing the observed and modeled emission
line fluxes will be returned. The list of lines are those defined by
``emline_list_dict`` in ``config.py``, and the columns returned are the
name of line, it observed monochromatic flux, the uncertainty in the
measurement, and the modeled flux. The default for this option is
``lineflux = True``.

.. _subsec:outputabsorption:

Observed and Modeled Absorption Indices (``absorption = True``)
---------------------------------------------------------------

For each object, a file containing the observed and modeled absorption
line indices will be returned. The list of features are those defined by
``absorption_index_dict`` in ``config.py``, and the columns returned are
the name of index, it observed value, the uncertainty in the
measurement, and the modeled value. The default for this option is
``absorption = True``.

Sample Plot (``sample plot = True``)
------------------------------------

Plot of all realizations of the posterior distribution for all model
parameters. Each MCMC walker is represented by an individual curve and
the :math:`x`-axis tracks the number of steps of each walker. The
default for this option is ``sample plot = False``.

SSP template spectra (``template spec = True``)
-----------------------------------------------

Save a figure of the SSP spectra used in the fit. The spectra are
weighted by the amount of time between each SSP age gridpoint. The SSP
grid will only carry ages that are relevant for the calculation: for
example, if the input galaxy sample is at redshifts :math:`z > 3`, the
SSP grid will exclude ages older than :math:`\sim 2` Gyr. The default
for this option is ``template spec = True``.

Triangle Plots (``triangle Plot = True``)
-----------------------------------------

For each object, a figure summarizing the results of the fits will be
returned. Each figure will contain of a number of different plots,
including illustrations showing 200 random MCMC solutions for the fitted
spectrum, the dust attenuation curve, and the galaxy’s star formation
rate history. Each figure will also contain contour plots illustrating
the co-variances of each parameter along with the parameters’ final
probability distribution functions. The user may select the image format
for these figures using the dictionary ``output_dict`` in ``config.py``,
with the available options being
``image format = .eps, .pdf, .pgf, .png, .ps, .raw, .rgba, .svg`` and
``’.svgz’``. By default, ``image format = .png``.