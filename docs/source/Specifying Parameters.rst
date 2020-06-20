.. _section-param-specs:

Specifying Parameters
=====================

The recorded SED of a galaxy is a function of many parameters: some
associated with the universe (i.e., object’s redshift and adopted
cosmology), some dependent on a galaxy’s stars (i.e., stellar masses,
ages, metallicities, initial mass functions, etc.), and some arising
from the galaxy’s interstellar medium (i.e., attenuation from dust and
emission from ionized gas, PAH molecules, and dust). For any individual
calculation, some of these parameters will be fixed and some will be
allowed to vary. For the “free” parameters, it is usually necessary to
provide some constraints on their values, either to ensure that the
final result are physical, or to limit the degrees of freedom in a
calculation.

The easiest way to specify parameters to ``MCSED`` is through its
configuration files. The file ``config.py`` contains the top-level
information about most of the important parameters that users will want
access to. Alternatively, several of the most frequently used parameters
can be modified via command-line arguments. If one wants to change some
of ``MCSED``’s preset options, such as the cosmology or the range of
values that a parameter is allowed to take on, this can be done in
individual modules, such as ``sfr.py`` and ``dust_abs.py``. A list of
the variables that a user might want to modify appears in :ref:`sec:parameterlist`.