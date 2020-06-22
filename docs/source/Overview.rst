.. _section-overview:

Overview
========

``MCSED`` is a program designed to model the spectral energy
distribution (SED) of galaxies. Creating a code to fit the light from
galaxies is not hard: in fact, there are a large number of such codes
already available to the public, including GRASIL (Panuzzo et al. 2005),
CIGALE (Noll et al. 2009), GalMC (Acquaviva et al. 2011), da Cunha
et al. (2012), BayeSED (Han & Han 2014), BEAGLE (Chevallard & Charlot
2016), Prospector (Leja et al. 2017), and FIREFLY (Wilkinson
et al. 2017). ``MCSED`` offers several advantages over most of these
codes:

* **Flexibility**: The physical properties of galaxies vary greatly over cosmic time and environment, and approximations which work well locally or in star-forming galaxies may be inappropriate at high redshift or quenched systems. Moreover, our understanding of stellar evolution, the details of mass loss, and the products of binary evolution are constantly improving. As a result, it is difficult to develop a general purpose SED-fitting code that works for all objects. ``MCSED`` is designed with flexibility in mind: while there are many options already built into the program, the code is built in modules, and one can easily substitute or add in new datasets or algorithms.

* **Multiwavelength Capability**: ``MCSED`` is built to fit a galaxy’s full SED, from the far-UV to the far-IR. Some of the physical processes it is designed to model include continuum emission from stars, continuum and line-emission from ionized gas, attenuation from dust, and mid- and far-IR emission from dust and polycyclic aromatic hydrocarbons (PAHs). That is not to say that one must have data covering the entire electromagnetic spectrum. For example, if only the UV through near-IR is available, ``MCSED`` can fit parameters for the starlight and dust attenuation, but ignore variables associated with the dust emissivity.

* **Variable Input Data**: Many SED codes fit data from a collection of broad- and intermediate-band filters. ``MCSED`` can do this, but it can also incorporate constraints provided by emission-lines fluxes and absorption-line spectral indices. One can easily reconfigure the code to accept virtually any type of constraint.

* **Efficiency**: ``MCSED`` performs its calculations by creating a complex stellar population (CSP) out of a linear combination of simple-stellar populations (SSPs) using an efficient Markov Chain Monte Carlo algorithm. It is very quick, and takes advantage of parallel processing.

At present, the major limitations of ``MCSED`` are:

* The fitted objects must have known redshifts. ``MCSED`` is not designed to estimate photometric redshifts or marginalize over a range of possible redshift values.

* AGN emission is not included in the modeling. ``MCSED`` only models the light from stars, the emission from gas which is photoionized by stars, and long-wavelength emission from energy re-radiated by dust.

* Chemical evolution is not followed. In its modeling, ``MCSED`` assumes only a single metallicity for all of a galaxy’s stellar populations.
