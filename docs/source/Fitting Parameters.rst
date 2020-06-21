.. _section:parameters:

Fitting Parameters
==================

``MCSED`` contains numerous options for fitting a galaxy’s spectral
energy distribution. Many of these options can be specified either on
the command line or entered directly in the ``config.py`` configuration
script, including 5 dust attenuation laws and 6 star formation rate
histories. These options (and others) are described below.

.. _subsec:attenuation:

Dust Attenuation Laws
---------------------

Several commonly used prescriptions for dust attenuation are included in
``MCSED``. These include the locally-derived laws of Calzetti
et al. (2000) and Conroy et al. (2010), the :math:`z \sim 2` relation
found by Reddy et al. (2015), and the generalization of the Calzetti
et al. (2000) prescription formulated by Noll et al. (2009) and Kriek &
Conroy (2013). Note that these laws describe the attenuation of the
*stellar continuum* of a galaxy. As shown by Calzetti et al. (2000) and
modeled by Charlot & Fall (2000), the extinction surrounding regions of
active star formation (including the ionized gas) is generally greater
than that for the older stellar populations. ``MCSED`` recognizes this
fact by allowing the user to model SEDs with two different values for
the attenuation, one for the ionized gas and stellar populations younger
than :math:`t_{\rm birth}`, and one for the galaxy’s older population
with

.. math:: E(B-V)_{\rm old} = \eta \, E(B-V)_{\rm young}

In ``MCSED``, the default values for :math:`t_{\rm birth}` and
:math:`\eta` are set in ``config.py`` to ``t_birth`` = 7 (log years) and
``EBV_old_young`` = 0.44. (The latter number comes the analysis of local
starburst galaxies by Calzetti et al. 2000). These values can easily be
reset in ``config.py``.

The dust attenuation law can be specified on the command line using the
``-dl`` flag, or by setting the variable ``dust_law`` in ``config.py``.
Within ``MCSED`` are four pre-coded attenuation laws (``’calzetti’``,
``’conroy’``, ``’reddy’``, and ``’noll’``) and one Galactic extinction
curve (``’cardelli’``). For each of these prescriptions, the amount of
extinction is parameterized by the differential reddening,
:math:`E(B-V)`, which is initialized in ``MCSED`` through the parameter
``EBV`` in ``dust_abs.py``. As a default, ``EBV`` is allowed to take on
values between ``[-0.05, 1.50]``.

Note that many attenuation laws are written in terms of :math:`R_V`,
which describes the ratio of total to differential extinction in the
:math:`V`-band. Different laws prefer different values of :math:`R_V`:
for example, in the prescription of Calzetti et al. (2000),
:math:`R_V = 4.05`, while the laws of Reddy et al. (2015) and Conroy
et al. (2010), prefer values of :math:`R_V = 2.505` and 2.0,
respectively. In ``config.py``, :math:`R_V` is set to a negative value
(``Rv = -1``), as this instructs ``MCSED`` to use the value of
:math:`R_V` most commonly associated with the attenuation curve chosen
by the user (which is specified in each class defined in
``dust_abs.py``). If the user wishes to use a different value of
:math:`R_V`, they can set ``Rv`` to the desired (positive) number in
``config.py``.

.. _subsubsec:calzetti:

``dust_law = ’calzetti’``
~~~~~~~~~~~~~~~~~~~~~~~~~

The Calzetti attenuation law, which is described (using slightly
different notations) in Calzetti et al. (2000) and Calzetti (2001), was
derived using :math:`\sim 40` local UV-bright starburst galaxies. It is
defined as

.. math::

   k(\lambda) =
     \begin{cases}
     \begin{split}
       2.659(-2.156+1.509/\lambda-0.198/\lambda^2\\+0.011/\lambda^3)+R_V \end{split}& \text{ for $0.12~\mu{\rm m} \leq\lambda\leq0.63~\mu$m} \\
       2.659(-1.857+1.040/\lambda)+R_V & \text{ for $0.63~\mu {\rm m} \leq\lambda\leq2.20~\mu$m}
     \end{cases}

where :math:`k(\lambda)=A(\lambda)/E(B-V)`, and :math:`R_V = 4.05`.
There is no 2175 Å attenuation excess in the Calzetti attenuation curve,
and no additional parameters are associated with its shape.
``MCSED`` uses ``’calzetti’`` as its default attenuation law.

.. _subsubsec:conroy:

``dust_law = ’conroy’``  or  ``’cardelli’``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Conroy et al. (2010) noticed that the attenuation of local disk galaxies
within the mass range :math:`9.5 < \log(M_*/M_\odot) < 10` could be
modeled by slightly modifying the Milky Way extinction law of Cardelli
et al. (1989). This law takes the functional

.. math:: A(\lambda)/A_V=a(x)+b(x)/R_V

where :math:`x=1/\lambda` (in units of :math:`\mu {\rm m}^{-1}`), and
the values of :math:`a(x)` and :math:`b(x)` are as follows: 6pt

:math:`\bullet` In the Extreme UV between
:math:`8.0~\mu{\rm m}^{-1} < x < 10.0~\mu{\rm m}^{-1}`
:math:`(0.1~\mu{\rm m} < \lambda < 0.125~\mu {\rm m})`:

.. math:: a(x) = -1.073 - 0.628(x-8)+0.137(x-8)^2-0.070(x-8)^3,

.. math:: b(x)=13.670+4.257(x-8)-0.420(x-8)^2+0.374(x-8)^3

:math:`\bullet` In the UV between
:math:`3.3~\mu{\rm m}^{-1} < x < 8.0~\mu{\rm m}^{-1}`
:math:`(0.125~\mu{\rm m} < \lambda < 0.17~\mu {\rm m})`:

.. math:: a(x) = 1.752 - 0.316x -\frac{0.104B}{(x-4.67)^2+0.341}+f_a,

.. math:: b(x)=-3.09+1.825x+\frac{1.206B}{(x-4.62)^2+0.263}+f_b.

where

.. math::

   \begin{cases}
   \begin{split}
       f_a = -0.04473(x-5.9)^2 - 0.009779(x-5.9)^3 \\
       f_b=0.213(x-5.9)^2 +0.121(x-5.9)^3 \end{split} &
       \text{ for $x > 5.9~\mu$m$^{-1}$} \\
   \begin{split}
       f_a=\left(\displaystyle\frac{3.3}{x}\right)^6(-0.0370+0.0469B-0.601B/R_V+0.542/R_V)  \\
       f_b = 0 \end{split} & \text{ for $x < 5.9~\mu$m$^{-1}$}
   \end{cases}

In the Optical between
:math:`1.1~\mu{\rm m}^{-1} < x <3.3~\mu{\rm m}^{-1}`:
:math:`(0.30~\mu{\rm m} < \lambda < 0.91~\mu {\rm m})`

.. math:: y=x-1.82

.. math::

   \begin{split}
       a(x)=1+0.177y-0.504y^2-0.0243y^3+0.721y^4+0.0198y^5\\-0.775y^6+0.330y^7,
       \end{split}

.. math::

   \begin{split}
       b(x)=1.413y+2.283y^2+1.072y^3-5.384y^4-0.622y^5+\\5.303y^6-2.090y^7
   \end{split}

In the Near-IR between
:math:`0.3~\mu{\rm m}^{-1} < x <1.1~\mu{\rm m}^{-1}`
:math:`(0.91~\mu{\rm m} < \lambda < 3.33~\mu {\rm m})`:

.. math:: a(x) = 0.574 x^{1.61}

.. math:: b(x) = -0.527 x^{1.61}

In the original formulation of Cardelli et al. (1989), :math:`B=1` and
:math:`f_a` = 0 between
:math:`3.3~\mu{\rm m} < x < 5.9~\mu`\ m\ :math:`^{-1}`, and over the
years, most extragalactic applications of the Cardelli law have adopted
the mean Milky Way total-to-differential extinction ratio of
:math:`R_V = 3.1`. In the Conroy et al. (2010) modification of this law,
the attenuation in nearby galaxies is best reproduced using
:math:`B = 0.8`, the expression for :math:`f_a` is a weak function of
:math:`B` and :math:`R_V`, and :math:`R_V \approx 2.0`.

In ``MCSED``, if ``dust_law = ’conroy’``, then by default, ``Rv = 2.0``
and ``B = 0.8``; if ``dust_law = ’cardelli’``, then ``Rv = 3.1`` and
``B = 1.0``. As noted above, ``Rv`` can be changed in ``config.py``; the
parameter ``B`` can be modified in ``dust_abs.py``.

.. _subsubsec:reddy:

``dust_law = ’reddy’``
~~~~~~~~~~~~~~~~~~~~~~

Reddy et al. (2015) derived an attenuation law for :math:`1.4 < z < 2.6`
galaxies selected on the basis of emission-line detections on HST grism
frames. The equations for this curve are

.. math::

   k(\lambda) =
     \begin{cases}
     \begin{split}
       -5.726+4.004/\lambda-0.525/\lambda^2\\+0.029/\lambda^3+R_V\end{split} & \text{ for $0.15~\mu {\rm m} \leq\lambda\leq0.60~\mu$m} \\
       \begin{split}
           -2.672-0.010/\lambda+1.532/\lambda^2\\
       -0.412/\lambda^3+R_V \end{split}& \text{ for $0.60~\mu {\rm m} \leq\lambda\leq2.20~\mu$m}
     \end{cases}

where :math:`k(\lambda)=A(\lambda)/E(B-V)` and :math:`R_V = 2.505`.

.. _subsubsec:noll:

``dust_law = ’noll’``
~~~~~~~~~~~~~~~~~~~~~

Noll et al. (2009) and then Kriek & Conroy (2013) developed a
generalization of the Calzetti attenuation curve using large samples of
:math:`0.5 < z < 2.0` galaxies. Their formalism perturbs the Calzetti
law, allows for a steeper/shallower extinction curve in the far-UV, and
admits the possibility of a 2175 Å bump. If we let :math:`\delta`
describe the deviation from the slope of the Calzetti attenuation law,
and let :math:`E_b` represent the strength of the 2175 Å bump, then:

.. math::

   \displaystyle k^\prime(\lambda) = \Big\{ k(\lambda) + D(\lambda) \Big\}
   \left( \frac{\lambda}{\rm 5500~Å} \right)^\delta 
   \label{eq:k}

where :math:`k(\lambda)` is the Calzetti wavelength dependence of
attenuation and :math:`D(\lambda)` is a Lorentzian-like “Drude” profile

.. math::

   D(\lambda) = E_b \, \frac{\left( \lambda \Delta\lambda \right)^2}{\left( \lambda^2 - \lambda_0^2
   \right)^2 + \left( \lambda \Delta\lambda \right)^2}
   \label{eq:bump}

Positive (negative) values of :math:`\delta` represent UV attenuation
curves that are shallower (steeper) than that modeled by Calzetti, while
positive values of :math:`E_b` indicate the presence of a 2175 Å bump.
In ``MCSED``, the default for :math:`\delta` is ``delta`` = 0 and the
variable is allowed to range between ``[-1., +1.]``. The default on the
bump strength is ``Eb`` = 0, with a default range from ``[-0.2, +6.0]``.
If the input data are not conducive for a measurement of :math:`E_b` (or
:math:`\delta`), these variables can be fixed by restricting their
ranges to that of a delta function.

.. _subsubsec:otherdust:

Other Attenuation Laws
~~~~~~~~~~~~~~~~~~~~~~

The attenuation laws pre-coded into ``MCSED`` are outlined above. Users
who wish to add a new dust law need only to define a new class in
``dust_abs.py``, following the same procedure used for the existing
laws, and call their new class via the ``dust_law`` keyword in
``config.py``. If the new class does not include a definition of
:math:`R_V`, one must be specified in ``config.py`` as described above.

.. _subsec:dustemission:

Dust Emission
-------------

``MCSED`` is capable of modeling the reprocessed radiation of warm and
cold dust in the interstellar medium. At present, ``MCSED`` contains
only one prescription for this long wavelength emission: the
Spitzer-based silicate-graphite-PAH model of Draine & Li (2007). However
it is relatively straightforward to include other models for this
component of the SED.

.. _subsubsec:draine-li:

``dust_em = ’DL07’``
~~~~~~~~~~~~~~~~~~~~

The Draine & Li (2007) prescription for mid- and far-IR dust emissivity
is built upon the models of Siebenmorgen & Kruegel (1992), Li & Draine
(2001), and Weingartner & Draine (2001). This model includes emission
lines from 25 PAH features between 3 and :math:`15~\mu`\ m (modeled as
Drude profiles), absorption and emission from PAH ions, and emission
from carbonaceous and silicate particles with a range of sizes. The
model, which has three free parameters, is laid out in full detail by
Draine & Li (2007).

Numerically, the model is defined via :math:`U`, which is the scale
factor between the interstellar radiation field found in the solar
neighborhood (Mathis et al. 1983) and that of the galaxy being modeled.
Dust emission is divided into two components: that produced from dust
which is heated by starlight with energies :math:`U < U_{\rm min}`, and
dust heated by starlight with :math:`U_{\rm min} < U < U_{\rm max}`. The
dust mass is then related to :math:`U` via

.. math::

   \frac{d{M_{\rm dust}}}{dU} = (1-\gamma){M_{\rm dust}}\delta(U-{U_{\rm min}}) + \gamma {M_{\rm dust}}\frac{\alpha-1}{U_{\rm min}^{1-\alpha}-U_{\rm max}^{1-\alpha}}U^{-\alpha}

where :math:`\alpha` and :math:`U_{\rm max}` are assigned their Milky
Way values of 2 and :math:`10^6`, respectively (Draine et al. 2007). The
three free parameters are therefore, :math:`U_{\rm min}`,
:math:`\gamma`, and :math:`q_{\rm PAH}`. This last variable represents
the percentage of the dust mass made up of PAH molecules containing less
than 1000 carbon atoms.

``MCSED`` creates two dust emission spectral components in the
wavelength range :math:`1~\mu{\rm m} \leq \lambda \leq 10,000~\mu`\ m by
performing a 2D interpolation within a grid of 22 unequally spaced
values of :math:`U_{\rm min}` between
:math:`0.1 < U_{\rm min} < 25.0` and 7
unequally spaced values of :math:`q_{\rm PAH}` between
:math:`0.47 < q_{\rm PAH} < 4.58` (Draine &
Li 2007).

The dust emission spectrum is computed in units of Jy cm\ :math:`^2`
sr\ :math:`^{-1}` H\ :math:`^{-1}`, where H represents a hydrogen
nucleon. If we assume that the emission is isotropic, this can be
converted into flux densities per solar mass of gas in units of
:math:`\mu`\ Jy M\ :math:`_\odot`\ :math:`^{-1}`. If :math:`M_{\rm gas}`
is the total (hydrogen) gas mass in the galaxy, :math:`p_\nu^{(0)}` is
the flux of/mass from dust heated by
:math:`U=U_{\rm min}`, and :math:`p_\nu` is the
flux/mass from dust heated by photons with
:math:`U_{\rm min}<U<U_{\rm max}`,
then the total dust emission from the galaxy is given by (Draine & Li 2007)

.. math::

   F_{\rm dust} = \frac{M_{\rm gas}}{4\pi} \Big[(1-\gamma)p_\nu^{(0)}(U_{\rm min},q_{\rm PAH}) + \gamma p_\nu(U_{\rm min},q_{\rm PAH},U_{\rm max}=10^6,\alpha=2) \Big]

This model assumes a constant dust-to-gas mass
ratio for each value of :math:`q_{\rm PAH}`, all very near :math:`0.01`.
When the dust-to-gas mass ratio is interpolated over
:math:`q_{\rm PAH}`, expressing the dust emission as a function of
:math:`M_{\rm dust}` is simple:

.. math::

   F_{\rm dust} = \frac{M_{\rm dust}}{4\pi}\frac{M_{\rm gas}}{M_{\rm dust}} \Big[(1-\gamma)p_\nu^{(0)} + \gamma p_\nu \Big]

As detailed in Draine et al. (2007) and in the equation above, the total dust mass,
:math:`M_{\rm dust}`, is used to normalize the dust emission spectrum.
This normalization can be a free parameter completely determined by
far-IR photometry, or linked to the amount of dust attenuation via the
assumption of energy balance. In theory, energy balance should
apply, as the energy attenuated should equal that emitted. However,
because the former measurement is sight-line dependent, while the latter
is generally isotropic, individual objects may appear to violate this
principle.

By default, the fitting of dust emission in ``MCSED`` is turned off, all
measurements from filters with rest-frame wavelengths longward of
``wave_dust_em`` = 2.5 \ :math:`\mu`\ m (set in ``config.py``) are
discarded, and the Draine & Li (2007) parameters are set to ``umin`` =
2.0, ``gamma`` = 0.05, and ``qpah`` = 2.5 with
``assume_energy_balance = False``. To instruct ``MCSED`` to fit the dust
emission, users can invoke the command-line option ``-fd`` or set
``fit_dust_em = True`` in ``config.py``. The additional command-line
option ``-aeb`` will then instruct ``MCSED`` to assume energy balance in
the SED calculation. (This last step can also be done by setting the
``assume_energy_balance = True`` in ``config.py``.)

If dust emission is fit, ``MCSED`` returns the 3-4 fitted parameters
(:math:`U_{\rm min}`, :math:`\gamma`, :math:`q_{\rm PAH}`, and :math:`M_{\rm dust}` if ``assume_energy_balance=False``) and 1-2 derived values: the dust mass (:math:`M_{\rm dust}`, if ``assume_energy_balance=True``) and the fraction of
the total dust luminosity that is radiated by dust grains in regions
where :math:`U > 100` (similar to :math:`\gamma` but for the hardest
radiation).

:math:`U_{\rm min}`, :math:`\gamma`, and :math:`q_{\rm PAH}` are defined
in ``dust_emission.py``, with defaults limits of ``[0.1, 25.0]`` for
``umin``, ``[0.0, 1.0]`` for ``gamma`` = 0.05, and ``[0.47, 4.58]`` for
``qpah``. The default on the dust mass is ``mdust`` = 7.0 (in log solar
units), and its range is ``[4.5, 10.0]``.

.. _subsubsec:otherdustemission:

Other Laws for Dust Emission
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Users who wish to add a different law will need to define a new class in
``dust_emission.py``, following the same procedure used for the existing
laws, and call their new class via the ``dust_em`` keyword in
``config.py``.

.. _subsec:SFRH:

Star Formation History
----------------------

The choice of star formation history can be specified at the command
line using the ``-sfh`` argument or by setting ``sfr`` in ``config.py``.
``MCSED`` contains 6 built-in options which describe how the star
formation rate in a galaxy evolves with time: five analytic expressions,
and one defined via a series of user-defined age bins. Both the
parameterized and binned versions of SFR history can be found in
``sfh.py``, along with the definitions of their parameters.

We detail the possible star formation rate histories below. Note that
``MCSED`` uses the time variable :math:`t` to represent the lookback
time from the epoch of observation in Gyr. In other words, :math:`t = 0`
is when the galaxy is being observed (redshift :math:`z`), and positive
values of :math:`t` represent times earlier in the history of the
universe.

``sfh = ’constant’``
~~~~~~~~~~~~~~~~~~~~

The simplest star formation history is ``sfr = ‘constant’``; the only
parameters are the Base-10 logarithm of the star formation rate in
(:math:`\log M_{\odot}` yr\ :math:`^{-1}`) and the age of the galaxy in
:math:`\log` Gyr. The default limits on the :math:`\log{\rm{SFR}}` are
``logsfr`` = ``[-3,3]``, while the age limits on the galaxy go from 1
Myr to the time difference between the epoch of observation (redshift
:math:`z`) and :math:`z = 20`. ``MCSED`` uses this as its default star
formation history.

``sfh = ’exponential’``
~~~~~~~~~~~~~~~~~~~~~~~

A popular parameterization of a galaxy’s SFR history is through an
exponential, i.e.,

.. math:: {\rm SFR} = A \exp^{-t /10^\tau}

The ``sfr = ’exponential’`` option has three parameters: the age of the
galaxy (in :math:`\log` Gyr), the e-folding timescale :math:`\tau` (in
:math:`\log` Gyr), and the normalization constant :math:`A` (in
:math:`\log M_{\odot}` yr\ :math:`^{-1}`). The default limits for age are
the same as for the ``’constant’`` case: 1 Myr to the time between the
observed redshift and :math:`z = 20`. The default limits on the
e-folding timescale, ``tau`` are ``[-3.5,3.5]`` (in log Gyr) and the
range of ``logsfr`` normalizations go from ``[-3,3.0]`` in
:math:`\log M_{\odot}` yr\ :math:`^{-1}`.

The default behavior is the SFR decaying exponentially with lookback time. This can be changed (i.e.
exponential growth with lookback time) by changing the ``sign`` variable to -1 in the ``exponential`` star
formation history class in ``sfh.py``.

``sfh = ’double_powerlaw’``
~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``’double_powerlaw’`` star formation history is a five-parameter
formulation of Behroozi et al. (2013), with

.. math:: {\rm SFR}(t) = A \left[\left(\frac{t}{\tau}\right)^B + \left(\frac{t}{\tau}\right)^{-C}\right]^{-1}

In this parameterization, :math:`B` describes the rate of increase in a
galaxy’s SFR early in its history, :math:`C` gives the rate at which the
SFR declines, :math:`\tau` is the lookback time corresponding to the
transition between these two phases, and :math:`A` is the normalization.
By default, the limits for the variables ``B`` and ``C`` range from
``[0,5]``, ``tau`` is allowed to vary from ``[-3.0,1.0]`` (in Gyr), and
the normalization ``a`` (in :math:`\log M_{\odot}` yr\ :math:`^{-1}`) can
go from ``[-1.0,5.0]``. As above, the age of the galaxy can vary between
1 Myr and the time between :math:`z = 20` and the epoch of observation.

``sfh = ’burst’``
~~~~~~~~~~~~~~~~~

The ``’burst’`` SFR history models a galaxy with a single burst of star
formation superposed on a constant star formation rate, i.e.,

.. math:: \text{SFR} = 10^{Sc} + \left\{ \frac{Sb  \times 10^{Sc}}{2 \pi \sigma} \exp{\left[-0.5\frac{(\log{t} - a_b)^2}{\sigma^2}   \right]} \right\} \\

The five parameters of this expression are the Base-10 log of the
quiescent star formation rate (:math:`Sc`), the burst strength
(:math:`Sb`) relative to the quiescent rate, the duration of the burst
as measured by the log-normal standard deviation (:math:`\sigma`), the
time since the burst (:math:`a_b`), and the age of the galaxy. In
``MCSED``, :math:`Sc` is represented by ``logsfr`` and has default
limits between ``[-3,3]`` in :math:`\log M_{\odot}` yr\ :math:`^{-1}`,
while the strength of the burst, ``burst_strength``, can vary from
``[1,10]``. The defaults on the burst length, defined via
``burst_sigma``, are ``[0.003,0.5]`` (in log Gyr), and the time since
the burst, ``burst_age``, goes from ``[6.0,7.5]`` in log years. (Note
that this last parameter is optimized for recent bursts of star
formation, and is given in log years, not Gyrs.) Once again, the age of
the galaxy can vary from 1 Gyr to the time between :math:`z = 20` and
the redshift of the galaxy.

``sfh = ’polynomial’``
~~~~~~~~~~~~~~~~~~~~~~

The most complex analytical expression for the SFR history of a galaxy
is ``’polynomial’``. This option stores as parameters the Base-10
logarithm of the SFR (in :math:`\log M_{\odot}` yr\ :math:`^{-1}`) at
certain (fixed) pivot points of lookback time (in :math:`\log` years),
with the degree of the polynomial being one less than the number of
pivot points. The pivot points in ``MCSED`` are stored in the array
``age_locs`` (which also defines the degree of the polynomial). The
default values for this array are ``[6.5,7.5,8.5]``, which implies a
quadratic relation. The user can change these default log ages in the
``sfh.py`` module.

For a given set of SFRs at the :math:`n+1` pivot points,
``MCSED`` determines the :math:`n`-th degree polynomial that goes
through those points using ``NumPy``\ ’s least-squares method. This
calculation is done in log space by setting up a matrix where each
column is the array of pivot points with respect to the median pivot
(:math:`\vec x=0`) raised to the power of the polynomial degree for the
leftmost column all the way down to :math:`0` for the rightmost column.
We then solve for the coefficients via

.. math::

   \begin{bmatrix}
   x_1^n & x_1^{n-1} & \dots & 1\\ 
   x_2^n & x_2^{n-1} & \dots & 1\\ 
   \vdots & \vdots & \ddots & \vdots\\
   x_{n+1}^n & x_{n+1}^{n-1} & \dots & 1
   \end{bmatrix}
   \begin{bmatrix}
   y_1\\ 
   y_2\\ 
   \vdots\\
   y_{n+1}
   \end{bmatrix}
   =
   \begin{bmatrix}
   v_1\\ 
   v_2\\ 
   \vdots\\
   v_{n+1}
   \end{bmatrix}

where the parameters of the model are labeled :math:`v_1`, :math:`v_2`,
:math:`\dots`, :math:`v_{n+1}`. If we then let :math:`a_{\textrm{mid}}`
be the median age of the pivot points and set the input lookback time
:math:`t` in Gyr, we can determine the SFR via

.. math::

   \begin{aligned}
       P(x) &= y_1 x^n + y_2 x^{n-1} + \dots + y_{n+1} \\
       \text{SFR} &= 10^{P\left(\log{t}+9-a_{\textrm{mid}}\right)}\end{aligned}

``sfh = ’binned_lsfr’``
~~~~~~~~~~~~~~~~~~~~~~~

In the ``'binned lsfr'`` option, the star formation rate history of a galaxy is divided into a series of age bins,
with the SFR internal to each bin assumed to be constant. ``MCSED`` fits for the log SFR within each time bin.
``MCSED`` has six (log) age bins defined by the ``ages`` array within ``sfh.py``; as a default, the bins are define as
ages = [8.0, 8.5, 9.0, 9.5, 9.8, 10.12]. These values are adopted from Leja et al. (2017) and are motivated
by physical considerations. The user can easily modify these age bins by editing 
the ages keyword defined in the ``binned_lsfr``
class in ``sfh.py``. While these ages extend to the age of the universe, only the SSP spectra
that are younger than the age of the galaxy will contribute to the nal SED model.
Since the SFR is assumed to be constant within each age bin, the computational efficiency can be
improved by collapsing the SSP grid and summing the spectra within each SFH time bin via a weighted
average, where the weights are determined by the amount of time between the SSP ages. This step is
automatically implemented and uses the ``bin_ssp_ages`` function defined in ``ssp.py``.

.. _subsubsection:otherSFR:

Other SFR Prescriptions
~~~~~~~~~~~~~~~~~~~~~~~

Users who wish to add a different star formation history prescription
need only to define a new class in ``sfh.py``, following the same
procedure used for the existing formulations, and call their new class
via the ``sfh`` keyword in ``config.py``.

.. _subsec:metallicity:

Stellar Metallicity
-------------------

``MCSED`` uses a library of SSP spectra which are layed out over a
two-dimensional grid in age and metallicity. In forming the CSP of a
galaxy, one can either fix the metallicity to some value, or allow
metallicity to be a free parameter. This choice is accomplished either
on the command line via the ``-z`` option, or by setting the variable
``metallicity`` in ``config.py``; a real value fixes the metallicity
:math:`Z` (where :math:`Z_\odot = 0.019`), while the boolean ``False``
instructs ``MCSED`` to treat the stellar metallicity as a free model
parameter and solve for the most likely abundance. If the stellar
metallicity is a free parameter, the model parameter is
:math:`\log(Z/Z_\odot)` with a uniform prior spanning the range
[:math:`-1.98`, 0.2]. The prior can be adjusted by editing the
``stellar_metallicity`` class in the ``metallicity.py`` file.

By default, ``MCSED`` sets ``metallicity`` to a fixed value of 0.0077
(:math:`\sim 40\%` solar). (A near-future option will allow the
metallicity to be tied to a galaxy’s stellar mass, as suggested by Peng
& Maiolino 2014 and Ma et al. 2016.) In either case,
``MCSED`` introduces a small metallicity scatter into the calculation
using a Gaussian kernel with dispersion
:math:`\sigma = 0.15 \, \log (Z / Z_{\odot})`. This minimizes potential
biases in stellar mass and other inferred quantities (see Mitchell
et al. 2013).

Note that in most cases, broad- and intermediate-band photometry will
provide (at best) only a weak constraint on metallicity. For stronger
constraints, one needs to include emission-line fluxes (or absorption
line indices) in the SED fits.

The current version of ``MCSED`` has no provision for following the
chemical evolution of the various stellar populations within a galaxy.
Only a single metallicity is used in the fits.

.. _subsec:ionization-param:

Ionization Parameter
--------------------

The ionization parameter, which measures the number ionizing photons per
atom, is important for modeling the nebular emission from galaxies. If
the ionization parameter is low, an electron will likely recombine into
a atom of ionization state :math:`i` before that atom encounters a
photon capable of producing an additional ionization. If the value of
:math:`U` is large, then an atom may undergo multiple ionizations before
a recombination occurs. Consequently, in order to properly model a
galaxy’s emission lines, ``MCSED`` needs some estimate of this important
parameter.

Currently, ``MCSED`` applies the same ionization parameter to all
galaxies in the input file. The default value of ``logU`` :math:`= -2.5`
can be changed in the command-line with the option ``-lu value`` or
modified directly in ``config.py``. The current limits of the nebular
models extend from :math:`-4 < \log U < -1` (Byler et al. 2017).

IGM Correction
--------------

``MCSED`` includes an option to apply a statistical correction to
wavelengths shortward of (rest frame) 1216 Å, in order to account of
absorption by the intergalactic medium (IGM). To maximize speed, this
correction is computed by means of 2-D linear interpolation in a grid of
IGM optical depths in observed-frame wavelength-redshift space. The
table itself was generated using the equations of Madau (1995), and
accounts for both Lyman line and continuum absorption. If
:math:`z_\mathrm{em}` is the redshift of the source and
:math:`\lambda_\mathrm{obs}` the observed wavelength, the correction for
continuum aborption is

.. math:: \tau_\mathrm{cont} \approx 0.25x_\mathrm{c}^3\left(x_\mathrm{em}^{0.46} - x_\mathrm{c}^{0.46} \right) + 9.4x_\mathrm{c}^{1.5}\left(x_\mathrm{em}^{0.18} - x_\mathrm{c}^{0.18} \right) - 0.7x_\mathrm{c}^3\left(x_\mathrm{em}^{-1.32} - x_\mathrm{c}^{-1.32} \right) - 0.023\left(x_\mathrm{em}^{1.68} - x_\mathrm{c}^{1.68} \right)

where :math:`x_\mathrm{em} = 1 + z_\mathrm{em}` and
:math:`x_c = \lambda_\mathrm{obs} / 911.75` Å.

Bound-bound absorptions are a bit more complex as the number of Lyman
lines that contribute to the opacity depends on the wavelength. If we
let :math:`A_n` be the optical depths coefficients of the Lyman lines
(derived from the Einstein-A coefficients and curve-of-growth analyses
with constant Doppler parameter :math:`b = 35` km s\ :math:`^{-1}`),
then the total optical depth is

.. math::

   \left\{ \begin{matrix}
       \sum_{j=2}^{n_\mathrm{max}} A_j \left( \frac{\lambda_\mathrm{obs}}{\lambda_j}\right)^{3.46} & \mathrm{ if }~\lambda_{\mathrm{obs}} < \lambda_{n_\mathrm{max}}(1+z_\mathrm{em})\\ 
       \sum_{j=2}^i A_j \left( \frac{\lambda_\mathrm{obs}}{\lambda_j}\right)^{3.46} & \mathrm{ if }~\lambda_{i+1}(1+z_\mathrm{em}) < \lambda_\mathrm{obs}<\lambda_{i}(1+z_\mathrm{em}) \\ 
       0 & \mathrm{ if }~\lambda_2(1+z_\mathrm{em}) < \lambda_{\mathrm{obs}}
       \end{matrix}\right.

``MCSED``\ ’s optical depth table includes all coefficients :math:`A_n`
from :math:`n=2` (Ly\ :math:`\alpha`) to :math:`n=40`.

The total optical depth is the sum of the line and continuum optical
depths. The user can opt to have IGM corrections applied with the
command-line option ``-igm`` (no other arguments) or by selecting
``IGM_correct = True`` in ``config.py``.

If statistical corrections for IGM absorption are insufficient, the user
can opt to exclude all filters with significant throughputs shortward of
a given wavelength from the SED fit. By default, no filters are re-
moved, but the user can specify a minimum wavelength used by ``MCSED`` by specifying ``blue_wave_cutoff = wavelength`` in ``config.py`` (rest-frame wavelength in Angstroms).

ISM Corrections
---------------

``MCSED`` uses the Python package
`dustmaps <https://dustmaps.readthedocs.io/en/latest/installation.html>`__
(Green 2018) to de-redden a computed SED for Milky Way extinction. To do
this, the program uses the Schlegel, Finkbeiner, & Davis (1998) 2-D dust
maps recalibrated by Schlafly & Finkbeiner (2011), and the Fitzpatrick
(1999) extinction curve with :math:`R_V = 3.1`.

The user can choose to let ``MCSED`` perform Milky Way dust corrections
by including the command-line option ``-ism coord_system`` or by setting
``ISM_correct_coords`` in ``config.py``. In either case, the user must
provide the coordinates of the objects in the input file; the available
options for ``coordinate_system`` include

.. table:: Coordinate Systems

   +---------------------------+-------------------------+--------------------+
   | Coordinate System         |   Coordinate System     | Coordinate System  |
   +===========================+=========================+====================+
   | altaz                     | barycentrictrueecliptic | cirs               |
   +---------------------------+-------------------------+--------------------+
   | fk4                       | fk4noeterms             | fk5                |
   +---------------------------+-------------------------+--------------------+
   | galactic                  | galacticlsr             | galactocentric     |
   +---------------------------+-------------------------+--------------------+
   | gcrs                      | geocentrictrueecliptic  | hcrs               |
   +---------------------------+-------------------------+--------------------+
   | heliocentrictrueeclipctic | **icrs**                | itrs               |
   +---------------------------+-------------------------+--------------------+
   | lsr                       | precessedgeocentric     | supergalactic      |
   +---------------------------+-------------------------+--------------------+

The coordinates should be provided in two columns labeled ‘``C1``’ and
‘``C2``,’ with the longitudinal coordinate first and the latitude second
(for example, RA :math:`\rightarrow` ‘``C1``’ and declination
:math:`\rightarrow` ‘``C2``’). The user is referred to the AstroPy `SkyCoord package
<https://docs.astropy.org/en/stable/api/astropy.coordinates.SkyCoord.html>`__ for details on the different coordinate systems and the coordinate orders.

If the user’s input sources are from Skelton et al. (2014), there is no
need to provide object coordinates; ``MCSED`` will retrieve the ICRS
coordinates directly from the catalog using the objects’ field names and
ID numbers.

Given the coordinates of each source, ``MCSED`` first generates a list
of the differential extinctions for all input galaxies. These
:math:`E(B-V)` values are converted to optical depths
:math:`\tau_\lambda=A_\lambda/1.086` where :math:`A_{\lambda} / A_V`
comes from the expressions given by Fitzpatrick (1999) and generated by
the Python package
`dust_extinction <https://dust-extinction.readthedocs.io/en/latest/>`__.

The user needs to install
`dustmaps <https://dustmaps.readthedocs.io/en/latest/installation.html>`__
only if they opt for letting ``MCSED`` perform ISM corrections. One of
the components required by ``dustmaps``, ``healpy``, is not available on
Windows (but please let us know if you are able to install it!), so currently, the ISM corrections will not work on Windows machines.

The Schlegel, Finkbeiner, & Davis (1998) 2-D dust maps must be downloaded and properly configured for the Milky Way corrections in ``MCSED`` to work. See :ref:`subsec:requirements` for details.