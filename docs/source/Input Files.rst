.. _section:inputs:

Input Files
===========

``MCSED`` is designed to model the SEDs of a set of galaxies with a
common set of photometric and spectroscopic data. Since ``MCSED`` can
accept a wide range of constraints, with dozens of possible filter
combinations, emission line fluxes, and/or absorption line spectra
indices, the format of the input file has some flexibility. But the
basic file structure is simple: the data are entered in a simple,
space-delimited ascii file, with the first line of the file labeling the
file’s columns. The rest of the file gives the measured flux densities
and/or emission-line fluxes and/or absorption line spectral indices for
each object, one object per line.

.. _subsec:columns:

Required Columns
----------------

The input data file has three required columns, which must have these
exact labels: ``Field``, ``ID``, and ``z``. In other words, an
object’s identification consists in two parts: a string which contains
the name of the field in which the object is found, and a unique integer
ID which is specific to that field. If both field and ID are
unnecessary, one can simply enter a placeholder for one of the entries.
Redshifts must be specified for every source.

The remaining columns in the input file should be pairs of numbers
representing photometric flux densities (and their :math:`1\,\sigma`
errors), emission line fluxes (and their errors), and/or absorption line
indices (and their errors). Since the quoted uncertainties for
photometric observations often do not include systematic and/or external
errors, ``MCSED`` also allows the user to specify a minimum fractional
uncertainty for any type of observation. The defaults for the minimum
errors can be found in ``config.py``, and by default are set to
``phot_floor_error`` = 0.05 (for photometric errors),
``emline_floor_error`` = 0.05 (for errors in emission-line fluxes), and
``absindx_floor_error`` = 0.05 (for errors in absorption line indices).
These defaults can be changed by editing the above parameters in
``config.py``.

.. _subsec:photometry:

Photometry
----------

.. _subsubsec:skelton:

Using Skelton et al. (2014)
~~~~~~~~~~~~~~~~~~~~~~~~~~~

``MCSED`` was originally written to analyze galaxies in the five CANDELS
fields (AEGIS, COSMOS, GOODS-N, GOODS-S, and UDS), hence there are
special commands built into the program to handle the PSF-matched
photometry from Skelton et al. (2014). If the user’s sources are in the
Skelton catalog, the objects can be specified by their field (i.e.,
``aegis``, ``cosmos``, ``goodsn``, ``goodss``, or ``uds``) and the
unique Skelton ID number. This links the input line directly to the
object’s photometry in the files provided by Skelton et al. (2014).
Momcheva et al. (2016) provide grism redshifts (and emission line
fluxes) for all Skelton photometry. Users with Skelton sources are
encouraged to use the grism redshifts for the redshift column.

Additional photometry not included in the Skelton catalogs can be
specified in the input file in the same way as general photometry as
discussed in :ref:`subsubsec:genphot`.

.. _subsubsec:genphot:

General Case
~~~~~~~~~~~~

If the input objects are not associated with the Skelton et al. (2014) catalog
(identified via the ``Field`` and ``ID`` columns described above), or if users
wish to supplement this catalog with additional photometry, the input file must
include additional columns. Photometric measurements should be given as flux
densities with :math:`1\,\sigma` uncertainties associated with each
measurement. The columns containing these data in the input file should be labeled
``f_filter_name`` and ``e_filter_name``, where ``filter_name`` is the
name of a ``.res`` file in the ``FILTERS`` directory. (In other words,
columns named ``f_hst_acsF606W`` and ``e_hst_acsF606W`` should refer to
the flux densities (not magnitudes!) and uncertainties taken through the
filter defined in ``FILTERS/hst_acsF606W.res``.) Following Skelton
et al. (2014), the units for flux density are scaled to an AB magnitude
of 25, so :math:`1.00 = 3.63 \times 10^{-30}` ergs cm\ :math:`^{-2}` s\ :math:`^{-1}` Hz\ :math:`^{-1}` (e.g., if the user’s flux densities are in :math:`\mu`\ Jy, the values must be multiplied by :math:`10^{0.4(25-23.9)} \approx 2.754`).

.. _subsec:emission-lines:

Emission Lines
--------------

``MCSED`` can include emission line fluxes in the likelihood function.
To do this, the user first specifies the line’s name (keyword ``Name``),
rest-frame wavelength (in Angstroms), and relative weight in the
``config.py`` emission-line dictionary. A weight of 1.0 means the line
contributes just as much weight to the likelihood function as a
photometric data point; a weight of 0.0 implies that the line is
ignored. The user then provides the objects’ emission line strengths and
:math:`1\,\sigma` error bars in exactly the same way as for the
photometry, i.e., by entering the data in the input file and labeling
the columns as ``Name_FLUX`` and ``Name_ERR``, where ``Name`` is the
line’s keyword as defined in ``config.py``. The emission line fluxes and
errors must be specified in units of :math:`10^{-17}` erg
cm\ :math:`^{-2}` s\ :math:`^{-1}`, unless a different multiplication
factor to the base unit of erg cm\ :math:`^{-2}` s\ :math:`^{-1}` is
specified by the keyword ``emline_factor`` in ``config.py``. The
emission lines currently included in ``config.py`` are given below.
Additional lines can be added by expanding the ``emline_list_dict`` in
``config.py``.

.. table:: Emission Lines Definitions

   +------------------------+----------+------------+--------+
   |  Line                  | Name     | Wavelength | Weight | 
   |                        |          | (Å)        |        |        
   +========================+==========+============+========+
   | H\ :math:`\beta`       | ``Hb``   | 4861       | 1.0    |
   +------------------------+----------+------------+--------+
   | H\ :math:`\alpha`      | ``Ha``   | 6563       | 1.0    |
   +------------------------+----------+------------+--------+
   | [O III]                | ``OIII`` | 5007       | 0.5    |
   +------------------------+----------+------------+--------+
   | [O II]                 | ``OII``  | 3727       | 0.5    |
   +------------------------+----------+------------+--------+
   | [N II]                 | ``NII``  | 6583       | 0.5    |
   +------------------------+----------+------------+--------+

Currently, ``MCSED`` can not fit blended emission lines.

.. _subsec:absorption-lines:

Absorption Line Indices
-----------------------

Absorption line indices can also be used in ``MCSED``’s likelihood
function. The way to do this is similar to that used for emission lines.
In the input file, the columns containing an absorption line index and
its uncertainty are labeled as ``Name_INDX`` and ``Name_Err``, where
``Name`` is the line’s keyword, as defined in ``config.py``. The indices
that are pre-defined in ``MCSED`` are listed in the table below. As one can see from the table,
the indices are defined via their wavelength ranges, the units they are
quoted in, and a relative weight similar to that defined for the
emission lines.

.. table:: Absorption Line Indices Definitions

   +-------+-------+--------+-------+------------+-------+-----------+-------+-------+
   |       | Index Band (Å) | Blue Continuum (Å) | Red Continuum (Å) |               |
   +=======+=======+========+=======+============+=======+===========+=======+=======+
   | Name  | Blue  | Red    | Blue  | Red        | Blue  | Red       | Units¹| Weight|
   +-------+-------+--------+-------+------------+-------+-----------+-------+-------+
   | Lick  | 4142. | 4177.  | 4080. | 4117.      | 4244. | 4284.     | 1     | 1.0   |
   | CN_1  | 125   | 125    | 125   | 625        | 125   | 125       |       |       |
   +-------+-------+--------+-------+------------+-------+-----------+-------+-------+
   | Lick  | 4142. | 4177.  | 4083. | 4096.      | 4244. | 4284.     | 1     | 1.0   |
   | CN_2  | 125   | 125    | 875   | 375        | 125   | 125       |       |       |
   +-------+-------+--------+-------+------------+-------+-----------+-------+-------+
   | Lick  | 4222. | 4234.  | 4211. | 4219.      | 4241. | 4251.     | 0     | 1.0   |
   | Ca4227| 250   | 750    | 000   | 750        | 000   | 000       |       |       |
   +-------+-------+--------+-------+------------+-------+-----------+-------+-------+
   | Lick  | 4281. | 4316.  | 4266. | 4282.      | 4318. | 4335.     | 0     | 1.0   |
   | G4300 | 375   | 375    | 375   | 625        | 875   | 125       |       |       |
   +-------+-------+--------+-------+------------+-------+-----------+-------+-------+
   | Lick  | 4369. | 4420.  | 4359. | 4370.      | 4442. | 4455.     | 0     | 1.0   |
   | Fe4383| 125   | 375    | 125   | 375        | 875   | 375       |       |       |
   +-------+-------+--------+-------+------------+-------+-----------+-------+-------+
   | Lick  | 4452. | 4474.  | 4445. | 4454.      | 4477. | 4492.     | 0     | 1.0   |
   | Ca4455| 125   | 625    | 875   | 625        | 125   | 125       |       |       |
   +-------+-------+--------+-------+------------+-------+-----------+-------+-------+
   | Lick  | 4514. | 4559.  | 4504. | 4514.      | 4560. | 4579.     | 0     | 1.0   |
   | Fe4531| 250   | 250    | 250   | 250        | 500   | 250       |       |       |
   +-------+-------+--------+-------+------------+-------+-----------+-------+-------+
   | Lick  | 4634. | 4720.  | 4611. | 4630.      | 4742. | 4756.     | 0     | 1.0   |
   | Fe4668| 000   | 250    | 500   | 250        | 750   | 500       |       |       |
   +-------+-------+--------+-------+------------+-------+-----------+-------+-------+
   | Lick  | 4847. | 4876.  | 4827. | 4847.      | 4876. | 4891.     | 0     | 1.0   |
   | Hb    | 875   | 625    | 875   | 875        | 625   | 625       |       |       |
   +-------+-------+--------+-------+------------+-------+-----------+-------+-------+
   | Lick  | 4977. | 5054.  | 4946. | 4977.      | 5054. | 5065.     | 0     | 1.0   |
   | Fe5015| 750   | 000    | 500   | 750        | 000   | 250       |       |       |
   +-------+-------+--------+-------+------------+-------+-----------+-------+-------+
   | Lick  | 5069. | 5134.  | 4895. | 4957.      | 5301. | 5366.     | 1     | 1.0   |
   | Mg1   | 125   | 125    | 125   | 625        | 125   | 125       |       |       |
   +-------+-------+--------+-------+------------+-------+-----------+-------+-------+
   | Lick  | 5154. | 5196.  | 4895. | 4957.      | 5301. | 5366.     | 1     | 1.0   |
   | Mg2   | 125   | 625    | 125   | 625        | 125   | 125       |       |       |
   +-------+-------+--------+-------+------------+-------+-----------+-------+-------+
   | Lick  | 5160. | 5192.  | 5142. | 5161.      | 5191. | 5206.     | 0     | 1.0   |
   | Mgb   | 125   | 625    | 625   | 375        | 375   | 375       |       |       |
   +-------+-------+--------+-------+------------+-------+-----------+-------+-------+
   | Lick  | 5245. | 5285.  | 5233. | 5248.      | 5285. | 5318.     | 0     | 1.0   |
   | Fe5270| 650   | 650    | 150   | 150        | 650   | 150       |       |       |
   +-------+-------+--------+-------+------------+-------+-----------+-------+-------+
   | Lick  | 5312. | 5352.  | 5304. | 5315.      | 5353. | 5363.     | 0     | 1.0   |
   | Fe5335| 125   | 125    | 625   | 875        | 375   | 375       |       |       |
   +-------+-------+--------+-------+------------+-------+-----------+-------+-------+
   | Lick  | 5387. | 5415.  | 5376. | 5387.      | 5415. | 5425.     | 0     | 1.0   |
   | Fe5406| 500   | 000    | 250   | 500        | 000   | 000       |       |       |
   +-------+-------+--------+-------+------------+-------+-----------+-------+-------+
   | Lick  | 5696. | 5720.  | 5672. | 5696.      | 5722. | 5736.     | 0     | 1.0   |
   | Fe5709| 625   | 375    | 875   | 625        | 875   | 625       |       |       |
   +-------+-------+--------+-------+------------+-------+-----------+-------+-------+
   | Lick  | 5776. | 5796.  | 5765. | 5775.      | 5797. | 5811.     | 0     | 1.0   |
   | Fe5782| 625   | 625    | 375   | 375        | 875   | 625       |       |       |
   +-------+-------+--------+-------+------------+-------+-----------+-------+-------+
   | Lick  | 5876. | 5909.  | 5860. | 5875.      | 5922. | 5948.     | 0     | 1.0   |
   | NaD   | 875   | 375    | 625   | 625        | 125   | 125       |       |       |
   +-------+-------+--------+-------+------------+-------+-----------+-------+-------+
   | Lick  | 5936. | 5994.  | 5816. | 5849.      | 6038. | 6103.     | 1     | 1.0   |
   | TiO1  | 625   | 125    | 625   | 125        | 625   | 625       |       |       |
   +-------+-------+--------+-------+------------+-------+-----------+-------+-------+
   | Lick  | 6189. | 6272.  | 6066. | 6141.      | 6372. | 6415.     | 1     | 1.0   |
   | TiO2  | 625   | 125    | 625   | 625        | 625   | 125       |       |       |
   +-------+-------+--------+-------+------------+-------+-----------+-------+-------+
   | Lick  | 4083. | 4122.  | 4041. | 4079.      | 4128. | 4161.     | 0     | 1.0   |
   | Hd_A  | 500   | 250    | 600   | 750        | 500   | 000       |       |       |
   +-------+-------+--------+-------+------------+-------+-----------+-------+-------+
   | Lick  | 4319. | 4363.  | 4283. | 4319.      | 4367. | 4419.     | 0     | 1.0   |
   | Hg_A  | 750   | 500    | 500   | 750        | 250   | 750       |       |       |
   +-------+-------+--------+-------+------------+-------+-----------+-------+-------+
   | Lick  | 4091. | 4112.  | 4057. | 4088.      | 4114. | 4137.     | 0     | 1.0   |
   | Hd_F  | 000   | 250    | 250   | 500        | 750   | 250       |       |       |
   +-------+-------+--------+-------+------------+-------+-----------+-------+-------+
   | Lick  | 4331. | 4352.  | 4283. | 4319.      | 4354. | 4384.     | 0     | 1.0   |
   | Hg_F  | 250   | 250    | 500   | 750        | 750   | 750       |       |       |
   +-------+-------+--------+-------+------------+-------+-----------+-------+-------+
   | D4000 | …     | …      | 3750. | 3950.      | 4050. | 4250.     | 2     | 1.0   |
   |       |       |        | 000   | 000        | 000   | 000       |       |       |
   +-------+-------+--------+-------+------------+-------+-----------+-------+-------+
   |¹Unit codes: 0 = Å; 1 = mag; 2 = ratio                                           |
   +-------+-------+--------+-------+------------+-------+-----------+-------+-------+

These definitions come from Bruzual (1983) and Worthey et al. (1994);
they are calculated by finding the average value of :math:`F_{\lambda}`
within the blue and red continuum bands, interpolating a line through
these values to estimate the continuum, :math:`F_C`, and then computing
equivalent width via

.. math:: {\rm EW} = \int_{\lambda_1}^{\lambda_2} \left( 1 - \frac{F_{\lambda}}{F_C} \right) d\lambda

**Important Note:** absorption line indices are defined for a specific
spectral resolution. ``MCSED`` makes no attempt to match this
resolution: it uses the SSP spectra as is. The user should consider this
carefully before deciding on the utility of this feature.
