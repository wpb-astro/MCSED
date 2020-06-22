.. _section-install:

Installation
============

To install ``MCSED``, move to a directory where you would like it to be
stored and type:

``git clone https://github.com/wpb-astro/MCSED.git``

A directory called ``MCSED`` will be created containing all of the
necessary files for the program.

.. _subsec:requirements:

Requirements
------------

``MCSED`` is a python 2.7 based code and requires a few standard python
based packages. All of these packages can be found in the Anaconda
distribution environment. Specific versions of some packages must be
specified:

 ``emcee`` version 2.1.0

 ``corner`` version 2.0.1

 ``seaborn`` version 0.8.1

 ``astropy`` version 2.0.6

 ``matplotlib`` version 2.1.2

 ``scipy`` version 1.0.0

 ``dustmaps`` (if a correction for foreground Milky Way dust extinction
is desired)

If a foreground Milky Way dust extinction correction is desired, one
must install the Schlegel, Finkbeiner, & Davis (1998) 2-D dust maps, as
described in the `dustmaps documentation <https://dustmaps.readthedocs.io/en/latest/installation.html>`__.
After installing ``dustmaps``, make a directory called ``sfd`` in the
location of the ``MCSED`` directory and then begin an interactive python
session by issuing the command ``ipython``, ``python -i``, or ``python`` in a
terminal session. Once the interactive session begins, issue:

``>> from dustmaps.config import config``

``>> config['data_dir'] = <MCSED_LOCATION>/sfd``

``>> import dustmaps.sfd``

``>> dustmaps.sfd.fetch()``

where ``<MCSED_LOCATION>`` is the location of the ``MCSED`` directory
that is created by the ``git clone`` call in :ref:`section-install`.

**Note**: One of the components required by ``dustmaps``, ``healpy``, is not available on Windows, so currently,
the ISM corrections will not work on Windows machines.

