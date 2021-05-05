.. _sec:parameterlist:

Summary of Important Variables
==============================

+---------------+------------------+-------------+-------------+-------------+
| Variable      | Location         | Default     | Range       | Description |
+===============+==================+=============+=============+=============+
| *Cosmology Parameters*           |             |             |             |
+---------------+------------------+-------------+-------------+-------------+
| H_0           | cosmology.py     | 69          | …           | Hubble Const|
|               |                  |             |             | (km/s/Mpc)  |
+---------------+------------------+-------------+-------------+-------------+
| omega_m       | cosmology.py     | 0.31        | …           | Matter      |
|               |                  |             |             | Density     |
+---------------+------------------+-------------+-------------+-------------+
| omega_l       | cosmology.py     | 0.69        | …           | Dark Energy |
|               |                  |             |             | Density     |
+---------------+------------------+-------------+-------------+-------------+
| *Input Data Parameters*          |             |             |             |
+---------------+------------------+-------------+-------------+-------------+
| use_input     | config.py        | True        | …           | Do not      |
| _data         |                  |             |             | ignore      |
|               |                  |             |             | additional  |
|               |                  |             |             | data in the |
|               |                  |             |             | input file  |
+---------------+------------------+-------------+-------------+-------------+
| phot_floor    | config.py        | 0.05        | …           | Minimum     |
| _error        |                  |             |             | fractional  |
|               |                  |             |             | error for   |
|               |                  |             |             | photometry  |
+---------------+------------------+-------------+-------------+-------------+
| emline_floor  | config.py        | 0.05        | …           | Minimum     |
| _error        |                  |             |             | fractional  |
|               |                  |             |             | error for   |
|               |                  |             |             | emission    |
|               |                  |             |             | line fluxes |
+---------------+------------------+-------------+-------------+-------------+
| absindx_floor | config.py        | 0.05        | …           | Minimum     |
| _error        |                  |             |             | fractional  |
|               |                  |             |             | error for   |
|               |                  |             |             | absorption  |
|               |                  |             |             | line        |
|               |                  |             |             | indices     |
+---------------+------------------+-------------+-------------+-------------+
| *SSP Parameters*                 |             |             |             |
+---------------+------------------+-------------+-------------+-------------+
| ssp           | config.py        | fsps        | …           | SSP Spectra |
+---------------+------------------+-------------+-------------+-------------+
| isochrone     | config.py        | padova      | …           |SSP Isochrone|
+---------------+------------------+-------------+-------------+-------------+
| metallicity   | config.py        | 0.0077      | [0.0001,    | Galaxy      |
|               |                  |             | 0.04]       | Metallicity |
+---------------+------------------+-------------+-------------+-------------+
| *Dust Attenuation Parameters*    |             |             |             |
+---------------+------------------+-------------+-------------+-------------+
| dust_law      | config.py        | calzetti    | …           | Attenuation |
|               |                  |             |             | Law (1 of 5 |
|               |                  |             |             | options)    |
+---------------+------------------+-------------+-------------+-------------+
| EBV_old_young | config.py        | 0.44        | …           | Extinction  |
|               |                  |             |             | of diffuse  |
|               |                  |             |             | component   |
|               |                  |             |             | compared to |
|               |                  |             |             | that of     |
|               |                  |             |             | birth clouds|
|               |                  |             |             | :math:`\eta`|
+---------------+------------------+-------------+-------------+-------------+
| t_birth       | config.py        | 7           | …           | Attenuation |
|               |                  |             |             | age divider |
|               |                  |             |             | (log yr)    |
+---------------+------------------+-------------+-------------+-------------+
| EBV           | dust_abs.py      | 0.15        | [-0.5,      | Reddening   |
|               |                  |             | 1.5]        |             |
+---------------+------------------+-------------+-------------+-------------+
| Rv            | config.py        | :math:`-1`  | …           | Specific to |
|               |                  |             |             | attenuation |
|               |                  |             |             | law         |
+---------------+------------------+-------------+-------------+-------------+
| delta         | dust_abs.py      | 0.0         | [-1.0,+1.0] | UV slope of |
|               |                  |             |             | Noll        |
|               |                  |             |             | attenuation |
|               |                  |             |             | law         |
+---------------+------------------+-------------+-------------+-------------+
| Eb            | dust_abs.py      | 2.5         | [-0.2,6.0]  | Bump        |
|               |                  |             |             | strength    |
|               |                  |             |             | for Noll    |
|               |                  |             |             | attenuation |
|               |                  |             |             | law         |
+---------------+------------------+-------------+-------------+-------------+
| *Dust Emission Parameters*       |             |             |             |
+---------------+------------------+-------------+-------------+-------------+
| dust_em       | config.py        | DL07        | …           | Prescription|
|               |                  |             |             | for dust    |
|               |                  |             |             | emission (1 |
|               |                  |             |             | of 1        |
|               |                  |             |             | option)     |
+---------------+------------------+-------------+-------------+-------------+
| wave_dust_em  | config.py        | 2.5         | …           | Maximum     |
|               |                  |             |             | wavelength  |
|               |                  |             |             | to fit      |
|               |                  |             |             | (microns)   |
+---------------+------------------+-------------+-------------+-------------+
| assume_energy | config.py        | False       | [T,F]       | Boolean for |
| _balance      |                  |             |             | balancing   |
|               |                  |             |             | dust        |
|               |                  |             |             | attenuation |
|               |                  |             |             | with        |
|               |                  |             |             | emission    |
+---------------+------------------+-------------+-------------+-------------+
| umin          | dust_emission.py | 2.0         | [0.1, 25.]  | `U_min`     |
|               |                  |             |             | for the     |
|               |                  |             |             | Draine & Li |
|               |                  |             |             | (2007) model|
+---------------+------------------+-------------+-------------+-------------+
| gamma         | dust_emission.py | 0.05        | [0.0, 1.0]  | γ           |
|               |                  |             |             | for the     |
|               |                  |             |             | Draine & Li |
|               |                  |             |             | (2007)      |
|               |                  |             |             | model       |
+---------------+------------------+-------------+-------------+-------------+
| qpah          | dust_emission.py | 2.5         | [0.47, 4.58]| `q_PAH` for |
|               |                  |             |             | the Draine  |
|               |                  |             |             | & Li (2007) |
|               |                  |             |             | model       |
+---------------+------------------+-------------+-------------+-------------+
| mdust         | dust_emission.py | 7.0         | [4.5, 10.0] | Log Dust    |
|               |                  |             |             | Mass (fit if|
|               |                  |             |             | energy      |
|               |                  |             |             | balance     |
|               |                  |             |             | not assumed)|
+---------------+------------------+-------------+-------------+-------------+
| *Star Formation History Parameters*            |             |             |
+---------------+------------------+-------------+-------------+-------------+
| sfh           | config.py        | binned_lsfr | …           | Star        |
|               |                  |             |             | Formation   |
|               |                  |             |             | Rate        |
|               |                  |             |             | History (1  |
|               |                  |             |             | of 7        |
|               |                  |             |             | options)    |
+---------------+------------------+-------------+-------------+-------------+
| age           | sfh.py           | :math:`-0.5`| [-3,…]      | Log Age of  |
|               |                  |             |             | the galaxy  |
|               |                  |             |             | (Gyr)       |
+---------------+------------------+-------------+-------------+-------------+
| logsfr        | sfh.py           | 1.0         | [-3.0, +3.0]| Current SFR |
|               |                  |             |             | (log solar  |
|               |                  |             |             | masses per  |
|               |                  |             |             | year)       |
+---------------+------------------+-------------+-------------+-------------+
| tau           | sfh.py           | :math:`-1.5`| [-3.5, +3.5]| E-folding   |
|               |                  |             |             | timescale   |
|               |                  |             |             | for         |
|               |                  |             |             | exponential |
|               |                  |             |             | SFH         |
+---------------+------------------+-------------+-------------+-------------+
| b             | sfh.py           | 2.0         | [0.0, +5.0] | Increasing  |
|               |                  |             |             | power-law   |
|               |                  |             |             | index for   |
|               |                  |             |             | Behroozi    |
|               |                  |             |             | et al. SFH  |
+---------------+------------------+-------------+-------------+-------------+
| c             | sfh.py           | 1.0         | [0.0, +5.0] | Decreasing  |
|               |                  |             |             | power-law   |
|               |                  |             |             | index for   |
|               |                  |             |             | Behroozi    |
|               |                  |             |             | et al. SFH  |
+---------------+------------------+-------------+-------------+-------------+
| taup          | sfh.py           | 2.4         | [-3.0, +1.0]| Timescale   |
|               |                  |             |             | for         |
|               |                  |             |             | Behroozi    |
|               |                  |             |             | et al. SFH  |
|               |                  |             |             | (Gyr)       |
+---------------+------------------+-------------+-------------+-------------+
| burst_age     | sfh.py           | 7.2         | [6.0, 7.5]  | Time since  |
|               |                  |             |             | the burst   |
|               |                  |             |             | (log yr)    |
+---------------+------------------+-------------+-------------+-------------+
| burst_sigma   | sfh.py           | 0.4         | [0.003,     | Log-normal  |
|               |                  |             | 0.5]        | burst       |
|               |                  |             |             | duration    |
|               |                  |             |             | (log Gyr)   |
+---------------+------------------+-------------+-------------+-------------+
| burst_strength| sfh.py           | 5.0         | [1.0, 10.0] | Burst       |
|               |                  |             |             | strength    |
|               |                  |             |             | relative to |
|               |                  |             |             | mean SFR    |
+---------------+------------------+-------------+-------------+-------------+
| age_locs      | sfh.py           | 6.5, 7.5,   | …           | Array of    |
|               |                  | 8.5         |             | pivot       |
|               |                  |             |             | points for  |
|               |                  |             |             | polynomial  |
|               |                  |             |             | SFH (log    |
|               |                  |             |             | yr)         |
+---------------+------------------+-------------+-------------+-------------+
| ages          | sfh.py           | 8.0, 8.5,   | …           | Age bins    |
|               |                  | 9.0, 9.5,   |             | limits for  |
|               |                  | 9.8, 10.12  |             | binned SFH  |
|               |                  |             |             | (log yr)    |
+---------------+------------------+-------------+-------------+-------------+
| *ISM and IGM Parameters*         |             |             |             |
+---------------+------------------+-------------+-------------+-------------+
| add_nebular   | config.py        | True        | …           | Boolean for |
|               |                  |             |             | including   |
|               |                  |             |             | nebular     |
|               |                  |             |             | emission    |
+---------------+------------------+-------------+-------------+-------------+
| logU          | config.py        | :math:`-2.5`| …           | Ionization  |
|               |                  |             |             | Parameter   |
+---------------+------------------+-------------+-------------+-------------+
| emline_factor | config.py        |1.0e-17      | …           | Emission    |
|               |                  |             |             | line        |
|               |                  |             |             | conversion  |
|               |                  |             |             | to          |
|               |                  |             |             | ergs/cm²/s  |
+---------------+------------------+-------------+-------------+-------------+
| ISM_correct   | config.py        | None        | Valid       | Coordinate  |
| _coords       |                  |             | Coord.      | system      |
|               |                  |             | System      | (string) of |
|               |                  |             |             | sources for |
|               |                  |             |             | Milky Way   |
|               |                  |             |             | dust        |
|               |                  |             |             | correction  |
+---------------+------------------+-------------+-------------+-------------+
| IGM_correct   | config.py        | False       | [T,F]       | Boolean for |
|               |                  |             |             | applying    |
|               |                  |             |             | correlation |
|               |                  |             |             | for IGM     |
|               |                  |             |             | absorption  |
+---------------+------------------+-------------+-------------+-------------+
| blue_wave     | config.py        | 1216        | …           | Minimum     |
| _cutoff       |                  |             |             | wavelength  |
|               |                  |             |             | to fit (Å)  |
+---------------+------------------+-------------+-------------+-------------+
| *Calculation Parameters*         |             |             |             |
+---------------+------------------+-------------+-------------+-------------+
| sigma_m       | mcsed.py         | 0.1         | …           | Error term  |
|               |                  |             |             | associated  |
|               |                  |             |             | with the    |
|               |                  |             |             | models      |
+---------------+------------------+-------------+-------------+-------------+
| nwalkers      | config.py        | 100         | …           | Number of   |
|               |                  |             |             | MCMC        |
|               |                  |             |             | walkers     |
+---------------+------------------+-------------+-------------+-------------+
| nsteps        | config.py        | 1000        | …           | Number of   |
|               |                  |             |             | steps for   |
|               |                  |             |             | each walker |
+---------------+------------------+-------------+-------------+-------------+
| progress_bar  | config.py        | False       | …           | Show bar    |
|               |                  |             |             | of fit      |
|               |                  |             |             | progress    |
+---------------+------------------+-------------+-------------+-------------+
| force_emcee   | config.py        | True        | …           | Force emcee |
| _finish       |                  |             |             | to finish   |
|               |                  |             |             | even if no  |
|               |                  |             |             | convergence |
+---------------+------------------+-------------+-------------+-------------+
| burnin        | config.py        | 0.25        | …           | Fraction of |
| _fraction     |                  |             |             | nsteps to   |
|               |                  |             |             | count as    |
|               |                  |             |             | burnin, if  |
|               |                  |             |             | force_      |
|               |                  |             |             | finish=True |
+---------------+------------------+-------------+-------------+-------------+
| reserved      | config.py        | 2           | …           | In          |
| _cores        |                  |             |             | parallel, # |
|               |                  |             |             | cores =     |
|               |                  |             |             | Total cores |
|               |                  |             |             | - Reserved  |
|               |                  |             |             | cores       |
+---------------+------------------+-------------+-------------+-------------+
| nobjects      | config.py        | 5           | …           | In test     |
|               |                  |             |             | mode,       |
|               |                  |             |             | number of   |
|               |                  |             |             | objects to  |
|               |                  |             |             | analyze     |
+---------------+------------------+-------------+-------------+-------------+
| test_zrange   | config.py        | …           | [1.0,       | Range of    |
|               |                  |             | 2.0]        | redshifts   |
|               |                  |             |             | for test    |
|               |                  |             |             | objects     |
+---------------+------------------+-------------+-------------+-------------+
| param         | config.py        | 5, 16, 50,  | …           | % of each   |
| _percentiles  |                  | 84, 95      |             | parameter   |
|               |                  |             |             | to report   |
+---------------+------------------+-------------+-------------+-------------+
| separate_     | config.py        | False       | …           | Return the  |
| stars_gas     |                  |             |             | stellar and |
|               |                  |             |             | nebular     |
|               |                  |             |             | components  |
+---------------+------------------+-------------+-------------+-------------+
| *Output Parameters*              |             |             |             |
+---------------+------------------+-------------+-------------+-------------+
| output_dict:  | config.py        | True        | [T,F]       | Produce     |
| parameters    |                  |             |             | Summary     |
|               |                  |             |             | File        |
+---------------+------------------+-------------+-------------+-------------+
| output_dict:  | config.py        | True        | [T,F]       | List        |
| settings      |                  |             |             | user-defined|
|               |                  |             |             | fitting     |
|               |                  |             |             | assumptions |
+---------------+------------------+-------------+-------------+-------------+
| output_dict:  | config.py        | False       | [T,F]       | Save all    |
| fitposterior  |                  |             |             | posterior   |
|               |                  |             |             | probability |
|               |                  |             |             | distribs    |
+---------------+------------------+-------------+-------------+-------------+
| output_dict:  | config.py        | False       | [T,F]       | Save each   |
| fullposterior |                  |             |             | walker,     |
|               |                  |             |             | step, and   |
|               |                  |             |             | parameter   |
|               |                  |             |             | information |
+---------------+------------------+-------------+-------------+-------------+
| output_dict:  | config.py        | True        | [T,F]       | Save        |
| bestfitspec   |                  |             |             | best-fitting|
|               |                  |             |             | SED         |
+---------------+------------------+-------------+-------------+-------------+
| output_dict:  | config.py        | True        | [T,F]       | Save        |
| fluxdensity   |                  |             |             | comparison  |
|               |                  |             |             | of observed |
|               |                  |             |             | and modeled |
|               |                  |             |             | photometry  |
+---------------+------------------+-------------+-------------+-------------+
| output_dict:  | config.py        | True        | [T,F]       | Save        |
| lineflux      |                  |             |             | comparison  |
|               |                  |             |             | of observed |
|               |                  |             |             | and modeled |
|               |                  |             |             | emission    |
|               |                  |             |             | lines       |
+---------------+------------------+-------------+-------------+-------------+
| output_dict:  | config.py        | True        | [T,F]       | Save        |
| absindx       |                  |             |             | comparison  |
|               |                  |             |             | of observed |
|               |                  |             |             | and modeled |
|               |                  |             |             | absorption  |
|               |                  |             |             | indices     |
+---------------+------------------+-------------+-------------+-------------+
| output_dict:  | config.py        | True        | [T,F]       | Produce     |
| triangle plot |                  |             |             | summary     |
|               |                  |             |             | figure for  |
|               |                  |             |             | each object |
+---------------+------------------+-------------+-------------+-------------+
| output_dict:  | config.py        | False       | [T,F]       | Save figure |
| sample plot   |                  |             |             | of parameter|
|               |                  |             |             | estimates   |
|               |                  |             |             | for all     |
|               |                  |             |             | MCMC chains |
+---------------+------------------+-------------+-------------+-------------+
| output_dict:  | config.py        | True        | [T,F]       | Save figure |
| template spec |                  |             |             | of the      |
|               |                  |             |             | age-weighted|
|               |                  |             |             | SSP spectra |
+---------------+------------------+-------------+-------------+-------------+
| output_dict:  | config.py        | .png        | …           | Image       |
| image format  |                  |             |             | format (1 of|
|               |                  |             |             | 9 options)  |
+---------------+------------------+-------------+-------------+-------------+
