""" MCSED - config.py

1) Configuration Settings

.. moduleauthor:: Greg Zeimann <gregz@astro.as.utexas.edu>

"""
# SSP code for models
ssp = 'fsps'           # options include: 'fsps'
isochrone = 'padova'   # options include: 'padova'
# SFH options include: 'constant', 'burst', 'polynomial', 'exponential', 
#                      'double_powerlaw', 'binned_lsfr'
sfh = 'constant' 
dust_law = 'calzetti'  # options include: 'calzetti', 'noll', 'reddy', 
                       #                  'conroy', 'cardelli'

# Dust emission parameters: 
#   if False, do not fit for dust emission component and remove all filters
#             redward of rest-frame wave_dust_em microns  (defined below)
#   else, set to string of desired dust emission class 
dust_em = False # options include: 'DL07', False 
# Assume energy balance or normalize the dust IR spectrum as a free parameter
assume_energy_balance = False

# Dust attenuation law parameters
#   extinction factor (if negative, use default value for dust law of choice)
Rv = -1 
# The relative attenuation between the birth cloud and the diffuse component
#   in the dust model, such that 
#   E(B-V)_diffuse = EBV_old_young * E(B-V)_birth
EBV_old_young = 0.44

t_birth = 7. # age of the birth cloud (log years)

# Ignore photometry, as appropriate
#   blue_wave_cutoff: ignore filters containing Lyman-alpha
#   wave_dust_em:     if not fitting dust emission component, ignore photometry
#                     dominated by dust emission
blue_wave_cutoff = 1216. # rest-frame wavelength in Angstroms 
wave_dust_em     = 30.   # rest-frame wavelength in microns 

# Stellar metallicity
#   If False, leave metallicity as a free model parameter
#   else, must be float: fixed metallicity of SSP models (Z_solar = 0.019)
metallicity = False

# Nebular Emission Properties
# The ionization parameter, logU, is held fixed
logU = -2.5

# EMCEE parameters
nwalkers = 100 
nsteps   = 1000 

# Number of test objects
nobjects = 5
test_zrange = (1.0, 2.0) # redshift range of test objects (uniform prior)

# Minimum fractional errors in observed photometry, 
#   emission line fluxes, and absorption line indices 
phot_floor_error    = 0.05
emline_floor_error  = 0.05
absindx_floor_error = 0.05

# Fractional error expected from the models, i.e., fractional error adopted
#   for model photometry, emission line fluxes, and absorption line indices
model_floor_error = 0.10

# Use input data (photometry, emission lines, absorption line indices)
#   If True, use additional data provided in the input file
#   else, ignore input data (in which case input Field,ID must match Skelton+14)
use_input_data = True 

# ISM/IGM correction
ISM_correct_coords = None # if None, do not apply an ISM correction
# Options for coords: 'altaz', 'barycentrictrueecliptic', 'cirs', 'fk4', 
#                     'fk4noeterms', 'fk5', 'galactic', 'galacticlsr', 
#                     'galactocentric', 'gcrs', 'geocentrictrueecliptic', 
#                     'hcrs', 'heliocentrictrueecliptic', 'icrs', 'itrs', 
#                     'lsr', 'precessedgeocentric', 'supergalactic'
IGM_correct = False

# Separate the stellar/nebular components
#   slower by factor of ~ 2, only needed if wish to return
#   best-fit spectrum for stellar and nebular components separately
separate_stars_gas = False 


# Output files
#   Supported image formats: eps, pdf, pgf, png, ps, raw, rgba, svg, svgz
output_dict = {'parameters'    : True,   # fitted parameters
               'settings'      : True,   # user-defined fitting assumptions
               'fitposterior'  : False,  # parameter posterior distributions
               'bestfitspec'   : True,   # best-fit SED model
               'fluxdensity'   : True,   # modeled and observed photometry
               'lineflux'      : True,   # modeled and observed emission lines
               'absindx'       : True,   # modeled, observed absorption indices
               'triangle plot' : True,   # summary diagnostic plot
               'sample plot'   : False,  # parameter estimates for MCMC chains
               'template spec' : True,   # save a plot of SSP spectra 
               'image format'  : 'png'}  # image type for plots

# Percentiles of each model parameter to report in the output file
param_percentiles = [5, 16, 50, 84, 95]

# When running in parallel mode, utilize (Total cores) - reserved_cores
reserved_cores = 2 # integer

# Input emission line strengths
#   keys are emission line name (str) corresponding to Name in the input file
#   values are two-element tuple: (rest-frame wavelength (Angstroms), weight)
#   measurements in input file have column names Name_FLUX, Name_ERR 
#   corresponding to line flux and error (null=-99) and lines will only
#   contribute to the model likelihood if they appear in the input file
emline_list_dict = {'OII' : (3727., 0.5), 'OIII' : (5007., 0.5),
                    'Hb'  : (4861., 1.),  'Ha' : (6563., 1.),
                    'NII' : (6583., 0.5), 'NeIII': (3869., 0.5)
                   }

emline_factor = 1e-17 # numerical conversion from input values to units ergs/cm2/s

# Input absorption line indices
#   keys are index name (str) corresponding to Name in input file
#   values are list of [weight, index band, blue continuum, red continuum, units]
#   units key: 0=Angstroms, 1=mag, 2=flux ratio (red/blue)
#   measurements in input file have column names Name_INDX, Name_Err
#   corresponding to index measurement and error (null=-99) and lines will only
#   contribute to the model likelihood if they appear in the input file
#
#                                        weight     index              blue continuum        red continuum     units
absorption_index_dict = {"Lick_CN1"    : [ 1., (4142.125, 4177.125), (4080.125, 4117.625), (4244.125, 4284.125), 1],
                         "Lick_CN2"    : [ 1., (4142.125, 4177.125), (4083.875, 4096.375), (4244.125, 4284.125), 1],
                         "Lick_Ca4227" : [ 1., (4222.25, 4234.75),   (4211.0, 4219.75),    (4241.0, 4251.0),     0],
                         "Lick_G4300"  : [ 1., (4281.375, 4316.375), (4266.375, 4282.625), (4318.875, 4335.125), 0],
                         "Lick_Fe4383" : [ 1., (4369.125, 4420.375), (4359.125, 4370.375), (4442.875, 4455.375), 0],
                         "Lick_Ca4455" : [ 1., (4452.125, 4474.625), (4445.875, 4454.625), (4477.125, 4492.125), 0],
                         "Lick_Fe4531" : [ 1., (4514.25, 4559.25),   (4504.25, 4514.25),   (4560.5, 4579.25),    0],
                         "Lick_Fe4668" : [ 1., (4634.0, 4720.25),    (4611.5, 4630.25),    (4742.75, 4756.5),    0],
                         "Lick_Hb"     : [ 1., (4847.875, 4876.625), (4827.875, 4847.875), (4876.625, 4891.625), 0],
                         "Lick_Fe5015" : [ 1., (4977.75, 5054.0),    (4946.5, 4977.75),    (5054.0, 5065.25),    0],
                         "Lick_Mg1"    : [ 1., (5069.125, 5134.125), (4895.125, 4957.625), (5301.125, 5366.125), 1],
                         "Lick_Mg2"    : [ 1., (5154.125, 5196.625), (4895.125, 4957.625), (5301.125, 5366.125), 1],
                         "Lick_Mgb"    : [ 1., (5160.125, 5192.625), (5142.625, 5161.375), (5191.375, 5206.375), 0],
                         "Lick_Fe5270" : [ 1., (5245.65, 5285.65),   (5233.15, 5248.15),   (5285.65, 5318.15),   0],
                         "Lick_Fe5335" : [ 1., (5312.125, 5352.125), (5304.625, 5315.875), (5353.375, 5363.375), 0],
                         "Lick_Fe5406" : [ 1., (5387.5, 5415.0),     (5376.25, 5387.5),    (5415.0, 5425.0),     0],
                         "Lick_Fe5709" : [ 1., (5696.625, 5720.375), (5672.875, 5696.625), (5722.875, 5736.625), 0],
                         "Lick_Fe5782" : [ 1., (5776.625, 5796.625), (5765.375, 5775.375), (5797.875, 5811.625), 0],
                         "Lick_NaD"    : [ 1., (5876.875, 5909.375), (5860.625, 5875.625), (5922.125, 5948.125), 0],
                         "Lick_TiO1"   : [ 1., (5936.625, 5994.125), (5816.625, 5849.125), (6038.625, 6103.625), 1],
                         "Lick_TiO2"   : [ 1., (6189.625, 6272.125), (6066.625, 6141.625), (6372.625, 6415.125), 1],
                         "Lick_Hd_A"   : [ 1., (4083.5, 4122.25),    (4041.6, 4079.75),    (4128.5, 4161.0),     0],
                         "Lick_Hg_A"   : [ 1., (4319.75, 4363.5),    (4283.5, 4319.75),    (4367.25, 4419.75),   0],
                         "Lick_Hd_F"   : [ 1., (4091.0, 4112.25),    (4057.25, 4088.5),    (4114.75, 4137.25),   0],
                         "Lick_Hg_F"   : [ 1., (4331.25, 4352.25),   (4283.5, 4319.75),    (4354.75, 4384.75),   0],
                         "D4000"       : [ 1., (False, False),       (3750., 3950.),       (4050., 4250.),       2]
                        }
 
