""" MCSED - config.py

1) Configuration Settings

.. moduleauthor:: Greg Zeimann <gregz@astro.as.utexas.edu>

"""
# SSP code for models
ssp = 'fsps'  # options include: 'fsps'
isochrone = 'padova'  # options include: 'padova'
# SFH options include: 'constant', 'burst', 'polynomial', 'exponential', 
#                      'double_powerlaw', 'empirical_direct', 'empirical',
sfh = 'constant' #'empirical_direct'
dust_law = 'calzetti' # options include: 'noll', 'calzetti'
dust_em = 'DL07'  # options include: 'DL07'

# Dust attenuation law parameters
# if set to a negative value, use the default value for dust law of choice
Rv = -1 # extinction factor
# The relative attenuation between the birth cloud and the diffuse component
# in the dust model 
# such that E(B-V)_diffuse = EBV_stars_gas * E(B-V)_birth
EBV_stars_gas = -1

t_birth = 7. # age of the birth cloud (log years)

# Fit dust emission parameters
# If True, fit the dust emission component. 
# If False, remove all filters redward of rest-frame wave_dust_em microns 
# and fix dust emission parameters to umin=2.0, gamma=0.05, qpah=2.5 
fit_dust_em  = False 
wave_dust_em = 2.5 # rest-frame wavelength in microns 
# assume energy balance or normalize the dust IR spectrum as a free parameter
assume_energy_balance = False 

# EMCEE parameters
nwalkers = 100 
nsteps   = 1000 

# Number of test objects
nobjects = 5
test_zrange = (1.9, 2.35) # redshift range of test objects (uniform prior)

# Nebular Emission Properties
# The ionization parameter, logU, is held fixed
logU = -2.5

# minimum fractional errors in observed photometry, 
# emission line fluxes, and absorption line indices 
phot_floor_error    = 0.10
emline_floor_error  = 0.10
absindx_floor_error = 0.10

# fractional error expected from the models, i.e., fractional error adopted
# for model photometry, emission line fluxes, and absorption line indices
model_floor_error = 0.10

# Use input data (photometry, emission lines, absorption line indices)
# If True, use additional data provided in the input file
# else, ignore input data (in which case input IDs must match Skelton+14 IDs)
use_input_data = True #False

# Input emission line strengths
# keys are emission line name (str) corresponding to name in input file
# values are two-element tuple: (rest-frame wavelength (Angstroms), weight)
# WPBWPB describe the weight
# see documentation XXXX for additional information
# WPB edit (e.g., OIII is not the blended feature, 
#           Balmer lines corrected for absorption, 
#           keys must match input columns of form _FLUX, _ERR...)
#           will only be used if present in input file,
#           must have null value = -99
#           must have both flux and error - cannot have flux with null error
#           can also set to {} or None, if preferred
emline_list_dict = {'OII' : (3727., 0.5), 'OIII' : (5007., 0.5),
                    'Hb' : (4861., 1.),   'Ha' : (6563., 1.),
                    'NII' : (6583., 0.5)
                   }

emline_factor = 1e-17 # numerical conversion from input values to units ergs/cm2/s

# Use metallicity-mass relationship from Ma et al. 2016
# NOTE: currently unavailable
metallicity_mass_relationship = False

# If False, leave metallicity as a free parameter
# else, must be float: fixed metallicity of SSP models
# if True, metallicity is fixed at 0.0077 (where Zsolar = 0.019)
metallicity = 0.0077 #False #0.0077 

# Output files
#   Supported image formats: eps, pdf, pgf, png, ps, raw, rgba, svg, svgz
output_dict = {'parameters'    : True,   # fitted parameters
               'settings'      : True,   # user-defined fitting assumptions
               'fitposterior'  : False,  # parameter posterior distributions
               'bestfitspec'   : True,   # best-fit SED model
               'fluxdensity'   : True,   # modeled and observed photometry
               'lineflux'      : True,  # modeled and observed emission lines
               'absorption'    : True,   # modeled, observed absorption indices
               'triangle plot' : True,   # summary diagnostic plot
               'sample plot'   : False,  # 
               'image format'  : 'png'}  # image type for plots

# percentiles of each parameter to report in the output file
param_percentiles = [5, 16, 50, 84, 95]

# When running in parallel mode, utilize (Total cores) - reserved_cores
reserved_cores = 2 # integer

# Input absorption line indices
# keys are index name (str) corresponding to name in input file
# values are list of [weight, index band, blue continuum, red continuum, units]
# weight: contribution to likelihood function (1 = equal to photometric point)
# units: 0=Angstroms, 1=mag, 2=flux ratio (red/blue)
# indices will only be used if they appear in the input file
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
 
