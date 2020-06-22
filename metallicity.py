""" MCSED - metallicity.py
  
module for handling the free or fixed stellar metallicity parameter

.. moduleauthor:: Will Bowman <wpb.astro@gmail.com>

"""

class stellar_metallicity:
    '''
    Define the stellar metallicity class

    stellar metallicity can be treated as a free model parameter 
    (set metallicity=False in config.py)

    or is fixed at a single value Z (where Zsolar=0.019)
    (set metallicity=Z in config.py)
    '''
    def __init__(self, met=-0.39, met_lims=[-1.98, 0.2], met_delta=0.3,
                 fix_met=False):
        ''' Initialize this class

        Parameters
        ----------
        met : float
            the stellar metallicity (in log solar units)
        fix_met : bool
            treat stellar metallicity as a fixed or free parameter
        '''
        self.met = met
        self.met_lims = met_lims
        self.met_delta = met_delta
        self.fix_met = fix_met

    def get_nparams(self):
        ''' Return number of parameters '''
        if self.fix_met:
            return 0
        else:
            return 1

    def get_params(self):
        ''' Return current parameters '''
        l = []
        if not self.fix_met:
            l.append(self.met)
        return l

    def get_param_lims(self):
        ''' Return current parameters limits '''
        l = []
        if not self.fix_met:
            l.append(self.met_lims)
        return l

    def get_param_deltas(self):
        ''' Return current parameter deltas '''
        l = []
        if not self.fix_met:
            l.append(self.met_delta)
        return l

    def get_names(self):
        ''' Return names of each parameter '''
        l = []
        if not self.fix_met:
            l.append('Log Z')
        return l

    def prior(self):
        ''' Uniform prior based on boundaries '''
        flag = (self.met >= self.met_lims[0]) * (self.met <= self.met_lims[1])
        return flag

    def set_parameters_from_list(self, input_list, start_value):
        ''' Set parameters from a list and a start_value

        Parameters
        ----------
        input_list : list
            list of input parameters (could be much larger than number of
            parameters to be set)
        start_value : int
            initial index from list to read out parameters
        '''
        if not self.fix_met:
            self.met = input_list[start_value]



