""" MCSED - dust.py


1) Dust Laws
    a) Calzetti:
    b) Noll:

.. moduleauthor:: Greg Zeimann <gregz@astro.as.utexas.edu>

"""

import numpy as np


def calzettilaw(wave):
    ''' Calzetti et al. (2000) dust attenuation curve, k(wave)

    Parameters
    ----------
    wave : numpy array (1 dim)
        wavelength

    Returns
    -------
    k : numpy array (1 dim)
        A(wave) = R_V * k(wave)
    '''
    invwv = 1/(wave/1e4)
    sel1 = np.nonzero(wave < 0.63e4)[0]
    sel2 = np.nonzero(np.logical_and(wave >= 0.63e4, wave < 2.2e4))[0]
    k1 = np.zeros(sel1.shape)
    k2 = np.zeros(sel2.shape)
    k1 = (2.659 * (-2.156 + 1.509 * invwv[sel1] - 0.198 * invwv[sel1]**2 +
          0.011 * invwv[sel1]**3) + 4.05)
    k2 = 2.659*(-1.857 + 1.040*invwv[sel2]) + 4.05
    k = np.zeros(wave.shape)
    k[sel1] = k1
    k[sel2] = k2
    return k


class calzetti:
    ''' Calzetti Dust Law
    '''
    def __init__(self, tau=0.7, tau_lims=[-0.2, 3.0], tau_delta=0.2):
        ''' Initialize Class

        Parameters
        -----
        tau : float
            Effective depth, e.g., Observed = True * exp**(-tau/4.05 * k(wave))
        '''
        self.tau = tau
        self.tau_lims = tau_lims
        self.tau_delta = tau_delta
        self.nparams = 1
        self.calz = None

    def get_params(self):
        ''' Return current parameters '''
        return [self.tau]

    def get_param_lims(self):
        ''' Return current parameter limits '''
        return [self.tau_lims]

    def get_param_deltas(self):
        ''' Return current parameter deltas '''
        return [self.tau_delta]

    def get_names(self):
        ''' Return names of each parameter '''
        return ['$tau_{dust}$']

    def prior(self):
        ''' Uniform prior based on boundaries '''
        tau_flag = (self.tau > self.tau_lims[0])*(self.tau < self.tau_lims[1])
        return tau_flag

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
        self.tau = input_list[start_value]

    def plot(self, ax, wave, color=[0/255., 175/255., 202/255.]):
        ''' Plot Dust Law for given set of parameters '''
        dust = self.evaluate(wave)
        ax.plot(wave, dust, color=color, alpha=0.4)

    def evaluate(self, wave):
        ''' Evaluate Dust Law

        Parameters
        ----------
        wave : numpy array (1 dim)
            wavelength

        Returns
        -------
        taulam : numpy array (1 dim)
            Effective optical depth as a function of wavelength
        '''
        if self.calz is None:
            self.calz = calzetti(wave)
        taulam = self.tau / 4.05 * self.calz
        return taulam


class noll:
    ''' Prescription for dust law comes from Noll et al. (2009), with constants
    defined in Kriek & Conroy (2013).  This dust attenuation law includes a
    bump at 2175A and a modified Calzetti et al. (2000) attenuation curve.

        A(wave) = frac{A_V}{4.05} (k'(wave) + D(wave))
                     left(frac{wave}{5500}right) ^delta
        D(wave) = frac{E_b (wave,dellam)^2 }{(wave^2-lam0^2)^2
                     + (wave,dellam)^2}
    '''
    def __init__(self, tau=0.7, delta=0.0, Eb=1.0, tau_lims=[-0.2, 3.0],
                 delta_lims=[-1., 1.], Eb_lims=[-0.2, 6.], tau_delta=0.2,
                 delta_delta=0.2, Eb_delta=0.4):
        ''' Initialize Class

        Parameters
        -----
        tau : float
            Effective depth, e.g., Observed = True * exp**(-tau/4.05 * k(wave))
        delta : float
            Power for powerlaw modification of Calzetti curve
        Eb : float
            Strength of 2175A bump.  See equation above for the Drude profile,
            D(wave)
        '''
        self.tau = tau
        self.delta = delta
        self.Eb = Eb
        self.tau_lims = tau_lims
        self.delta_lims = delta_lims
        self.Eb_lims = Eb_lims
        self.tau_delta = tau_delta
        self.delta_delta = delta_delta
        self.Eb_delta = Eb_delta
        self.nparams = 3
        self.calz = None

    def get_params(self):
        ''' Return current parameters '''
        return [self.tau, self.delta, self.Eb]

    def get_param_lims(self):
        ''' Return current parameter limits '''
        return [self.tau_lims, self.delta_lims, self.Eb_lims]

    def get_param_deltas(self):
        ''' Return current parameter deltas '''
        return [self.tau_delta, self.delta_delta, self.Eb_delta]

    def get_names(self):
        ''' Return names of each parameter '''
        return ['$tau_{dust}$', '$\delta$', '$E_b$']

    def prior(self):
        ''' Uniform prior based on boundaries '''
        tau_flag = (self.tau > self.tau_lims[0])*(self.tau < self.tau_lims[1])
        delta_flag = ((self.delta > self.delta_lims[0]) *
                      (self.delta < self.delta_lims[1]))
        Eb_flag = (self.Eb > self.Eb_lims[0])*(self.Eb < self.Eb_lims[1])
        return tau_flag * delta_flag * Eb_flag

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
        self.tau = input_list[start_value]
        self.delta = input_list[start_value+1]
        self.Eb = input_list[start_value+2]

    def plot(self, ax, wave, color=[0/255., 175/255., 202/255.]):
        ''' Plot Dust Law for given set of parameters '''
        dust = self.evaluate(wave)
        ax.plot(wave, dust, color=color, alpha=0.4)

    def evaluate(self, wave):
        ''' Evaluate Dust Law

        Parameters
        ----------
        wave : numpy array (1 dim)
            wavelength

        Returns
        -------
        taulam : numpy array (1 dim)
            Effective optical depth as a function of wavelength
        '''
        dellam = 350.
        lam0 = 2175.
        if self.calz is None:
            self.calz = calzettilaw(wave)
        Dlam = (self.Eb * (wave*dellam)**2 /
                          ((wave**2-lam0**2)**2+(wave*dellam)**2))
        taulam = (self.tau / 4.05 *
                  (self.calz+Dlam)*(wave/5500)**(self.delta))
        return taulam
