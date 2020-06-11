""" MCSED - dust.py


1) Dust Laws
    a) Calzetti:
    b) Noll:

.. moduleauthor:: Greg Zeimann <gregz@astro.as.utexas.edu>

"""

import numpy as np


def calzettilaw(wave, Rv=4.05):
    ''' Calzetti et al. (2000) dust attenuation curve, k(wave)

    Parameters
    ----------
    wave : numpy array (1 dim)
        wavelength in Angstroms
    Rv : float
        extinction factor
    Returns
    -------
    k : numpy array (1 dim)
        A(wave) / Av = k(wave) / Rv
    '''
# WPB DELETE
#    print('this is Rv:  '+str(Rv))
    invwv = 1/(wave/1e4)
    sel1 = np.nonzero(wave < 0.63e4)[0]
    sel2 = np.nonzero(np.logical_and(wave >= 0.63e4, wave < 2.2e4))[0]
    k1 = np.zeros(sel1.shape)
    k2 = np.zeros(sel2.shape)
    k1 = (2.659 * (-2.156 + 1.509 * invwv[sel1] - 0.198 * invwv[sel1]**2 +
          0.011 * invwv[sel1]**3) + Rv)
    k2 = 2.659*(-1.857 + 1.040*invwv[sel2]) + Rv
    k = np.zeros(wave.shape)
    k[sel1] = k1
    k[sel2] = k2
    return k


class calzetti:
    ''' Calzetti Dust Law
    '''
    def __init__(self, EBV=0.15, EBV_lims=[-0.05, 1.50], EBV_delta=0.02, 
                 Rv=4.05, EBV_old_young=None):
        ''' Initialize Class

        Parameters
        -----
        EBV : float
            Color excess
            The observed (B-V) color minus the intrinsic (B-V) color
            Relates to Av, the magnitude of extinction in the V band (5500A) via
                Av = Rv * EBV
        Rv : float (held fixed throughout the fitting)
            Extinction factor
        EBV_old_young : float (held fixed throughout the fitting, defined in config.py)
            coefficient between the attenuation applied to the
            diffuse and birth cloud components, such that
            E(B-V)_diffuse = EBV_old_young * E(B-V)_birth
        '''
        self.EBV = EBV
        self.EBV_lims = EBV_lims
        self.EBV_delta = EBV_delta
        self.calz = None
        self.Rv = Rv
        self.EBV_old_young = EBV_old_young

    def get_nparams(self):
        ''' Return number of parameters '''
        return 1

    def get_params(self):
        ''' Return current parameters '''
        return [self.EBV]

    def get_param_lims(self):
        ''' Return current parameter limits '''
        return [self.EBV_lims]

    def get_param_deltas(self):
        ''' Return current parameter deltas '''
        return [self.EBV_delta]

    def get_names(self):
        ''' Return names of each parameter '''
        return ['E(B-V)']

    def prior(self):
        ''' Uniform prior based on boundaries '''
        EBV_flag = (self.EBV > self.EBV_lims[0])*(self.EBV < self.EBV_lims[1])
        return EBV_flag

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
        self.EBV = input_list[start_value]

    def plot(self, ax, wave, color=[0/255., 175/255., 202/255.], alpha=0.2):
        ''' Plot Dust Law for given set of parameters '''
        dust = self.evaluate(wave)
        ax.plot(wave, dust, color=color, alpha=alpha)

    def evaluate(self, wave, new_wave=False):
        ''' Evaluate Dust Law

        Parameters
        ----------
        wave : numpy array (1 dim)
            wavelength
        new_wave : bool
            recompute k(wave), even if already present
            only used in emission line strengths: finer wavelength grid

        Returns
        -------
        Alam : numpy array (1 dim)
            Effective optical depth as a function of wavelength
            Observed = True * 10**(-0.4 * Av / Rv * k(wave))
            A(wave) = E(B-V) * k(wave) = Av / Rv * k(wave)
        '''
# WPB DELETE
#        print('this is self.Rv:  '+str(self.Rv))
        if self.calz is None:
            self.calz = calzettilaw(wave, self.Rv)
        if new_wave:
            kwave = calzettilaw(wave, self.Rv)
        else:
            kwave = self.calz
        Alam = self.EBV * kwave
        return Alam


class noll:
    ''' Prescription for dust law comes from Noll et al. (2009), with constants
    defined in Kriek & Conroy (2013).  This dust attenuation law includes a
    bump at 2175A and a modified Calzetti et al. (2000) attenuation curve.

        A(wave) = frac{A_V}{R_V} (k'(wave) + D(wave))
                     left(frac{wave}{5500}right) ^delta
        D(wave) = frac{E_b (wave,dellam)^2 }{(wave^2-lam0^2)^2
                     + (wave,dellam)^2}
    '''
    def __init__(self, EBV=0.15, delta=0.0, Eb=2.5, EBV_lims=[-0.05, 1.50],
                 delta_lims=[-1., 1.], Eb_lims=[-0.2, 6.], EBV_delta=0.02,
                 delta_delta=0.3, Eb_delta=1.0, 
                 Rv=4.05, EBV_old_young=None):
        ''' Initialize Class

        Parameters
        -----
        EBV : float
            Color excess
            The observed (B-V) color minus the intrinsic (B-V) color
            Relates to Av, the magnitude of extinction in the V band (5500A) via
                Av = Rv * EBV
        delta : float
            Power for powerlaw modification of Calzetti curve
        Eb : float
            Strength of 2175A bump.  See equation above for the Drude profile,
            D(wave)
        Rv : float (held fixed throughout the fitting)
            Extinction factor
        EBV_old_young : float (held fixed throughout the fitting, defined in config.py)
            coefficient between the attenuation applied to the
            diffuse and birth cloud components, such that
            E(B-V)_diffuse = EBV_old_young * E(B-V)_birth
        '''
        self.EBV = EBV
        self.delta = delta
        self.Eb = Eb
        self.EBV_lims = EBV_lims
        self.delta_lims = delta_lims
        self.Eb_lims = Eb_lims
        self.EBV_delta = EBV_delta
        self.delta_delta = delta_delta
        self.Eb_delta = Eb_delta
        self.calz = None
        self.Rv = Rv
        self.EBV_old_young = EBV_old_young

    def get_nparams(self):
        ''' Return number of parameters '''
        return 3

    def get_params(self):
        ''' Return current parameters '''
        return [self.EBV, self.delta, self.Eb]

    def get_param_lims(self):
        ''' Return current parameter limits '''
        return [self.EBV_lims, self.delta_lims, self.Eb_lims]

    def get_param_deltas(self):
        ''' Return current parameter deltas '''
        return [self.EBV_delta, self.delta_delta, self.Eb_delta]

    def get_names(self):
        ''' Return names of each parameter '''
        return ['E(B-V)', '$\delta$', '$E_b$']

    def prior(self):
        ''' Uniform prior based on boundaries '''
        EBV_flag = (self.EBV > self.EBV_lims[0])*(self.EBV < self.EBV_lims[1])
        delta_flag = ((self.delta > self.delta_lims[0]) *
                      (self.delta < self.delta_lims[1]))
        Eb_flag = (self.Eb > self.Eb_lims[0])*(self.Eb < self.Eb_lims[1])
        return EBV_flag * delta_flag * Eb_flag

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
        self.EBV = input_list[start_value]
        self.delta = input_list[start_value+1]
        self.Eb = input_list[start_value+2]

    def plot(self, ax, wave, color=[0/255., 175/255., 202/255.], alpha=0.2):
        ''' Plot Dust Law for given set of parameters '''
        dust = self.evaluate(wave)
        ax.plot(wave, dust, color=color, alpha=alpha)

    def evaluate(self, wave, new_wave=False):
        ''' Evaluate Dust Law

        Parameters
        ----------
        wave : numpy array (1 dim)
            wavelength
        new_wave : bool
            recompute k(wave), even if already present
            only used in emission line strengths: finer wavelength grid

        Returns
        -------
        Alam : numpy array (1 dim)
            Effective optical depth as a function of wavelength
        '''
# WPB DELETE
#        print('this is self.Rv:  '+str(self.Rv))
        dellam = 350.
        lam0 = 2175.
        if self.calz is None:
            self.calz = calzettilaw(wave, self.Rv)
        if new_wave:
            kwave = calzettilaw(wave, self.Rv)
        else:
            kwave = self.calz

        Dlam = (self.Eb * (wave*dellam)**2 /
                          ((wave**2-lam0**2)**2+(wave*dellam)**2))
        Alam = (self.EBV * (kwave+Dlam)*(wave/5500)**(self.delta))
        return Alam


class reddy:
    '''
    Reddy et al. (2015) derived an attenuation law for 1.4 < z < 2.6 galaxies 
    selected on the basis of emission-line detections on HST grism frames. 
    '''
    def __init__(self, EBV=0.15, EBV_lims=[-0.1, 2.40], EBV_delta=0.02, 
                 Rv=2.505, EBV_old_young=None):
        ''' Initialize Class

        Parameters
        -----
        EBV : float
            Color excess
            The observed (B-V) color minus the intrinsic (B-V) color
            Relates to Av, the magnitude of extinction in the V band (5500A) via
                Av = Rv * EBV
        Rv : float (held fixed throughout the fitting)
            Extinction factor
        EBV_old_young : float (held fixed throughout the fitting, defined in config.py)
            coefficient between the attenuation applied to the
            diffuse and birth cloud components, such that
            E(B-V)_diffuse = EBV_old_young * E(B-V)_birth
        '''
        self.EBV = EBV
        self.EBV_lims = EBV_lims
        self.EBV_delta = EBV_delta
        self.klam = None
        self.Rv = Rv
        self.EBV_old_young = EBV_old_young

    def get_nparams(self):
        ''' Return number of parameters '''
        return 1

    def get_params(self):
        ''' Return current parameters '''
        return [self.EBV]

    def get_param_lims(self):
        ''' Return current parameter limits '''
        return [self.EBV_lims]

    def get_param_deltas(self):
        ''' Return current parameter deltas '''
        return [self.EBV_delta]

    def get_names(self):
        ''' Return names of each parameter '''
        return ['E(B-V)']

    def prior(self):
        ''' Uniform prior based on boundaries '''
        EBV_flag = (self.EBV > self.EBV_lims[0])*(self.EBV < self.EBV_lims[1])
        return EBV_flag

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
        self.EBV = input_list[start_value]

    def plot(self, ax, wave, color=[0/255., 175/255., 202/255.], alpha=0.2):
        ''' Plot Dust Law for given set of parameters '''
        dust = self.evaluate(wave)
        ax.plot(wave, dust, color=color, alpha=alpha)

    def reddylaw(self, wave):
        ''' Reddy et al. (2015) dust attenuation curve, k(wave)
    
        Parameters
        ----------
        wave : numpy array (1 dim)
            wavelength in Angstroms
        Returns
        -------
        k : numpy array (1 dim)
            A(wave) / Av = k(wave) / Rv
        '''
# WPB DELETE
#        print('this is self.Rv:  '+str(self.Rv))
        Rv = self.Rv
        invwv = 1/(wave/1e4)
        sel1 = np.nonzero(wave < 0.60e4)[0]
        sel2 = np.nonzero(np.logical_and(wave >= 0.60e4, wave < 2.2e4))[0]
        k1 = np.zeros(sel1.shape)
        k2 = np.zeros(sel2.shape)
        k1 = (-5.726 + 4.004 * invwv[sel1] - 0.525 * invwv[sel1]**2 +
              0.029 * invwv[sel1]**3 + Rv)
        k2 = (-2.672 - 0.010 * invwv[sel2] + 1.532 * invwv[sel2]**2 -
              0.412 * invwv[sel2]**3 + Rv)
        k = np.zeros(wave.shape)
        k[sel1] = k1
        k[sel2] = k2
        return k

    def evaluate(self, wave, new_wave=False):
        ''' Evaluate Dust Law

        Parameters
        ----------
        wave : numpy array (1 dim)
            wavelength
        new_wave : bool
            recompute k(wave), even if already present
            only used in emission line strengths: finer wavelength grid

        Returns
        -------
        Alam : numpy array (1 dim)
            Effective optical depth as a function of wavelength
            Observed = True * 10**(-0.4 * Av / Rv * k(wave))
            A(wave) = E(B-V) * k(wave) = Av / Rv * k(wave)
        '''
        if self.klam is None:
            self.klam = self.reddylaw(wave)
        if new_wave:
            kwave = self.reddylaw(wave)
        else:
            kwave = self.klam
        Alam = self.EBV * kwave
        return Alam


class conroy:
    '''

    '''
    def __init__(self, EBV=0.15, B=0.5, EBV_lims=[-0.06, 2.00], 
                 B_lims=[-0.5, 1.5], EBV_delta=0.02, B_delta=0.2,
                 Rv=3.1, EBV_old_young=None):
        ''' Initialize Class

        Parameters
        -----
        EBV : float
            Color excess
            The observed (B-V) color minus the intrinsic (B-V) color
            Relates to Av, the magnitude of extinction in the V band (5500A) via
                Av = Rv * EBV
        B : float
            Strength of 2175A bump. 
            B=1 represents the standard Milky Way bump described 
            by Cardelli et al. (1989) for R_V = 3.1, 
            and B=0 corresponds to no excess attenuation at 2175A
        Rv : float (held fixed throughout the fitting)
            Extinction factor
        EBV_old_young : float (held fixed throughout the fitting, defined in config.py)
            coefficient between the attenuation applied to the
            diffuse and birth cloud components, such that
            E(B-V)_diffuse = EBV_old_young * E(B-V)_birth
        '''
        self.EBV = EBV
        self.B = B
        self.EBV_lims = EBV_lims
        self.B_lims = B_lims
        self.EBV_delta = EBV_delta
        self.B_delta = B_delta
        self.axlam = None
        self.bxlam = None
        self.Rv = Rv
        self.EBV_old_young = EBV_old_young

    def get_nparams(self):
        ''' Return number of parameters '''
        return 2

    def get_params(self):
        ''' Return current parameters '''
        return [self.EBV, self.B]

    def get_param_lims(self):
        ''' Return current parameter limits '''
        return [self.EBV_lims, self.B_lims]

    def get_param_deltas(self):
        ''' Return current parameter deltas '''
        return [self.EBV_delta, self.B_delta]

    def get_names(self):
        ''' Return names of each parameter '''
        return ['E(B-V)', 'B']

    def prior(self):
        ''' Uniform prior based on boundaries '''
        EBV_flag = (self.EBV > self.EBV_lims[0])*(self.EBV < self.EBV_lims[1])
        B_flag = (self.B > self.B_lims[0])*(self.B < self.B_lims[1])
        return EBV_flag * B_flag

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
        self.EBV = input_list[start_value]
        self.B = input_list[start_value+1]

    def plot(self, ax, wave, color=[0/255., 175/255., 202/255.], alpha=0.2):
        ''' Plot Dust Law for given set of parameters '''
        dust = self.evaluate(wave)
        ax.plot(wave, dust, color=color, alpha=alpha)

    def conroylaw(self, wave):
        ''' Conroy et al. (2010) dust attenuation curve, a(x) and b(x)
        where A(wave) / Av = a(x) + b(x) / Rv
        and x = 1 / wave (units inverse microns)
    
        Parameters
        ----------
        wave : numpy array (1 dim)
            wavelength in Angstroms
        Returns
        -------
        a : numpy array (1 dim)
        b : numpy array (1 dim)
        '''
        Rv = self.Rv
        x = 1/(wave/1e4)
# WPBWPB delete
#        sel1 = np.nonzero(np.logical_and(wave >= 0.125e4, wave < 0.17e4))[0]
#        sel2 = np.nonzero(np.logical_and(wave >= 0.17e4,  wave < 0.30e4))[0]
#        sel3 = np.nonzero(np.logical_and(wave >= 0.30e4,  wave < 0.91e4))[0]
#        sel4 = np.nonzero(np.logical_and(wave >= 0.91e4,  wave < 3.33e4))[0]

        sel1 = np.nonzero(np.logical_and(x >= 5.9, x < 8.0))[0]
        sel2 = np.nonzero(np.logical_and(x >= 3.3, x < 5.9))[0]
        sel3 = np.nonzero(np.logical_and(x >= 1.1, x < 3.3))[0]
        sel4 = np.nonzero(np.logical_and(x >= 0.3, x < 1.1))[0]

        a1, b1 = np.zeros(sel1.shape), np.zeros(sel1.shape)
        a2, b2 = np.zeros(sel2.shape), np.zeros(sel2.shape)
        a3, b3 = np.zeros(sel3.shape), np.zeros(sel3.shape)
        a4, b4 = np.zeros(sel4.shape), np.zeros(sel4.shape)

        x1 = x[sel1]
        fa1 = (-0.0447 * (x1-5.9)**2. 
               - 0.00978 * (x1-5.9)**3.)
        fb1 = (0.213 * (x1-5.9)**2. 
               + 0.121 * (x1-5.9)**3.)
        a1 = (fa1 + 1.752 - 0.316 * x1 
              - (0.104*self.B) / ((x1 - 4.67)**2. + 0.341) )
        b1 = (fb1 - 3.09 + 1.825 * x1
              + (1.206*self.B) / ((x1 - 4.62)**2. + 0.263) )

        x2 = x[sel2]
        fa2 = ( (3.3/x2)**6. * (-0.0370 + 0.0469 * self.B
                - 0.601 * self.B / Rv + 0.542 / Rv) )
        a2 = (fa2 + 1.752 - 0.316 * x2
              - (0.104*self.B) / ((x2 - 4.67)**2. + 0.341) )
        b2 = (-3.09 + 1.825 * x2
              + (1.206*self.B) / ((x2 - 4.62)**2. + 0.263) )

        x3 = x[sel3]
        y = x3 - 1.82
        a3 = (1. + 0.177*y - 0.504*y**2. - 0.0243*y**3. + 0.721*y**4.
              + 0.0198*y**5. - 0.775*y**6. + 0.330*y**7.)
        b3 = (1.413*y + 2.283*y**2. + 1.072*y**3. - 5.384*y**4.
              - 0.622*y**5. + 5.303*y**6. - 2.090*y**7.)

        x4 = x[sel4]
        a4 = 0.574 * x4 ** 1.61
        b4 = -0.527 * x4 ** 1.61


        a, b = np.zeros(wave.shape), np.zeros(wave.shape)
        a[sel1] = a1
        b[sel1] = b1
        a[sel2] = a2
        b[sel2] = b2
        a[sel3] = a3
        b[sel3] = b3
        a[sel4] = a4
        b[sel4] = b4

        return a, b

    def evaluate(self, wave, new_wave=False):
        ''' Evaluate Dust Law

        Parameters
        ----------
        wave : numpy array (1 dim)
            wavelength
        new_wave : bool
            recompute k(wave), even if already present
            only used in emission line strengths: finer wavelength grid

        Returns
        -------
        Alam : numpy array (1 dim)
            Effective optical depth as a function of wavelength
            A(wave) / Av = a(x) + b(x) / Rv
            and x = 1 / wave (units inverse microns)
            and Av = Rv * EBV 
        '''
        if (isinstance(wave, float)) | (isinstance(wave, int)):
            wave = np.array([wave])

        if True: #(self.axlam is None) | (new_wave):
            a, b = self.conroylaw(wave)
            self.axlam = a
            self.bxlam = b
        else:
            # update terms that depend on B (free parameter)
            # the others only depend on wavelength
            x = 1./(wave/1e4)

            sel1 = np.nonzero(np.logical_and(x >= 5.9, x < 8.0))[0]
            sel2 = np.nonzero(np.logical_and(x >= 3.3, x < 5.9))[0]

            a1, b1 = np.zeros(sel1.shape), np.zeros(sel1.shape)
            a2, b2 = np.zeros(sel2.shape), np.zeros(sel2.shape)

            x1 = x[sel1]
            fa1 = (-0.0447 * (x1-5.9)**2.
                   - 0.00978 * (x1-5.9)**3.)
            fb1 = (0.213 * (x1-5.9)**2.
                   + 0.121 * (x1-5.9)**3.)
            a1 = (fa1 + 1.752 - 0.316 * x1
                  - (0.104*self.B) / ((x1 - 4.67) + 0.341) )
            b1 = (fb1 - 3.09 + 1.825 * x1
                  + (1.206*self.B) / ((x1 - 4.62)**2. + 0.263) )

            x2 = x[sel2]
            fa2 = ( (3.3/x2)**6. * (-0.0370 + 0.0469 * self.B
                    - 0.601 * self.B / self.Rv + 0.542 / self.Rv) )
            a2 = (fa2 + 1.752 - 0.316 * x2
                  - (0.104*self.B) / ((x2 - 4.67) + 0.341) )
            b2 = (-3.09 + 1.825 * x2
                  + (1.206*self.B) / ((x2 - 4.62)**2. + 0.263) )

            self.axlam[sel1] = a1
            self.bxlam[sel1] = b1
            self.axlam[sel2] = a2
            self.bxlam[sel2] = b2

        axlam = self.axlam
        bxlam = self.bxlam
        Rv = self.Rv
        Av = Rv * self.EBV

        Alam = Av * ( axlam + bxlam / Rv )
        return Alam

