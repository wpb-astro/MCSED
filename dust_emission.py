# -*- coding: utf-8 -*-
"""
Created on Tue Apr  3 11:09:49 2018

@author: gregz
"""

import numpy as np
from scipy.interpolate import LinearNDInterpolator
import matplotlib.pyplot as plt


class DL07:
    ''' Prescription for dust emission comes from Draine & Li (2007) '''
    def __init__(self, umin=2.0, gamma=0.05, qpah=2.5, mdust=7.0, umin_lims=[0.1, 25.0],
                 gamma_lims=[0, 1.], qpah_lims=[0.47, 4.58], mdust_lims=[4.5,10.0],umin_delta=0.4,
                 gamma_delta=0.02, qpah_delta=1.0, mdust_delta=0.3, fixed=True, assume_energy_balance=False):
        ''' Initialize Class

        Parameters
        -----
        umin : float
            Minimum in distribution of incident intensity
        gamma : float
            Fraction of distribution in incident intensity greater than umin
        qpah : float
            Percentage of PAH emission
        mdust: float
            Log10(M_dust/M_sun), where M_dust is the total mass of dust in the galaxy
        assume_energy_balance: 
            If true, use energy balance between dust attenuation and emission to normalize dust IR spectrum
            If false, dust IR spectrum normalization is a free parameter
        '''
        self.umin = umin
        self.gamma = gamma
        self.qpah = qpah
        self.mdust = mdust
        self.umin_lims = umin_lims
        self.gamma_lims = gamma_lims
        self.qpah_lims = qpah_lims
        self.mdust_lims = mdust_lims
        self.umin_delta = umin_delta
        self.gamma_delta = gamma_delta
        self.qpah_delta = qpah_delta
        self.mdust_delta = mdust_delta
        self.fixed = fixed
        self.assume_energy_balance = assume_energy_balance
        self.get_dust_emission_matrix()

    def get_dust_emission_matrix(self):
        '''Draine & Li (2007) provide ascii files giving dust emissivity at various values of umin and qpah.
        We specifically use the MW3.1 models with Umax = 1e6.
        This function creates a 2D linear interpolator for dust emissivity (units of uJy*(10pc)^2/M_sun)
        The dust emissivity has two components: emission of dust heated by starlight at Umin (DE1) and emission of dust
        heated by starlight with U>Umin (DE2)
        '''
        DE1, DE2 = (np.zeros((7, 22, 1001)), np.zeros((7, 22, 1001)))
        #DE1mod, DE2mod = (np.zeros((7, 22, 1001)), np.zeros((7, 22, 1001)))
        self.qpaharray = np.array([0.47, 1.12, 1.77, 2.50, 3.19, 3.90, 4.58]) #Set of qpah values in Draine & Li tables
        self.htodarray = np.array([100., 100., 99.010, 98.039, 98.039, 97.087, 96.154]) #Ratio of total hydrogen to dust mass in the model for each qpah value
        self.uminarray = np.array([0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.7, 0.8,
                                   1.0, 1.2, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 7.0,
                                   8.0, 12.0, 15.0, 20.0, 25.0]) #Set of umin values in Draine & Li tables
        Xgrid, Ygrid = np.meshgrid(self.qpaharray, self.uminarray)
        X = np.vstack([Xgrid.ravel(), Ygrid.ravel()]).swapaxes(0, 1)

        for j, i in enumerate(np.arange(0, 70, 10)):
            de = np.loadtxt('DUSTEMISSION/DL07/DL07_MW3.1_%02d.dat' % i)
            self.wave = de[:, 0] * 1e4
            dnue = np.abs(np.hstack([0, np.diff(3e18 / self.wave)]))
            for k, v in enumerate(np.arange(1, de.shape[1], 2)):
                #norm = np.dot(dnue, de[:, v])
                DE1[j, k, :] = de[:, v]
                #DE1mod[j,k,:] = de[:,v]
                #norm = np.dot(dnue, de[:, v+1])
                DE2[j, k, :] = de[:, v+1]
                #DE2mod[j,k,:] = de[:,v+1]
                #print "DE1[%d,%d,0] = %.3e"%(j,k,DE1[j,k,0])
                #print "DE2[%d,%d,0] = %.3e"%(j,k,DE2[j,k,0])
        shape = DE1.shape
        self.interpumin = LinearNDInterpolator(X, DE1.reshape(shape[0] *
                                                              shape[1],
                                                              shape[2],order='F'))
        self.interpumax = LinearNDInterpolator(X, DE2.reshape(shape[0] *
                                                              shape[1],
                                                              shape[2],order='F'))
    
    def get_nparams(self):
        ''' Return number of parameters '''
        if self.fixed:
            return 0
        else:
            if self.assume_energy_balance:
                return 3
            else:
                return 4 

    def get_params(self):
        ''' Return current parameters '''
        if self.fixed:
            return []
        else:
            if self.assume_energy_balance:
                return [self.umin, self.gamma, self.qpah]
            else:
                return [self.umin, self.gamma, self.qpah, self.mdust]

    def get_param_lims(self):
        ''' Return current parameter limits '''
        if self.fixed:
            return []
        else:
            if self.assume_energy_balance:
                return [self.umin_lims, self.gamma_lims, self.qpah_lims]
            else:
                return [self.umin_lims, self.gamma_lims, self.qpah_lims, self.mdust_lims]

    def get_param_deltas(self):
        ''' Return current parameter deltas '''
        if self.fixed:
            return []
        else:
            if self.assume_energy_balance:
                return [self.umin_delta, self.gamma_delta, self.qpah_delta]
            else:
                return [self.umin_delta, self.gamma_delta, self.qpah_delta, self.mdust_delta]

    def get_names(self):
        ''' Return names of each parameter '''
        if self.fixed:
            return []
        else:
            if self.assume_energy_balance:
                return ['Umin', '$\gamma$', '$q_{pah}$']
            else:
                return ['Umin', '$\gamma$', '$q_{pah}$', '$M_{dust}$']

    def prior(self):
        ''' Uniform prior based on boundaries '''
        umin_flag = ((self.umin > self.umin_lims[0]) *
                     (self.umin < self.umin_lims[1]))
        gamma_flag = ((self.gamma > self.gamma_lims[0]) *
                      (self.gamma < self.gamma_lims[1]))
        qpah_flag = ((self.qpah > self.qpah_lims[0]) *
                     (self.qpah < self.qpah_lims[1]))
        mdust_flag = ((self.mdust > self.mdust_lims[0]) *
                     (self.mdust < self.mdust_lims[1]))
        if self.assume_energy_balance:
            return umin_flag * gamma_flag * qpah_flag
        else:
            return umin_flag * gamma_flag * qpah_flag * mdust_flag

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
        if not self.fixed:
            self.umin = input_list[start_value]
            self.gamma = input_list[start_value+1]
            self.qpah = input_list[start_value+2]
            if not self.assume_energy_balance:
                self.mdust = input_list[start_value+3]

    def plot(self, ax, wave, color=[0/255., 175/255., 202/255.], alpha=0.2):
        ''' Plot Dust Law for given set of parameters '''
        dustem = self.evaluate(wave)
        ax.plot(wave, dustem, color=color, alpha=alpha)

    def plotpractice(self, wave, z):
        #Distance in Mpc
        from cosmology import Cosmology
        C = Cosmology()
        D = C.luminosity_distance(z)
        dustem1 = self.evaluate(wave)
        dustem = np.interp(wave, wave * (1. + z),
                        dustem1 * (1. + z))
        dustem/=D**2 #Divide by distance ^ 2 (to get flux)
        dustem*=1.0e-29 #Go from uJy to erg/cm^2/s/Hz
        nu = 3.0e18/wave
        plt.figure()
        plt.loglog(wave/1.0e4, nu*dustem, 'b-')
        xlims = np.array([2.0,1000.])/(1.+z)
        plt.xlim(xlims)
        #plt.ylim(nu[-1]*dustem[-1],nu[-1]*dustem[-1]/3.0 * 2.0e5)
        #yl = 5.0e20*1.0e-29/(1.0e5*D)**2
        #yu = 3.0e25*1.0e-29/(1.0e5*D)**2
        #plt.ylim(yl,yu)
        plt.xlabel(r"$\lambda$ ($\mu m$)")
        #plt.ylabel(r"Dust emissivity ($\mu$Jy (10pc)$^2$ Hz M$_{sun}^{-1}$) ")
        plt.ylabel(r"Dust emissivity (erg cm$^{-1}$ s$^{-1}$) ")
        plt.savefig("DustTesting/%.2f_%.4f_%.2f_%.3f.png"%(self.umin,self.gamma,self.qpah,self.mdust))
        plt.show()
        #plt.close()

    def evaluate(self, wave):
        ''' Evaluate Dust Law

        Parameters
        ----------
        wave : numpy array (1 dim)
            wavelength in Angstroms

        Returns
        -------
        DustE : numpy array (1 dim)
            Dust emission spectrum (flux in uJy 10 pc)--linear combination of emission from dust heated by starlight at Umin and dust heated by starlight at U>Umin
        '''
        DustE = (self.interpumin(self.qpah, self.umin) * (1. - self.gamma) +
                 self.interpumax(self.qpah, self.umin) * self.gamma) #Units Jy*cm^2/sr/H (hydrogen atom)--luminosity/solid angle/mass

        DustE *= 1.249e24 * np.interp(self.qpah,self.qpaharray,self.htodarray) #Converting from Jy*cm^2/sr/H to units uJy/M_sun (flux/mass unit) at 10 pc assuming isotropic emission and then multiplying by Hydrogen to dust mass ratio so that the factor to get quantity so that it just needs to be multiplied by dust mass to get the proper normalization
        if not self.assume_energy_balance:
            DustE *= 10**self.mdust #Multiplying by total dust mass to get to flux in uJy at 10 pc

        return np.interp(wave, self.wave, DustE)

    # def evaluate2(self, wave):
    #     ''' Evaluate Dust Law

    #     Parameters
    #     ----------
    #     wave : numpy array (1 dim)
    #         wavelength in Angstroms

    #     Returns
    #     -------
    #     DustE : numpy array (1 dim)
    #         Dust emission spectrum (units Jy cm^2 sr^-1 H^-1)--linear combination of emission from dust heated by starlight at Umin and dust heated by starlight at U>Umin
    #     '''
    #     DustE = (self.interpumin2(self.qpah, self.umin) * (1. - self.gamma) +
    #              self.interpumax2(self.qpah, self.umin) * self.gamma)

    #     return np.interp(wave, self.wave, DustE)
