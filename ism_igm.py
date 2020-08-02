"""

ISM/IGM corrections

"""

import numpy as np
from scipy.interpolate import interp1d, interp2d
from astropy.io import fits
from astropy.table import Table
import os.path as op
import sys
sys.path.insert(0, 'ISM_IGM')

def get_MW_EBV(args):
    ''' 
    This function reads the input file to obtain the coordinates for 
    the source objects. It is called only when ISM_correct_coords is not None, 
    since the only thing that depends on coordinates is the Milky Way dust 
    correction. Note: The program will assume degrees; this can be manually 
    changed here if the user wishes.

    The columns for coordinates should be named "C1" and "C2" if the sources 
    are not in the Skelton Catalog. For sources in Skelton, coordinates are 
    optional (since they will be automatically read from the Skelton catalog).

    Note: The program will assume degrees; this can be manually 
    changed here if the user wishes.

    Parameters
    ----------
    args : class
        The args class is carried from function to function with information
        from command line input and config.py
    Returns
    -------
    E(B-V) for Milky Way  : 1D array
    '''
    F = Table.read(args.filename, format='ascii')
    Fcols = F.colnames
    nobj = len(F['Field'])
    try:
        C1 = F['C1']
        C2 = F['C2']
    except:
        # Skelton catalogs
        fields = ['aegis', 'cosmos', 'goodsn', 'goodss', 'uds']
        name_base = '_3dhst.v4.1.cat.FITS'
        if F['Field'][0].lower() in fields:
            for fd in F['Field']: 
                #Make sure the input isn't a mix of Skelton and non-Skelton
                assert fd.lower() in fields, "%s not in Skelton"%(fd) 
            C1, C2 = np.zeros(nobj), np.zeros(nobj)
            field_dict = {}
            for field in fields:
                field_dict[field] = fits.open(op.join('3dhst_catalogs',
                                                    field+name_base))[1]
            for i, datum in enumerate(F):
                # assumes first element is field name and second is ID
                loc = datum[0].lower()
                C1[i] = field_dict[loc].data['ra'][int(datum[1])-1]
                C2[i] = field_dict[loc].data['dec'][int(datum[1])-1]
            args.ISM_correct_coords = 'FK5' #Skelton coordinates are RA and Dec
        else:
            print("No coordinates given and no match to Skelton Catalog")
            return np.array([np.nan]*nobj)

    from dustmaps.sfd import SFDQuery
    from astropy.coordinates import SkyCoord
    coords = SkyCoord(C1, C2, unit='deg', frame=args.ISM_correct_coords.lower())
    sfd = SFDQuery()
    ebv = sfd(coords)
    return ebv


def get_tauIGMf():
    """
    Statistical IGM correction using 2-D interpolation on results 
    from Madau (1995) equations
    
    Parameters
    ----------
    None

    Returns
    -------
    tauIGMf: 2-D interpolation optical depth function
    (of observed wavelength and redshift from 0 to 4)
    """
    fulltauIGM = np.loadtxt("ISM_IGM/IGM_tau_z_0_4.dat")
    assert len(fulltauIGM[0,1:])>len(fulltauIGM[1:,0]) 
    tauIGMf = interp2d(fulltauIGM[0,1:],fulltauIGM[1:,0],fulltauIGM[1:,1:],
                       bounds_error=False,fill_value=0.0)
    return tauIGMf

def get_tauISMf():
    """
    Milky Way Dust correction using 1-D interpolation on wavelength given 
    2-D maps originally by Schlegel, Finkbeiner, and Davis (1998) and updated 
    by Schlafly and Finkbeiner (2011); 
    a Rv=3.1 Fitzpatrick (1999) extinction law is used

    Parameters
    ----------
    None

    Returns
    -------
    tauISMf: 1-D interpolation optical depth function 
    (of [observed] wavelength) using Rv=3.1 Fitzpatrick (1999) extinction law
    """
    wavISM, Rlam = np.loadtxt("ISM_IGM/Rlam_31.dat", unpack=True)
    tauISMf = interp1d(wavISM, Rlam,
                       bounds_error=False,fill_value=(max(Rlam),0.0))
    return tauISMf
