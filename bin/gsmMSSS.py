#!/usr/bin/env python
import sys
import warnings
import argparse
from copy import copy
import numpy as np
from astropy.table import Table
try:
    import pyvo
except:
    print("Cannot import the pyvo module. You will need to provide a catalogue file to the gsmMSSS function to work.")

## The following variables are for convenience
ARCSECTORAD = float('4.8481368110953599358991410235794797595635330237270e-6')
RADTOARCSEC = float('206264.80624709635515647335733077861319665970087963')
SECTORAD    = float('7.2722052166430399038487115353692196393452995355905e-5')
RADTOSEC    = float('13750.987083139757010431557155385240879777313391975')
RADTODEG    = float('57.295779513082320876798154814105170332405472466564')
DEGTORAD    = float('1.7453292519943295769236907684886127134428718885417e-2')
DEGTOHRS    = float('6.6666666666666666666666666666666666666666666666e-2')
RADTOHRS    = float('3.8197186342054880584532103209403446888270314977710')
HRSTORAD    = float('2.6179938779914943653855361527329190701643078328126e-1')
PI          = float('3.1415926535897932384626433832795028841971693993751')
TWOPI       = float('6.2831853071795864769252867665590057683943387987502')
PIBYTWO     = float('1.5707963267948966192313216916397514420985846996876')
SECPERDAY   = float('86400.0')
SECPERJULYR = float('31557600.0')



##### #####
#####  Log
##### #####
'''
 2014-09-16,    v0.1,    R. Breton,    Created
'''



##--------------------------------------------------------- Some useful functions
def Cartesian_to_polar(x, y, radians=False):
    """
    Transforms cartesian coordinates to polar.

    Parameters
    ----------
    x, y : float, ndarray
        Cartesian coordinates.
    radians : bool
        If true, angles are in radians, otherwise degrees
 
    Examples   
    ----------
        r, theta = Cartesian_to_polar(x, y)
    """
    r = np.sqrt(x**2+y**2)
    theta = np.arctan2(y,x)
    if not radians:
        theta *= RADTODEG
    return r, theta

def Deg_to_dms(val):
    """
    Converts degrees to degrees, minutes, and seconds of arc.
    """
    if (val < 0.0): sign = -1
    else: sign = 1
    arc = np.fmod(np.fabs(val), 180)
    d = int(arc)
    arc = (arc - d) * 60.0
    m = int(arc)
    s = (arc - m) * 60.0
    if sign==-1 and d==0:
        return (sign * d, sign * m, sign * s)
    else:
        return (sign * d, m, s)

def Deg_to_hms(val):
    """
    Converts degrees to hours, minutes, and seconds of arc.
    """
    val = np.fmod(val, 360)
    if (val < 0.0): val = val + 360
    arc = DEGTOHRS * val
    h = int(arc)
    arc = (arc - h) * 60.0
    m = int(arc)
    s = (arc - m) * 60.0
    return (h, m, s)

def DEC_to_str(val, precision=2, sign=True):
    """
    Formats a numerical DEC float (degrees) to a nice
    string representation.

    Parameters
    ----------
    val : float
        DEC value to convert to a string representation.
    precision : int
        Number of decimal points to display for the seconds of arc.
    sign : bool
        If true, will force the display of +/- at the front of the string.
    """
    if sign:
        sign = '+'
    else:
        sign = '-'
    d, m, s = Deg_to_dms(val)
    return "{0:{sign}}.{1}.{2:.{prec}f}".format(d, m, s, sign=sign, prec=precision)

def Distance_on_sky(ra1, dec1, ra2, dec2, radians=False, one_to_one=False):
    """
    Returns the angular distance between two points in the sky.
    
    Parameters
    ----------
    ra1, dec1 : float, ndarray
        Coordinates of first point(s) (degree or radians)
    ra2, dec2 : float, ndarray
        Coordinates of second point(s) (degree or radians)
    radians : bool
        If true angles are in radians, otherwise degrees
    one_to_one : bool
        If true, will calculate the distance for each possible pair of
        (ra1,dec1) and (ra2,dec2). Hence the output will be shaped as
        ra1.size,ra2.size.

    Note
    ----------
    The dimensions of ra1/dec1 and ra2/dec2 must be compatible in the
    broadcasting sense. So either float vs ndarray, ndarray vs float or
    ndarray vs ndarray of same shapes for element-wise calculation.
    The exception would be in the case one_to_one = True.
    """
    if type(ra1) == type(str()):
        ra1 = psr_utils.ra_to_rad(ra1)
        dec1 = psr_utils.dec_to_rad(dec1)
        ra2 = psr_utils.ra_to_rad(ra2)
        dec2 = psr_utils.dec_to_rad(dec2)
        distance = Distance_on_sphere(ra1, PIBYTWO-dec1, ra2, PIBYTWO-dec2, radians=True, one_to_one=one_to_one)
        if not radians:
            distance *= RADTODEG
    elif not radians:
        distance = Distance_on_sphere(ra1, 90.-dec1, ra2, 90.-dec2, radians=radians, one_to_one=one_to_one)
    else:
        distance = Distance_on_sphere(ra1, PIBYTWO-dec1, ra2, PIBYTWO-dec2, radians=radians, one_to_one=one_to_one)
    return distance

def Distance_on_sphere(phi1, theta1, phi2, theta2, radians=False, one_to_one=False):
    """
    Returns the angular distance between two points on a sphere.
    
    Parameters
    ----------
    phi1, theta1 : float, ndarray
        Coordinates of first point(s) (degree or radians)
    phi2, theta2 : float, ndarray
        Coordinates of second point(s) (degree or radians)
    radians : bool
        If true angles are in radians, otherwise degrees
    one_to_one : bool
        If true, will calculate the distance for each possible pair of
        (phi1,theta1) and (phi2,theta2). Hence the output will be shaped as
        phi1.size,phi2.size.

    Note
    ----------
    The dimensions of phi1/theta1 and phi2/theta2 must be compatible in the
    broadcasting sense. So either float vs ndarray, ndarray vs float or
    ndarray vs ndarray of same shapes for element-wise calculation.
    The exception would be in the case one_to_one = True.
    """
    if one_to_one:
        phi1 = phi1[:,np.newaxis]
        theta1 = theta1[:,np.newaxis]
        phi2 = phi2[np.newaxis]
        theta2 = theta2[np.newaxis]
    x1, y1, z1 = Spherical_to_cartesian(1., theta1, phi1, radians=radians)
    x2, y2, z2 = Spherical_to_cartesian(1., theta2, phi2, radians=radians)
    cs = x1*x2 + y1*y2 + z1*z2
    xc = y1*z2 - z1*y2
    yc = z1*x2 - x1*z2
    zc = x1*y2 - y1*x2
    sn = np.sqrt(xc*xc + yc*yc + zc*zc)
    r, distance = Cartesian_to_polar(cs, sn, radians=radians)
    return distance

def gsmWriter(msssTable, outfile):
    """
    Takes an MSSS table and print it as gsm sky model.

    Parameters
    ----------
    msssTable : TableMSSS
        MSSS table object to dump the content for into a file
    outfile : str
        Output filename to store the sky model
    """
    with open(outfile, 'w') as f:
        ## Writing the preambule
        f.write("FORMAT = Name, Type, Ra, Dec, I, Q, U, V, "
            "ReferenceFrequency='150e6', SpectralIndex='[0.0]', MajorAxis, "
            "MinorAxis, Orientation\n\n# the next lines define the sources\n")

        ## Iterating over the rows
        for i,t in enumerate(msssTable):
            name = t['ID']
            ra = RA_to_str(t['RA'], precision=8)
            dec = DEC_to_str(t['DEC'], precision=8)
            I = t['A0_HBA']
            Q = ''
            U = ''
            V = ''
            freq = ''
            index = t['A1']
            majax = msssTable.MAJAX(i, hba=True, lba=False, sfflag=[0,1])
            minax = msssTable.MINAX(i, hba=True, lba=False, sfflag=[0,1])
            pa = msssTable.PA(i, hba=True, lba=False, sfflag=[0,1])
            if np.any(majax > 0.0) and np.any(minax > 0.0):
                sourcetype = 'GAUSSIAN'
                majax = majax.mean()
                minax = minax.mean()
                pa = pa.mean()
            else:
                sourcetype = 'POINT'

            entry = ("{name}, {sourcetype}, {ra}, {dec}, {I:.4f}, {Q}, {U}, {V}, "
                "{freq}, [{index:.4f}]").format(
                    name = name,
                    sourcetype = sourcetype,
                    ra = ra,
                    dec = dec,
                    I = I,
                    Q = Q,
                    U = U,
                    V = V,
                    freq = freq,
                    index = index
                    )
            if sourcetype == 'GAUSSIAN':
                entry += (", {majax:.1f}, {minax:.1f}, {pa:.1f}\n").format(
                    majax = majax,
                    minax = minax,
                    pa = pa
                    )
            else:
                entry += "\n"

            f.write(entry)
    return

def RA_to_str(val, precision=2, sign=False):
    """
    Formats a numerical RA float (degrees) to a nice
    string representation.

    Parameters
    ----------
    val : float
        RA value to convert to a string representation.
    precision : int
        Number of decimal points to display for the seconds of arc.
    sign : bool
        If true, will force the display of +/- at the front of the string.
    """
    if sign:
        sign = '+'
    else:
        sign = '-'
    h, m, s = Deg_to_hms(val)
    return "{0:{sign}}:{1}:{2:.{prec}f}".format(h, m, s, sign=sign, prec=precision)

def Spherical_to_cartesian(r, theta, phi, radians=False):
    """
    Transforms spherical coordinates to cartesian.
    
    Parameters
    ----------
    r, theta, phi : float, ndarray
        radius, inclination (polar), azimuth
    radians : bool
        if true, angles are in radians, otherwise degrees

    Examples
    ----------
        x, y, z = Spherical_to_cartesian(r, theta, phi)
    """
    if not radians:
        theta = theta*DEGTORAD
        phi = phi*DEGTORAD
    x = r*np.sin(theta)*np.cos(phi)
    y = r*np.sin(theta)*np.sin(phi)
    z = r*np.cos(theta)
    return x, y, z

class TableMSSS(Table):
    """
    This is a subclass of the astropy Table class tailored to the MSSS data.
    It contains convenience functions (Sint, e_Sint, Freqs, etc) that allow
    to retrieve the table value for a certain index. Since these quantities
    comprise mutliple columns, it is easier this way. Each of these functions
    takes the row index, and an optional sfflag parameters which is the type
    of data to return (None: the whole masked array entries, int: a single
    type of entries, list: a list of entries).

    For example, to get all Sint fluxes for SFFLAG 0 and 1 (i.e. PyBDSM) for
    row 15:
    >>> flux_pybdsm = Catalogue.Sint(15, [0,1])
    """
    warnings.simplefilter("ignore")
    
    def __init__(self, *args, **kwargs):
        ## Instantiate the object from the parent class
        super(TableMSSS, self).__init__(*args, **kwargs)
        self._init_cols()

    def _init_cols(self):
        ## Preparing the list of columns that pertain to each category of entries
        cols = copy(self.colnames)
        cols.sort()
        self._col_e_MAJAX = [s for s in cols if s.startswith('e_MAJAX')]
        self._col_e_MINAX = [s for s in cols if s.startswith('e_MINAX')]
        self._col_e_PA = [s for s in cols if s.startswith('e_PA')]
        self._col_e_Sint = [s for s in cols if s.startswith('e_Sint')]
        self._col_e_Spk = [s for s in cols if s.startswith('e_Spk')]
        self._col_MAJAX = [s for s in cols if s.startswith('MAJAX')]
        self._col_MINAX = [s for s in cols if s.startswith('MINAX')]
        self._col_PA = [s for s in cols if s.startswith('PA')]
        self._col_SFFLAG = [s for s in cols if s.startswith('SFFLAG')]
        self._col_Sint = [s for s in cols if s.startswith('Sint')]
        self._col_Spk = [s for s in cols if s.startswith('Spk')]
        self._freqs = np.array([float(s[-3:]) for s in self._col_SFFLAG])
        self._lba = np.flatnonzero(self._freqs < 100.)
        self._hba = np.flatnonzero(self._freqs >= 100.)
        if self._lba.size >= 1:
            self.has_lba = True
        else:
            self.has_lba = False
        self._hba = np.flatnonzero(self._freqs >= 100.)
        if self._hba.size >= 1:
            self.has_hba = True
        else:
            self.has_hba = False

    def __getitem__(self, item):
        out = Table.__getitem__(self, item)
        #out = super(TableMSSS, self).__getitem__(item)
        if isinstance(out, self.__class__):
            out._init_cols()
        return out

    def _getitem_special(self, index, items, sfflag=None, hba=True, lba=True):
        if hba and not lba:
            items = [items[i] for i in self._hba]
        elif not hba and lba:
            items = [items[i] for i in self._lba]
        t = np.ma.array([self.__getitem__(item)[index] for item in items])
        if isinstance(sfflag, int):
            flag = self.SFFLAG(index, hba=hba, lba=lba).data
            inds = flag == sfflag
            t = t[inds].data
            return t
        elif isinstance(sfflag, (tuple,list,np.ndarray)):
            flag = self.SFFLAG(index, hba=hba, lba=lba).data
            inds = sum([flag == isfflag for isfflag in sfflag]).astype(bool)
            t = t[inds].data
            return t
        else:
            return t

    @classmethod
    def read(cls, *args, **kwargs):
        return cls(Table.read(*args, **kwargs))

    def write(self, *args, **kwargs):
        tbl = Table(data=self.columns)
        tbl.write(*args, **kwargs)

    def e_MAJAX(self, index, sfflag=None, hba=True, lba=True):
        return self._getitem_special(index, self._col_e_MAJAX, sfflag=sfflag, hba=hba, lba=lba)

    def e_MINAX(self, index, sfflag=None, hba=True, lba=True):
        return self._getitem_special(index, self._col_e_MINAX, sfflag=sfflag, hba=hba, lba=lba)

    def e_PA(self, index, sfflag=None, hba=True, lba=True):
        return self._getitem_special(index, self._col_e_PA, sfflag=sfflag, hba=hba, lba=lba)

    def e_Sint(self, index, sfflag=None, hba=True, lba=True):
        return self._getitem_special(index, self._col_e_Sint, sfflag=sfflag, hba=hba, lba=lba)

    def e_Spk(self, index, sfflag=None, hba=True, lba=True):
        return self._getitem_special(index, self._col_e_Spk, sfflag=sfflag, hba=hba, lba=lba)

    def Freqs(self, index, sfflag=None, hba=True, lba=True):
        if isinstance(sfflag, int):
            flag = self.SFFLAG(index, hba=hba, lba=lba).data
            inds = flag == sfflag
            return self._freqs[inds]
        elif isinstance(sfflag, (tuple,list,np.ndarray)):
            flag = self.SFFLAG(index, hba=hba, lba=lba).data
            inds = sum([flag == isfflag for isfflag in sfflag]).astype(bool)
            return self._freqs[inds]
        else:
            if hba and not lba:
                freqs = [self._freqs[i] for i in self._hba]
            elif not hba and lba:
                freqs = [self._freqs[i] for i in self._lba]
            else:
                freqs = self._freqs
            return freqs

    def MAJAX(self, index, sfflag=None, hba=True, lba=True):
        return self._getitem_special(index, self._col_MAJAX, sfflag=sfflag, hba=hba, lba=lba)

    def MINAX(self, index, sfflag=None, hba=True, lba=True):
        return self._getitem_special(index, self._col_MINAX, sfflag=sfflag, hba=hba, lba=lba)

    def PA(self, index, sfflag=None, hba=True, lba=True):
        return self._getitem_special(index, self._col_PA, sfflag=sfflag, hba=hba, lba=lba)

    def SFFLAG(self, index, sfflag=None, hba=True, lba=True):
        return self._getitem_special(index, self._col_SFFLAG, sfflag=sfflag, hba=hba, lba=lba)

    def Sint(self, index, sfflag=None, hba=True, lba=True):
        return self._getitem_special(index, self._col_Sint, sfflag=sfflag, hba=hba, lba=lba)

    def Spk(self, index, sfflag=None, hba=True, lba=True):
        return self._getitem_special(index, self._col_Spk, sfflag=sfflag, hba=hba, lba=lba)



##--------------------------------------------------------- Main function
def gsmMSSS(outfile, ra, dec, radius, cutoff=0.1, assoc=None, patchname=None, cat=None, ndetections=5, verbose=False):
    """
    Return a gms.py-like sky model based on the MSSS catalogue.

    Parameters
    ----------
    outfile : str
        Output filename to store the sky model
    ra : float
        Right ascension of the cone centre to extract (J2000, degrees)
    dec : float
        Declination of the cone centre to extract (J2000, degrees)
    radius : float
        Radius of the sky cone to extract (degrees)
    cutoff : float
        Minimum flux of the MSSS sources to extract (Jy, at 150 MHz)
    assoc : float -----> Not implemented
        Uncertainty in source matching (degrees).
        Default should be 0.00278.
    patchname : str -----> Not implemented
        All sources belong to this single patch
    cat : str
        Filename to read the MSSS catalogue from. By default, will try
        to connect to the MSSS catalogue server.
    ndetections : int
        Minimum number of HBA bands for which the source needs to be detected
        in order to include it in the sky model output.
    """
    if verbose:
        print("-"*80)
        print("Running gsmMSSS at location RA:{0:.5}, DEC:{1:.5}".format(ra,dec))
        print("  Extracting cone of size {0} degres, with flux cutoff {1} Jy".format(radius,cutoff))
        print("  Required minimum number of band detections {}".format(ndetections))
        print("  Using catalogue {0}".format('*default*' if cat is None else cat))
        print("  Output to be saved in file {0}".format(outfile))
        print("")

    ## Patchname functionality is not implemented, display warning
    if assoc is not None:
        warnings.warn("gsmMSSS Warning: The assoc functionality has not been implemented yet. The value you that passed will have no effect.")

    ## Patchname functionality is not implemented, display warning
    if patchname is not None:
        warnings.warn("gsmMSSS Warning: The patchname functionality has not been implemented yet. The value that you passed will have no effect.")

    ## Connecting to the database table if no catalogue filename is provided
    if cat is None:
        try:
            query = vo.scs.SCSQuery('http://vo.astron.nl:8080/msss/q/cone/scs.xml')
            query.pos = (ra, dec)
            query.radius = radius
            msss = TableMSSS.read(query.execute_votable())
        except:
            raise Exception("A problem happened while try to connect to the MSSS database.")
    else:
        msss = TableMSSS.read(cat, format='votable')

    if verbose:
        print("  Total number of sources in the MSSS catalogue: {0}".format(len(msss)))
    ## Returning the list of indices that fall within the sky cone, only needed if we query from a file
    if cat is not None:
        dist = Distance_on_sky(ra, dec, msss['RA'].view(np.ndarray), msss['DEC'].view(np.ndarray), radians=False, one_to_one=False)
        inds = dist <= radius
        ## Cutting down the catalogue to the cone selection
        msss = msss[inds]
        if verbose:
            print("  Total number of sources after selecting the sky cone: {0}".format(len(msss)))

    ## Cutting down the catalogue based on the number of band detections
    inds = msss['NDET_HBA'] >= ndetections
    msss = msss[inds]
    if verbose:
        print("  Total number of sources after applying the number of detections threshold: {0}".format(len(msss)))

    ## Cutting down the catalogue based on the flux cutoff
    inds = msss['A0_HBA'] >= cutoff
    msss = msss[inds]
    if verbose:
        print("  Total number of sources after applying the flux cutoff: {0}".format(len(msss)))

    ## Write the selected sources into a file using the gsm format
    gsmWriter(msss, outfile)

    return



# This is the main entry.
if __name__ == "__main__":

    ## The command-line argument parser
    parser = argparse.ArgumentParser(description='Return a sky model for a given RA/DEC')
    parser.add_argument('outfile', help='Output filename to store the sky model')
    parser.add_argument('RA', type=float, help='Right Ascension of the cone centre to extract (J2000, degrees)')
    parser.add_argument('DEC', type=float, help='Declination of the cone centre to extract (J2000, degrees)')
    parser.add_argument('radius', type=float, help='Radius of the sky cone to extract (degrees)')
    parser.add_argument('-c', '--cutoff', type=float, help='Minimum flux of the MSSS sources to extract (Jy, at 150 MHz). Default 0.1 Jy.', default=0.1)
    parser.add_argument('--cat', type=str, help='Filename to read the MSSS catalogue from. By default, will use msss.xml from the location of that script.')
    parser.add_argument('-n', '--ndetections', type=int, help='Minimum number of HBA bands for which the source needs to be detected in order to include it in the sky model output. Default 5.', default=5)
    parser.add_argument('-v' ,'--verbose', action='count', help='Verbosity')

    args = parser.parse_args()

    gsmMSSS(args.outfile, args.RA, args.DEC, args.radius, cutoff=args.cutoff, cat=args.cat, ndetections=args.ndetections, verbose=args.verbose)







