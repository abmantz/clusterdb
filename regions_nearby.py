import argparse
from astropy.coordinates import SkyCoord
import astropy.units as u
from glob import glob
import numpy as np
import pandas as pd

class Match:
    def __init__(self, coords, distance, name=None, radius=None, redshift=None):
        self.coords = coords
        self.radius = radius
        self.distance = distance
        self.name = name
        self.redshift = redshift
    def __str__(self):
        if self.radius is None or not np.isfinite(self.radius):
            s = 'point('
        else:
            s = 'circle('
        s += str(self.coords.ra.degree) + "," + str(self.coords.dec.degree)
        if self.radius is not None and np.isfinite(self.radius):
            s += "," + str(self.radius.to_value(u.degree))
        s += ") # "
        if self.name is not None:
            s += "text={" + self.name + '} '
        if self.redshift is not None:
            s += 'z=' + str(self.redshift) + ' '
        s += 'distance: ' + str(self.distance.to(u.arcmin))
        return s

class generic_from_csv:
    def __init__(self):
        # create required fields, not that things will function without being
        # properly constructed
        self.pos = None
        self.z = None
        self.pos_err = None
        self.r500 = None
        self.name = None
    def read(self, fname, ra_col, dec_col, catname='', z_col=None, name_col=None, pos_err_col=None, r500_col=None, coord_unit='deg', pos_err_unit=u.degree, r500_unit=u.degree):
        x = pd.read_csv(fname)
        self.pos = SkyCoord(ra=x[ra_col], dec=x[dec_col], unit=coord_unit)
        if z_col is None:
            self.z = [None] * len(self.pos)
        else:
            self.z = x[z_col]
        if pos_err_col is None:
            self.pos_err = [None] * len(self.pos)
        else:
            self.pos_err = np.array(x[pos_err_col]) * pos_err_unit
        if r500_col is None:
            self.r500 = [None] * len(self.pos)
        else:
            self.r500 = np.array(x[r500_col]) * r500_unit
        if name_col is None:
            self.name = ["???"+catname] * len(self.pos)
        else:
            self.name = [n+catname for n in x[name_col]]
        return x
    def search(self, from_coords, radius):
        d = from_coords.separation(self.pos)
        selection = np.where(d<=radius)[0]
        return [Match(self.pos[i], d[i], name=self.name[i], radius=self.r500[i], redshift=self.z[i]) for i in selection[np.argsort(d[selection])]]

class BCS(generic_from_csv):
    def __init__(self):
        generic_from_csv.__init__(self) # unnecessarily
        self.read('data/bcs.csv.gz', ra_col='RAdeg', dec_col='DEdeg', catname=' (BCS)', z_col='z', name_col='Name')

class eBCS(generic_from_csv):
    def __init__(self):
        generic_from_csv.__init__(self) # unnecessarily
        self.read('data/ebcs.csv.gz', ra_col='RAdeg', dec_col='DEdeg', catname=' (eBCS)', z_col='z', name_col='Name')

class REFLEX(generic_from_csv):
    def __init__(self):
        generic_from_csv.__init__(self) # unnecessarily
        self.read('data/reflex.csv.gz', ra_col='RAdeg', dec_col='DEdeg', catname=' (REFLEX)', z_col='z', name_col='RXC')

class CIZA(generic_from_csv):
    def __init__(self):
        generic_from_csv.__init__(self) # unnecessarily
        self.read('data/ciza.csv.gz', ra_col='RAdeg', dec_col='DEdeg', catname=' (CIZA)', z_col='z', name_col='CIZA')
        
class CIZA2(generic_from_csv):
    def __init__(self):
        generic_from_csv.__init__(self) # unnecessarily
        self.read('data/ciza2.csv.gz', ra_col='RAdeg', dec_col='DEdeg', catname=' (CIZA2)', z_col='z', name_col='CIZA')
        
class MCXC(generic_from_csv):
    def __init__(self):
        generic_from_csv.__init__(self) # unnecessarily
        catname = ' (MCXC)'
        x = self.read('data/mcxc.csv.gz', ra_col='RAdeg', dec_col='DEdeg', catname=catname, z_col='z', name_col='MCXC', r500_col='r500')
        for i in range(len(self.name)):
            if x['Aname'][i] != "":
                self.name[i] = str(x['Aname'][i])+catname
            elif x['Oname'][i] != "":
                self.name[i] = str(x['Oname'][i])+catname

class PSZ2(generic_from_csv):
    def __init__(self):
        generic_from_csv.__init__(self) # unnecessarily
        self.read('data/psz2.csv.gz', ra_col='RA', dec_col='dec', catname='', z_col='redshift', name_col='name', pos_err_col='pos_err', r500_col='r500')

class SPTSZ(generic_from_csv):
    def __init__(self):
        generic_from_csv.__init__(self) # unnecessarily
        self.read('data/sptsz.csv.gz', ra_col='RA', dec_col='dec', catname='', z_col='redshift', name_col='name', r500_col='r500')

class SPTECS(generic_from_csv):
    def __init__(self):
        generic_from_csv.__init__(self) # unnecessarily
        self.read('data/sptecs.csv.gz', ra_col='RA', dec_col='dec', catname='', z_col='redshift', name_col='name', r500_col='r500')

class SPTPOL100d(generic_from_csv):
    def __init__(self):
        generic_from_csv.__init__(self) # unnecessarily
        self.read('data/sptpol100d.csv.gz', ra_col='RA', dec_col='dec', catname='', z_col='redshift', name_col='name', r500_col='r500')

class AdvancedACTPol(generic_from_csv):
    def __init__(self):
        generic_from_csv.__init__(self) # unnecessarily
        self.read('data/aactpol.csv.gz', ra_col='RA', dec_col='dec', catname='', z_col='redshift', name_col='name', r500_col='r500')

class YetAnotherClusterSurvey(generic_from_csv):
    def __init__(self, csvname):
        generic_from_csv.__init__(self) # unnecessarily
        self.read(csvname, ra_col='RA', dec_col='dec', catname='', z_col='redshift', name_col='name')

        
parser = argparse.ArgumentParser(description="Look up known clusters within some angular distance of a given position.")
parser.add_argument(
    'ra',
    type=float,
    help='J2000 right ascension in degrees'
)
parser.add_argument(
    'dec',
    type=float,
    help='J2000 declination in degrees'
)
parser.add_argument(
    '--radius',
    type=float,
    help="search radius in arcmin (default=15')",
    default=15.0
)
parser.add_argument(
    '--old-xray',
    dest="oldxray",
    help="check BCS, REFLEX, etc. in addition to MCXC",
    action="store_true"
)


args = parser.parse_args()
radius = args.radius * u.arcmin

target = SkyCoord(ra=args.ra, dec=args.dec, unit='deg')

cats = [MCXC(), PSZ2(), SPTSZ(), SPTECS(), SPTPOL100d(), AdvancedACTPol()]
if args.oldxray:
    cats += [BCS(), eBCS(), REFLEX(), CIZA(), CIZA2()]
for f in glob('local/*.csv.gz'):
    cats += [YetAnotherClusterSurvey(f)]

print('# Region file format: DS9 version 4.1')
print('global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1')
print('fk5')

for cat in cats:
    for m in cat.search(target, radius):
        print(str(m))
