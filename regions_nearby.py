import argparse
from astropy.coordinates import SkyCoord
import astropy.units as u
from glob import glob
import os
import numpy as np
import pandas as pd
import yaml

# get the root of the code repo
thisdir = os.path.dirname(os.path.realpath(__file__)) + '/'

# class that knows how to print output
class Match:
    def __init__(self, coords, distance, name=None, radius=None, redshift=None, bonus=None, bonus_col=None):
        self.coords = coords
        self.radius = radius
        self.distance = distance
        self.name = name
        self.redshift = redshift
        if bonus_col is not None and bonus is not None:
            self.bonus = ' '+bonus_col+'='+str(bonus)
        else:
            self.bonus = ''
    def __str__(self):
        if self.radius is None or not np.isfinite(self.radius):
            s = 'point('
        else:
            s = 'circle('
        s += str(self.coords.ra.degree) + "," + str(self.coords.dec.degree)
        if self.radius is not None and np.isfinite(self.radius):
            s += "," + str(self.radius.to_value(u.arcmin)*radius_boost) + "'"
        s += ") # "
        if self.name is not None:
            s += "text={" + self.name + '} '
        if self.redshift is not None:
            s += 'z=' + str(self.redshift) + ' '
        s += 'distance: ' + str(self.distance.to_value(u.arcmin)) + "'"
        s += self.bonus
        return s
    def notregion(self):
        if self.name is None:
            s = ""
        else:
            s = self.name
        s += ' | ' + str(self.coords.ra.degree) + " | " + str(self.coords.dec.degree) + ' | ' + str(self.distance.to_value(u.arcmin)) + ' | '
        if self.radius is not None:
            s += str(self.radius.to_value(u.arcmin)*radius_boost)
        s += ' | '
        if self.redshift is not None:
            s += str(self.redshift)
        return s

# class to read catalogs from a csv
class generic_from_csv:
    def __init__(self, *args, **kwargs):
        self.data = self._read(*args, **kwargs)
    def _read(self, fname, ra_col, dec_col, catname='', z_col=None, name_col=None, pos_err_col=None, r500_col=None, coord_unit='deg', pos_err_unit=u.degree, r500_unit=u.degree, bonus_col=None):
        x = pd.read_csv(thisdir+fname)
        self.pos = SkyCoord(ra=x[ra_col], dec=x[dec_col], unit=coord_unit)
        if z_col is None:
            self.z = [None] * len(self.pos)
        else:
            self.z = x[z_col]
        if pos_err_col is None:
            self.pos_err = [None] * len(self.pos)
        else:
            self.pos_err = np.array(x[pos_err_col]) * pos_err_unit
        self.r500 = [None] * len(self.pos)
        if r500_col is None:
            pass
        else:
            try:
                self.r500 = np.array(x[r500_col]) * r500_unit
            except:
                pass
        if name_col is None:
            self.name = ["???"+catname] * len(self.pos)
        else:
            self.name = [str(n)+catname for n in x[name_col]]
        self.catname = catname
        if bonus_col is None:
            self.bonus = [None] * len(self.pos)
        else:
            self.bonus = x[bonus_col]
        self.bonus_col = bonus_col
        return x
    def search(self, from_coords, radius):
        d = from_coords.separation(self.pos)
        selection = np.where(d<=radius)[0]
        return [Match(self.pos[i], d[i], name=self.name[i], radius=self.r500[i], redshift=self.z[i], bonus=self.bonus[i], bonus_col=self.bonus_col) for i in selection[np.argsort(d[selection])]]

# specialization to deal with a catalog with multiple name columns
class MCXC(generic_from_csv):
    def __init__(self, *args, **kwargs):
        generic_from_csv.__init__(self, *args, **kwargs)
        x = self.data
        for i in range(len(self.name)):
            if str(x['Aname'][i]) != "nan":
                self.name[i] = str(x['Aname'][i])+self.catname
            elif str(x['Oname'][i]) != "nan":
                self.name[i] = str(x['Oname'][i])+self.catname

# these are the kinds of catalog classes we know of
classes = {'generic_from_csv':generic_from_csv, 'MCXC':MCXC}

# read properties of built-in catalogs
with open(thisdir+'data/catalogs.yaml', 'r') as thisfile:
    info = yaml.safe_load(thisfile.read())
# same for local catalogs, if any, overwriting any common entries
try:
    with open(thisdir+'local/catalogs.yaml', 'r') as thisfile:
        loc = yaml.safe_load(thisfile.read())
        info.update(loc)
        local_cats = list(loc.keys())
except FileNotFoundError:
    local_cats = []
    pass

# these are the built-in catalogs used by default, plus locals, sorted
default_str = 'MCXC PSZ2 SPTSZ SPTECS SPTPOL100d AdvancedACTPol SDSS_RM'
if not 'DES_Y3' in local_cats:
    default_str += 'DES_Y1'
default_cats = list(set(default_str.split() + local_cats))
default_cats.sort()
default_cats = ' '.join(default_cats)

# list of all catalogs
all_cats = list(info.keys())
all_cats.sort()


# argument parsing
parser = argparse.ArgumentParser(description="Look up known clusters within some angular distance of a given position.\nDefault catalogs: "+default_cats)
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
    help="search radius in arcmin (default=15)",
    default=15.0
)
parser.add_argument(
    '--cats',
    help="catalogs to search, a subset of " + ' '.join(all_cats),
    default=default_cats
)
parser.add_argument(
    '--old-xray',
    dest="oldxray",
    help="add BCS, eBCS, REFLEX, CIZA, CIZA2 to the list of cats",
    action="store_true"
)
parser.add_argument(
    '--notregion',
    dest="notregion",
    help="write in non-region format",
    action="store_true"
)
parser.add_argument(
    '--r200',
    dest="r200",
    help="use r200 instead of r500 for region radii",
    action="store_true"
)

args = parser.parse_args()
radius = args.radius * u.arcmin
radius_boost = 1.0
radius_name = 'r500'
if args.r200:
    radius_boost = 1.53
    radius_name = 'r200'


# where to search
target = SkyCoord(ra=args.ra, dec=args.dec, unit='deg')


# finalize the list of catalogs to look in
catss = args.cats
if args.oldxray:
    catss += ' BCS eBCS REFLEX CIZA CIZA2'
catsl = list(set(catss.split())) # eliminate repeats

# read the catalogs
cats = [classes[info[k].pop('class')](**info[k]) for k in catsl]


# write an appropriate header
if args.notregion:
    print("name | ra | dec | distance[arcmin] | "+radius_name+"[arcmin] | redshift")
else:
    print('# Region file format: DS9 version 4.1')
    print('global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1')
    print('fk5')

# write matches found in each catalog
for cat in cats:
    for m in cat.search(target, radius):
        if args.notregion:
            print(m.notregion())
        else:
            print(str(m))
