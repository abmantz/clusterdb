import argparse
from astropy.coordinates import SkyCoord
import astropy.units as u
import numpy as np
import pandas as pd

class Match:
    def __init__(self, coords, distance, name=None, radius=None):
        self.coords = coords
        self.radius = radius
        self.distance = distance
        self.name = name
    def __str__(self):
        if self.radius is None:
            s = 'point('
        else:
            s = 'circle('
        s += str(self.coords.ra.degree) + "," + str(self.coords.dec.degree)
        if self.radius is not None:
            s += "," + str(self.radius.to_value(u.degree))
        s += ") # "
        if self.name is not None:
            s += "text={" + self.name + '} '
        s += 'distance: ' + str(self.distance.to(u.arcmin))
        return s

class PSZ2:
    def __init__(self):
        x = pd.read_csv('data/psz2.csv.gz')
        self.pos = SkyCoord(ra=x['RA'], dec=x['dec'], unit='deg')
        self.z = x['redshift']
        self.pos_err = x['pos_err'] * u.degree
        self.r500 = x['r500'] * u.degree
        self.name = x['name']
    def search(self, from_coords, radius):
        d = from_coords.separation(self.pos)
        selection = np.where(d<=radius)[0]
        return [Match(self.pos[i], d[i], name=self.name[i], radius=self.r500[i]) for i in selection[np.argsort(d[selection])]]


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


args = parser.parse_args()
radius = args.radius * u.arcmin

target = SkyCoord(ra=args.ra, dec=args.dec, unit='deg')

print('# Region file format: DS9 version 4.1')
print('global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1')
print('fk5')

for cat in [PSZ2()]:
    for m in cat.search(target, radius):
        print(str(m))
