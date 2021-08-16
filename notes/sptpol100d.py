import astropy.io.fits as pyfits
import numpy as np
from pandas import DataFrame

# https://pole.uchicago.edu/public/data/sptsz-clusters/
fits = pyfits.open('sptpol100d_catalog_huang19.fits.gz')
cat = fits[1].data

def criticalDensity(z):
    return 1.360e+11 + 1.224e+11*z + 1.224e+11*z**2 + 4.080e+10*z**3
def comovingDistance(x):
    return 4272.3498*x - 1093.4370*x**2 + 126.0632*x**3
def dA(z):
    return comovingDistance(z) / (1.0+z)

f = DataFrame({'name':cat['SPT_ID'],
           'RA':cat['RA'],
           'dec':cat['Dec'],
           'redshift':np.where(cat['redshift']>0.0, cat['redshift'], np.nan),
           'r500':np.where(cat['redshift']>0.0, (cat['M500']*1e14 / (4./3. * np.pi * 500. * criticalDensity(cat['redshift'])))**0.3333 / dA(cat['redshift']) * 180./np.pi, np.nan)
          })

f.to_csv('sptpol100d.csv.gz', index=False)
