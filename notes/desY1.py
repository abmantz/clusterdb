import astropy.io.fits as pyfits
import numpy as np
from pandas import DataFrame

fits = pyfits.open('redmapper_y1a1_public_v6.4_catalog.fits.gz')
cat = fits[1].data

def criticalDensity(z):
    return 1.360e+11 + 1.224e+11*z + 1.224e+11*z**2 + 4.080e+10*z**3
def comovingDistance(x):
    return 4272.3498*x - 1093.4370*x**2 + 126.0632*x**3
def dA(z):
    return comovingDistance(z) / (1.0+z)

f = DataFrame({'ID':cat['ID'],
           'RA':cat['RA'],
           'dec':cat['DEC'],
           'redshift':np.where(cat['Z_LAMBDA']>0.0, cat['Z_LAMBDA'], np.nan),
           'lambda': cat['LAMBDA'],
           'r500':np.where(cat['Z_LAMBDA']>0.0, ( (cat['LAMBDA']*np.exp(-4.80))**(1./0.73)*8e14 / (4./3. * np.pi * 500. * criticalDensity(cat['Z_LAMBDA'])))**0.3333 / dA(cat['Z_LAMBDA']) * 180./np.pi, np.nan)
          })

f.to_csv('desY1a1.csv.gz', index=False)
