import astropy.io.fits as pyfits
import numpy as np
from pandas import DataFrame

# https://irsa.ipac.caltech.edu/data/Planck/release_2/catalogs/
#  -> Sunyaev-Zeldovich [SZ] Catalog
union = pyfits.open('HFI_PCCS_SZ-union_R2.08.fits.gz')
unioncat = union[1].data

def criticalDensity(z):
    return 1.360e+11 + 1.224e+11*z + 1.224e+11*z**2 + 4.080e+10*z**3
def comovingDistance(x):
    return 4272.3498*x - 1093.4370*x**2 + 126.0632*x**3
def dA(z):
    return comovingDistance(z) / (1.0+z)

f = DataFrame({'INDEX':unioncat['INDEX'],
           'name':unioncat['NAME'],
           'RA':unioncat['RA'],
           'dec':unioncat['DEC'],
           'redshift':np.where(unioncat['REDSHIFT']>0.0, unioncat['REDSHIFT'], np.nan),
           'pos_err':unioncat['POS_ERR'] / 60.,
           'r500':np.where(unioncat['REDSHIFT']>0.0, (unioncat['MSZ']*1e14 / (4./3. * np.pi * 500. * criticalDensity(unioncat['REDSHIFT'])))**0.3333 / dA(unioncat['REDSHIFT']) * 180./np.pi, np.nan)
          })

f.to_csv('psz2.csv.gz', index=False)
