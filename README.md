# clusterdb

Look up galaxy clusters near a given position on the sky. Output is formatted as a region file for [DS9](https://sites.google.com/cfa.harvard.edu/saoimageds9).
Matching sources will be either points or circles (depending on the information in a given catalog), with the source name as a label.
The angular distance between the search location and source position appears as a comment.

## Running

`python regions_nearby.py -h` gives instructions. Output goes to the terminal.

## Dependencies
* astropy
* numpy
* pandas

## Implemented catalogs
* ROSAT Brightest Cluster Survey and Extended Brightest Cluster Survey: (e)BCS - [Ebeling et al. 1998](http://adsabs.harvard.edu/cgi-bin/nph-bib_query?bibcode=1998MNRAS.301..881E&db_key=AST) and [2000](http://adsabs.harvard.edu/abs/2000MNRAS.318..333E)
* ROSAT-ESO Flux Limited X-ray (REFLEX): [Bohringer et al. 2004](http://adsabs.harvard.edu/cgi-bin/nph-bib_query?bibcode=2004A%26A...425..367B&db_key=AST)
* Clusters In the Zone of Avoidance (CIZA I and II): [Ebeing et al. 2002](http://adsabs.harvard.edu/abs/2002ApJ...580..774E) and [Kocevski et al. 2007](http://adsabs.harvard.edu/abs/2007ApJ...662..224K)

Note: the catalogs above are in principle included in the MCXC, and are only searched separately from it if the `--old-xray` flag is used.

* Meta-Catalog of X-ray detected Clusters (MCXC): [Piffaretti et al. 2011](http://adsabs.harvard.edu/abs/2011A%26A...534A.109P)
* Second Planck catalog of Sunyaev-Zel'dovich sources (PSZ2): [Planck Collaboration 2016](http://adsabs.harvard.edu/abs/2016A%26A...594A..27P)
* South Pole Telescope
  * 2500 square degree SPT-SZ survey: [Bleem et al. 2015](http://adsabs.harvard.edu/abs/2015ApJS..216...27B)
  * SPTPol Extended Cluster Survey (SPTECS): [Bleem et al. 2020](https://ui.adsabs.harvard.edu/abs/2020ApJS..247...25B)
  * SPTPol 100 square degree: [Huang et al. 2020](https://ui.adsabs.harvard.edu/abs/2020AJ....159..110H)

## Adding a new catalog

One can either write a new class to deal with a new data format, or inherit the existing class, which reads a csv file. The only required data are
* positions in a format that `astropy.coordinates` can interpret
* source names

Optionally, one can include a radius (in angular units), in which case the corresponding region written will be a circle. Otherwise, the region will be a point.
