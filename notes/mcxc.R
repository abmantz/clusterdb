# data source: https://cdsarc.cds.unistra.fr/viz-bin/cat/J/A+A/534/A109
f = read.fwf("mcxc.dat", 
 widths= c(12,-1,18,-1,54,-1,2,-1,2,-1,4,-1,1,2,-1,2,-1,2,-1,7,-1,7,-1,7,-1,7,-1,6,-1,12,-1,12,-1,5,-1,9,-1,7,-1,7,-1,42,-1,12,-1,12,-1,12,-1,12,-1,5,-1,5,-1,5,-1,5),
 col.names=c('MCXC', 'Oname', 'Aname', 'RAh', 'RAm', 'RAs', 'DE-', 'DEd', 'DEm', 'DEs', 'RAdeg', 'DEdeg', 'GLON', 'GLAT', 'z', 'Cat', 'Sub-Cat', 'Scale', 'L500', 'M500', 'R500', 'Notes', 'Cat1', 'Cat2', 'Cat3', 'Cat4', 'L500r1', 'L500r2', 'L500r3', 'L500r4'),
 colClasses="character",
 strip.white=TRUE)

f$r500 = as.numeric(f$R500)*1e3 / as.numeric(f$Scale) / 3600

write.csv(f, 'mcxc.csv', row.names=FALSE)
