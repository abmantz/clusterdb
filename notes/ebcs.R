# data source: https://cdsarc.cds.unistra.fr/viz-bin/cat/J/MNRAS/318/333
f = read.fwf("table1.dat.gz", 
 widths= c(5,-1,14,-1,7,-1,7,-1,4,-1,4,-1,4,-1,3,-1,4,-1,4,-1,4,-1,6,-1,3,-1,5,-2,2),
 col.names=c('Notes', 'Name', 'RAdeg', 'DEdeg', 'nH20', 'Texp', 'CRVTP', 'RadVTP', 'CR', 'e_CR', 'kT', 'z', 'FX', 'LX', 'r_z'),
 colClasses="character",
 strip.white=TRUE)

write.csv(f, 'ebcs.csv', row.names=FALSE)
