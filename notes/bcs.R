# data source: https://cdsarc.cds.unistra.fr/viz-bin/cat/J/MNRAS/301/881
f = read.fwf("table3.dat.gz", 
 widths= c(4,-1,15,-1,2,-1,2,-1,4,-1,1,2,-1,2,-1,2,-1,7,-1,7,-1,4,-1,4,-1,5,-1,4,-1,6,-1,5,-1,4,1,-1,5,-1,5,-1,6,-1,2),
 col.names=c('n_Name', 'Name', 'RAh', 'RAm', 'RAs', 'DE-', 'DEd', 'DEm', 'DEs', 'RAdeg', 'DEdeg', 'NH', 'Texp', 'CR(VTP)',
            'R(VTP)', 'CR', 'e_CR', 'kT', 'n_kT', 'z', 'FX', 'LX', 'r_z'),
 colClasses="character",
 strip.white=TRUE)

write.csv(f, 'bcs.csv', row.names=FALSE)
