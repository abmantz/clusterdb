# data source: https://cdsarc.cds.unistra.fr/viz-bin/cat/J/A+A/425/367
f = read.fwf("table1.dat", 
 widths= c(4,-1,12,-1,22,-2,2,-1,2,-1,4,-1,1,2,-1,2,-1,2,-1,7,-1,7,-1,7,-1,6,-1,5,-1,6,-1,6,1,-1,5,-1,10,1),
 col.names=c('prefix', 'CIZA', 'Name', 'RAh', 'RAm', 'RAs', 'DE-', 'DEd', 'DEm', 'DEs', 'GLON', 'GLAT', 'CR', 'e_CR', 'NH', 'FX', 'z', 'u_z', 'LX', 'r_z', 'n_z'),
 colClasses="character",
 strip.white=TRUE)

f$RAdeg = 15.0*(as.numeric(f$RAh) + as.numeric(f$RAm)/60.0 + as.numeric(f$RAs)/3600.0)

sign = rep(1, nrow(f))
sign[which(f[,'DE.']=='-')] = -1

f$DEdeg = sign*(as.numeric(f$DEd) + as.numeric(f$DEm)/60.0 + as.numeric(f$DEs)/3600.0)

write.csv(f, 'ciza.csv', row.names=FALSE)
