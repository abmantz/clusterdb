# data source: https://cdsarc.cds.unistra.fr/viz-bin/cat/J/A+A/425/367
f = read.fwf("reflex70.dat", 
 widths= c(3,12,-2,10,-1,1,-1,2,-1,2,-1,4,-3,1,2,-1,2,-1,2,-1,6,-1,3,-1,7,-1,6,-1,7,-1,5,-1,7,-1,5,-1,3,-1,32),
 col.names=c('prefix', 'RXC', 'Name', 'n_Name', 'RAh', 'RAm', 'RAs', 'DE-', 'DEd', 'DEm', 'DEs', 'z', 'Ngal', 'Flux', 'e_Flux', 'ObsLum', 'Aper', 'LumCor', 'NH', 'Notes', 'r_z'),
 colClasses="character",
 strip.white=TRUE)

f$RAdeg = 15.0*(as.numeric(f$RAh) + as.numeric(f$RAm)/60.0 + as.numeric(f$RAs)/3600.0)

sign = rep(1, nrow(f))
sign[which(f[,'DE.']=='-')] = -1

f$DEdeg = sign*(as.numeric(f$DEd) + as.numeric(f$DEm)/60.0 + as.numeric(f$DEs)/3600.0)

write.csv(f, 'reflex.csv', row.names=FALSE)
