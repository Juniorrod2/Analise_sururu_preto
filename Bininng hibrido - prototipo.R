
hibrid_bucketing <- function(spectra,integrals,bin_width=0.01){
real_spectra <- pepsMatrixToDF(spectra)
integrated_signals <- NMR_integration(real_spectra,integrals)

integrated_regions <- list()
for(i in 1:dim(regioes_integrais)[1]){
  integrated_regions[[i]] <- c(regioes_integrais[[i,2]],regioes_integrais[[i,3]])
}
names(integrated_regions) <- regioes_integrais[[1]]

unintegrated_spectra <- RegionRemoval(spec,fromto.rr = integrated_regions,
                                      typeofspectra = "manual")

unintegrated_spectra <- Bucketing(unintegrated_spectra,width = T,mb=bin_width)
unintegrated_spectra <- pepsMatrixToDF(unintegrated_spectra)

unintegrated_spectra <- CleanDataMatrix(as.matrix(unintegrated_spectra[-1]))

hibrid_binning_data <- cbind(integrated_signals,unintegrated_spectra)
}

# ------------------------------------------------------------------------------

sururu_data <- tidyr::separate(hibrid_binning_data[c(-1,-10:-15),],Sample,into=c("Group"),sep = " - ",remove = F)


