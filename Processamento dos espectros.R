library(PepsNMR); library(NAPRMN)


Sururu_fid <- ReadFids("Espectros sururu preto")

fid <- Sururu_fid$Fid_data
info <- Sururu_fid$Fid_info

fid_proc <- GroupDelayCorrection(fid,info)

fid_proc <- SolventSuppression(fid_proc)

fid_proc <- Apodization(fid_proc,info)

fid_proc <- ZeroFilling(fid_proc) 

spec <- FourierTransform(fid_proc,info)

spec <- ZeroOrderPhaseCorrection(spec)


spec <- InternalReferencing(spec,info)


spec <- BaselineCorrection(spec)

spec <- NegativeValuesZeroing(spec)


spec <- WindowSelection(spec,from.ws = 9.2,to.ws = 0.5)

saveRDS(spec,"Sururu_spectra")

fullres <- pepsMatrixToDF(spec)


saveRDS(fullres,"Sururu_preto_espectro_fullRes")

spec_binning <- Bucketing(spec,width = TRUE,
                          mb=0.01,intmeth = "t")


bins <- pepsMatrixToDF(spec_binning)

writexl::write_xlsx(bins,"bins_sururu.xlsx")
