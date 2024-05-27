#inicializacao dos pacotes
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

# Metado de alinhamento via PTW
#align_spec2 <- Warping(spec,ptw.wp = T,optim.crit = "WCC",verbose = T)

spec <- BaselineCorrection(spec)

spec <- NegativeValuesZeroing(spec)

plot_interactive_Spectra(spec)

#align_spec2 <- Warping(spec,ptw.wp = T,optim.crit = "WCC",verbose = T)

plot_interactive_Spectra(align_spec2)

y.spec <- WindowSelection(Y,from.ws = 9.2,to.ws = 0.5)

#Gera a previsualizacao do espectro, pode ser chamada apos qualquer um dos intermediarios
#acima. Conferir a documentacao da funcao para detalhes da representacao do espectro
#e parametros opcionais
plot_interactive_Spectra(y.spec)

# Salva o objeto em um arquivo binario (no computador) que pode ser lido no futuro
saveRDS(spec,"Sururu_spectra")

# Converte a matrix em um dataframe e converte os numeros complexos em numeros reais
fullres <- pepsMatrixToDF(spec)

# Salva o objeto em um arquivo binario (no computador) que pode ser lido no futuro
saveRDS(fullres,"Sururu_preto_espectro_fullRes")

#Produz o binning dos espectros
spec_binning <- Bucketing(y.spec,width = TRUE,
                          mb=0.01,intmeth = "t")

# Converte a matrix em um dataframe e converte os numeros complexos em numeros reais
bins <- pepsMatrixToDF(spec_binning)

#Salva/grava o dataframe com os bins na forma de uma planilha do excel (xlsx)
#o "writexl::" indica o namespace da funcao e substutui um "library(writexl)" que
#produziria o mesmo resultado 
writexl::write_xlsx(bins,"bins_sururu.xlsx")
