
spec <- readRDS("Sururu_spectra")

sururu_spectra <- readRDS("Sururu_preto_espectro_fullRes")

sururu_data <- cbind("group"=c(rep("P4C1",8),rep("P6C3",5)),sururu_spectra)

sururu_bins <- cbind("group"=c(rep("P4C1",8),rep("P6C3",5)),bins)


library(ropls)

PLS <- opls(sururu_data[-1:-2],sururu_data$group,predI = 2)

PLS_data <- extract_ropls_data(PLS)

library(ggplot2)

ggplot(PLS_data$Scores,aes(p1,p2,color=Group)) + geom_point(size=3) + 
  stat_ellipse(geom="polygon",aes(fill=Group),alpha=0.2)+theme_bw() + 
  geom_text(aes(label=Samples),nudge_x = 20)

loading <- ggplot(PLS_data$Loadings,aes(p1,p2,label=bins,color=as.character(Vip)))+geom_text()
plotly::ggplotly(loading)


NAPRMN::plot_interactive_Spectra(spec)
