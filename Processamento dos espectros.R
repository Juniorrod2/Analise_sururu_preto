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

# Metodo de alinhamento via PTW
#align_spec2 <- Warping(spec,ptw.wp = T,optim.crit = "WCC",verbose = T)

#alinhamento
library(speaq)

spec_peaks <- detectSpecPeaks(Re(spec),baselineThresh = 5E6)

refInd<- findRef(spec_peaks)$refInd

spec_aling <- dohCluster(Re(spec),peakList = spec_peaks,
                         refInd = refInd,maxShift = 100)

plot_interactive_Spectra(spec_aling)

spec <- BaselineCorrection(spec_aling)

spec <- NegativeValuesZeroing(spec_aling)

plot_interactive_Spectra(spec)

spec <- WindowSelection(spec,from.ws = 9.2,to.ws = 0.5)

spec <- RegionRemoval(spec,fromto.rr = list("water"=c(4.75,5.2)))
#Produz o binning dos espectros
spec_binning <- Bucketing(spec,width = TRUE,
                          mb=0.01,intmeth = "t")

# Converte a matrix em um dataframe e converte os numeros complexos em numeros reais
bins <- pepsMatrixToDF(spec_binning)

bins <- tidyr::separate(bins,Sample,into=c("Group"),sep = " - ",remove = F)

#bins <- dplyr::filter(bins,Group%in%c("P8C1","P5C2","P6C3"))
#bins <- dplyr::filter(bins,!Group%in%c("P4C2","P4C2- 2","P4C2- 3","P4C2- 4","P4C2- 5","P4C2- 6"))
bins <- dplyr::filter(bins,Group%in%c("P8C1","P8C3"))

library(ropls);library(ggplot2);library(ggrepel);library(patchwork);library(forcats)

regioes_integrais <- readxl::read_excel("regioes_integrais.xlsx")

norm_bins <- CleanDataMatrix(bins[-1:-2])
norm_bins <- DataNormalization(norm_bins,rowNorm = "MedianNorm",scaleNorm = "ParetoNorm")
bins <- cbind(bins[1:2],norm_bins)

pls_sururu <- opls(bins[-1:-2],bins$Group,predI=3)

pls_sururu_data <- extract_ropls_data(pls_sururu)

gplot <- ggplot(pls_sururu_data$Loadings,aes(p1,p2,label=bins,color=as.character(Vip)))+geom_text()
plotly::ggplotly(gplot)


sururu_hbins <- hibrid_bucketing(spectra = spec,integrals = regioes_integrais,bin_width = 0.01)

hbins <- tidyr::separate(sururu_hbins,Sample,into=c("Group"),sep = " - ",remove = F)

hbins <- dplyr::filter(hbins,Group%in%c("P8C1","P5C2","P6C3"))

norm_hbins <- CleanDataMatrix(hbins[-1:-2])
norm_hbins <- DataNormalization(norm_hbins,rowNorm = "MedianNorm",scaleNorm = "ParetoNorm")
hbins <- cbind(hbins[1:2],norm_hbins)

identified_metabolites <- hbins[1:44]


pca_sururu <- opls(hbins[-1:-2])

pca_sururu_data <- extract_ropls_data(pca_sururu)

pca_sururu_data$Loadings$bins <- rownames(pca_sururu_data$Loadings)

pca_sururu_loadings <- pca_sururu_data$Loadings[1:42,]

pca_sururu_scores <- cbind("Group"=hbins$Group,pca_sururu_data$Scores)

scores <- ggplot(pca_sururu_scores,aes(as.numeric(p1),as.numeric(p2),color=Group,fill=Group,label=Samples))+
  geom_point(size=4)+#geom_text(nudge_x=0.2,size=5)+
  stat_ellipse(geom="polygon",alpha=0.2)+theme_bw()+
  scale_color_manual(breaks=c("P6C3","P5C2","P8C1"),labels=c("Hg 0,01 (mg/kg)","Hg 0,16 (mg/kg)","Hg 0,35 (mg/kg)"),
                     values=c("#30E50C","#E5A30C","#FB0303"))+
  scale_fill_manual(breaks=c("P6C3","P5C2","P8C1"),labels=c("Hg 0,01 (mg/kg)","Hg 0,16 (mg/kg)","Hg 0,35 (mg/kg)"),
                    values=c("#30E50C","#E5A30C","#FB0303"))+
  labs(color="Group",fill="Group",x="Comp 1 (41%)",y="Comp 1 (17%)",title="PCA (M. charruana)")

loading <- ggplot(pca_sururu_loadings,aes(as.numeric(p1),as.numeric(p2),label=bins)) +
  geom_text_repel(seed = 123,fontface="bold",size=5,force_pull = 0.8)+geom_point(color="black",size=3)+
  theme_bw()+scale_color_gradient(low="gray",high="red")+#xlim(-0.25,0.6)+
  labs(x="Comp 1 (41%)",y="Comp 1 (17%)")

png("PCA - Sururu.png",res = 300,height = 3600,width = 6200,units = "px")
scores+loading
dev.off()

pls_sururu <- opls(hbins[-1:-2],hbins$Group,predI=3)

pls_sururu_data <- extract_ropls_data(pls_sururu)

pls_sururu_data$Loadings$bins <- rownames(pls_sururu_data$Loadings)

sururu_loadings <- pls_sururu_data$Loadings[1:42,]

scores <- ggplot(pls_sururu_data$Scores,aes(as.numeric(p1),as.numeric(p2),color=Group,fill=Group,label=Samples))+
  geom_point(size=4)+#geom_text(nudge_x=0.2,size=5)+
  stat_ellipse(geom="polygon",alpha=0.2)+theme_bw()+
  scale_color_manual(breaks=c("P6C3","P5C2","P8C1"),labels=c("Hg 0,01 (mg/kg)","Hg 0,16 (mg/kg)","Hg 0,35 (mg/kg)"),
                     values=c("#30E50C","#E5A30C","#FB0303"))+
  scale_fill_manual(breaks=c("P6C3","P5C2","P8C1"),labels=c("Hg 0,01 (mg/kg)","Hg 0,16 (mg/kg)","Hg 0,35 (mg/kg)"),
                    values=c("#30E50C","#E5A30C","#FB0303"))+
  labs(color="Group",fill="Group",x="Comp 1 (41%)",y="Comp 1 (16%)",title="PLS-DA (M. charruana)")
  
loading <- ggplot(sururu_loadings,aes(as.numeric(p1),as.numeric(p2),label=bins,color=Vip)) +
  geom_text_repel(seed = 123,fontface="bold",size=5,force_pull = 0.8)+geom_point(color="black",size=3)+
  theme_bw()+scale_color_gradient(low="gray",high="red")+#xlim(-0.25,0.6)+
  labs(x="Comp 1 (41%)",y="Comp 1 (16%)",color="V.I.P.")

png("PLS-DA - Sururu.png",res = 300,height = 3600,width = 6200,units = "px")
scores+loading
dev.off()

library(ggpubr);library(tidyr)

svg("Sururu-anova.svg",height = 12,width = 18)
gather(identified_metabolites,-Sample,-Group,key="Metabolites",value="Intensities")%>%
  mutate(Group=fct_relevel(Group,c("P6C3","P5C2","P8C1")))%>%
  ggboxplot(x="Group",y="Intensities",color="Group")+facet_wrap(~Metabolites,scales = "free")+
  stat_compare_means(vjust = -2)+geom_pwc(hide.ns = T,label = "p.signif",step.increase = 0.35,bracket.nudge.y = 0.1,vjust = 0.1)+
  scale_y_continuous(expand = expansion(mult = c(0.05,0.4)))+
  scale_color_manual(values=c("#30E50C","#E5A30C","#FB0303"),labels=c("Hg 0,01 (mg/kg)","Hg 0,16 (mg/kg)","Hg 0,35 (mg/kg)"))+theme_bw()
dev.off()
