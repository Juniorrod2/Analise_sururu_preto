bins_sururu <- cbind("Group"=c(rep("P4C1",8),rep("P6C3",5)),bins)


library(ropls)

pca_sururu <- opls(bins_sururu[,c(-1:-2)],predI=4)

pls_sururu <- opls(bins_sururu[,c(-1:-2)],bins_sururu$Group,predI=3,orthoI=1)



plot(pca_sururu,parCompVi=c(1,2),parAsColFcVn=bins_sururu$Group)

plot(pls_sururu,typeVc=c("x-score","x-loading"))

pls_sururu_data <- extract_ropls_data(pls_sururu)

library(ggplot2)

scores <- ggplot(pls_sururu_data$Scores,aes(p1,o1,color=Group,fill = Group))+
  geom_point(size=3)+stat_ellipse(geom="polygon",alpha=0.4)+theme_bw()+
  labs(x="Predictive component",y="Orthogonal component",title = "OPLS-DA")


loadings <- ggplot(pls_sururu_data$Loadings,aes(p1,o1,color=as.character(Vip),label=bins)) +
  geom_text()

library(plotly)

ggplotly(plot)

library(patchwork)

scores+loadings

spectra <- pepsMatrixToDF(spec)

integrated_spectra <- NMR_integration(spectra,regioes_integrais)
integrated_spectra <- cbind("Group"=c(rep("P4C1",8),rep("P6C3",5)),integrated_spectra)


library(ggpubr)

library(tidyr)

boxplot_spec <- gather(integrated_spectra,-Group:-Sample,key="Metabolites",value="Int")

ggboxplot(boxplot_spec,x="Group",y="Int",color="Group")+
  facet_wrap(~Metabolites,scales = "free")+geom_pwc(bracket.nudge.y = -0.5)



