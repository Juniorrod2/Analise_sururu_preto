#leitura dos dados dos bins a partir de um arquivo binario gerado durante o processamento
#dos espectros
spec <- readRDS("Sururu_spectra")

##leitura dos dados dos espectros processados em sua resolucao original
#a partir de um arquivo binario gerado durante o processamento dos espectros
sururu_spectra <- readRDS("Sururu_preto_espectro_fullRes")

#Criacao do vetor (coluna) de grupos de amostras
bins_sururu <- cbind("Group"=c(rep("P4C1",8),rep("P6C3",5)),bins)

#Inicializacao do pacote para treinamento/ajuste dos dados nos modelos multivariados
library(ropls)

#Inicializacao do pacote NAPRMN contendo funcoes auxiliares para metabolomica
library(NAPRMN)

#Criacao de um modelo de PCA com 4 PCs, armazenando no objeto "pca_sururu"
pca_sururu <- opls(bins_sururu[,c(-1:-2)],predI=4)

#Neste exemplo foi fornecido como segundo argumento o vetor/coluna de grupos (Y), convertendo
#o modelo em uma PLS-DA com 3 componentes. 
pls_sururu <- opls(sururu_data[-1:-2],sururu_data$Group,predI=3)

pls_sururu <- opls(sururu_bin_data[-1:-2],sururu_bin_data$Group,predI=3)

#Adicionando o argumento orthoI incluimos uma componente ortogonal no modelo, que agora
#Sera convertido em um modelo de OPLS-DA
opls_sururu <- opls(bins_sururu[,c(-1:-2)],bins_sururu$Group,predI=3,orthoI=1)

#A partir do objeto armazenando o modelo, podemos usar a função plot() para personalizar
#os graficos de saida e alterar outros paramentros, como quais componentes serão representadas
#(argumento parCompVi), adicionar um vetor de grupos para gerar a coloracao das amostras no
#scores plot (Argumento parAsColFcVn), entre outros. 
plot(pca_sururu,parCompVi=c(1,2),parAsColFcVn=bins_sururu$Group)

#Aqui o modelo de PLS-DA foi plotado apenas com seu scores e loading plot
plot(pls_sururu,typeVc=c("x-score","x-loading"))

#Extracao dos eixos das componentes e dados adcionais (VIP, etc...) a partir do modelo 
#para permitir a plotagem dos dados com uma biblioteca grafica externa
pls_sururu_data <- extract_ropls_data(pls_sururu)
pls_sururu_data$Loadings$bins <- rownames(pls_sururu_data$Loadings)

#Inicializacao do pextract_ropls_data()#Inicializacao do pacote ggplot2, um pacote grafico que pode ser utilzado para produzir 
#diversos tipos de graficos a partir de uma planilha/dataframe
library(ggplot2)

#Plotagem do scores plot com o pacote ggplot2 a partir dos dados extraidos do modelo
#de PLS-DA
scores <- ggplot(pls_sururu_data$Scores,aes(p1,p2,color=Group,fill = Group))+
  geom_point(size=3)+stat_ellipse(geom="polygon",alpha=0.4)+theme_bw()+
  labs(x="Predictive component 1",y="Predictive component 2",title = "PLS-DA")

#Plotagem do loading plot com o pacote ggplot2 a partir dos dados extraidos do modelo
#de PLS-DA
loadings <- ggplot(pls_sururu_data$Loadings[1:24,],aes(p1,p2,color=bins,label=bins)) +
  geom_text()#+scale_color_gradient(high = "red",low="darkgray")

#para visualizar os graficos basta "chamar" os objetos que armazenaram os plots
scores
loadings

#Se preferir o graficos tambem podem ser gerados sem serem armazenados em um objeto e
#serem plotados imediatamente
ggplot(pls_sururu_data$Scores,aes(p1,p2,color=Group,fill = Group))+
  geom_point(size=3)+stat_ellipse(geom="polygon",alpha=0.4)+theme_bw()+
  labs(x="Predictive component 1",y="Predictive component 2",title = "PLS-DA")

ggplot(pls_sururu_data$Loadings,aes(p1,p2,color=Vip,label=bins)) +
  geom_text()+scale_color_gradient(high = "red",low="darkgray")


#Inicializacao do pacote plotly que permite converter os graficos estaticos do ggplot2
#em figuras interativas que permite zoom e alteracoes dinamicas no grafico
library(plotly)

#esta funcao do pacote plotly faz a conversao do ggplot em um objeto interativo
ggplotly(scores)
ggplotly(loadings)

#inicializado do pacote para organizacao de diferentes plots em um unico painel
library(patchwork)

#Com o pacote patchwork carregado e possivel "adicionar" os graficos uns aos outros 
#e eles serao unidos em um unico painel. Os operadores / e () tambem produzem outros formatos
scores+loadings
scores/loadings
scores+(loadings/loadings)


#Inicalizacao do pacote para leitura de planilhas no formato xlsx 
library(readxl)
#Leitura da planilha contendo as regioes de integracao de cada pico em um objeto do tipo
#datafrem
regioes_integrais <- read_excel("regioes_integrais.xlsx")



# a funcao abaixo recebe um objeto contendo a matriz de espectros processados e uma planilha
#contedo as regioes de corte de cada sinal, permindo integrar regioes especificas do espectro
integrated_spectra <- NMR_integration(sururu_spectra,regioes_integrais)
integrated_spectra <- cbind("Group"=c(rep("P4C1",8),rep("P6C3",5)),integrated_spectra)

#Normalizacao e escalonamento dos dados 

normalized_integrals <- DataNormalization(integrated_spectra[-1:-2],rowNorm = "MedianNorm",
                                          transNorm ="LogNorm",scaleNorm = "ParetoNorm")

integrated_spectra <- cbind(integrated_spectra[1:2],normalized_integrals)

# Extensao do ggplot2 que permite gerar graficos contendo testes estatisticos de hipotese
library(ggpubr)

#Pacote utilizado para reorganizacao de dados, aqui a funcao gather() que reoganiza o formato
#dos dados para permitir gerar o painel com os boxplots vem deste pacote
library(tidyr)

boxplot_spec <- gather(integrated_spectra,-Group:-Sample,key="Metabolites",value="Int")

#funcao do pacote ggpubr que produz boxplots dos dados na planilha/dataframe
#boxplot_spec e com a adicao do layer/camada geom_pwc() calcula e insere na figura
#os resultados dos testes de hipoteses, ver a documentacao para as variacoes de testes
#de hipotese permitidos
ggboxplot(boxplot_spec,x="Group",y="Int",color="Group")+
  facet_wrap(~Metabolites,scales = "free")+geom_pwc(bracket.nudge.y = -0.15,hide.ns = T,
                                                    label = "p.signif")
#Ou em outro estilo
ggboxplot(boxplot_spec,x="Group",y="Int",color="Group")+
  facet_wrap(~Metabolites,scales = "free")+stat_compare_means(vjust = 1)


