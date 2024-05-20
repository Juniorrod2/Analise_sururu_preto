#instalacao do pacote pepsNMR (So necessario no primero uso, a instalacao permance na maquina)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("PepsNMR")

#instalacao do pacote ropls
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ropls")

#Instacao do pacote NAPRMN
if (!require("devtools", quietly = TRUE))
  install.packages("devtools")

devtools::install_github("Juniorrod2/NAPRMN_install")

#No script de analise estatistica pode ser necessario instalar outros pacotes
#Usar o exemplo abaixo so alterando o nome do pacote a ser instalado
#os 3 listados acima sao excecoes por pertencerem a repositorios diferentes do cran
install.packages("ggplot2")





