# 00. Pacotes para instalação.----
# Comando para instalação:

# install.packages(c("ggspatial", "raster", "sdm", "tidyverse", "spThin"), 
# clean = TRUE, dependencies = TRUE)

# Algumas dependências não puderam ser resolvidas pelo R, no caso do meu 
# sistema atual, Archlinux, foi preciso instalá-las através do repositório
# AUR, com os seguintes comandos: 

# yay -S r-classint
# yay -S udunits

# 01. Bibliotecas Requeridas ----
library(ggspatial)
library(raster)
library(sdm)
library(tidyverse)
library(spThin)


# 02. Paleta de cores LASEP :) ----
pal <- c("#B4A54B",
         "#D7FCA7",
         "#A3FF84",
         "#56D54A",
         "#00A537")

# 03. Dados de Ocorrência da espécie ----
# 04. Dados de variáveis bioclimáticas ----
# a. Presente ----
# b. Futuro ----
# 05. Modelagem
# 06. Testes estatísticos do(s) modelo(s)
# 07. Plotagem de mapas com Ensamble

