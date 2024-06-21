# ---- 00. Pacotes para instalação.----

# COMECE AQUI ~ Comando para instalação:  ---- 
# install.packages(c("ggspatial", "raster", "sdm", "tidyverse", "spThin","geobr") 

# Algumas dependências não puderam ser resolvidas pelo R, no caso do meu 
# sistema atual, Archlinux, foi preciso instalá-las através do repositório
# AUR, com os seguintes comandos: 

# yay -S r-classint
# yay -S udunits

# Em outras distribuições Linux ou no Windows o método de instalação 
# vai ser diferente.



# ---- 01. Bibliotecas Requeridas ----

library(ggspatial)
library(raster)
library(sdm)
library(tidyverse)
library(geobr)# Necessário para baixar shapefiles
#library(spThin) # Para ajustar e selecionar coordenadas espaçadas em intervalos regulares
library(scales)



# Pacotes complementares do sdm precisam ser instalados com o comando:
#sdm::installAll() # RODAR 1 VEZ após instalar o sdm


# ---- 03. Shapefiles com o pacote geobr ----

# Baixar os shapes separadamente
shape <- read_state(
  code_state = "all", 
  year = 2020, 
  simplified = FALSE
)

shape = shape %>% filter(name_region == "Nordeste")

# Salvar o shape (Opicional)
sf::st_write(shape, "shapefiles/nordeste.shp", append = FALSE)
sf::st_write(shape, "shapefiles/nordeste.geojson", append = FALSE)
#sf::st_write(br, "shapefiles/br.geojson")

# Visualizar o Shape
plot(shape$geom) 

ggplot() + 
  geom_sf(data = shape_country, fill = "lightgray", color = "black") +
  theme_minimal() +
  labs(title = "Title", 
       captio = "Fonte: geobr")

# 04. Dados de Ocorrência da espécie ----

# GBIF
b_apod_ocorrencia_gbif <- dismo::gbif( genus = "Borreria", 
                                            species = "apodiensis", 
                                            geo = TRUE, 
                                            removeZeros = TRUE,
                                            download = TRUE
)

gbifcsv <- read.csv("dados/ocorrencia/b_apodiensis_ocorrencia_gbif.csv")
gbifcsv %>% select(species, lat, lon) 
write.csv(gbifcsv, "dados/pontosParaGeocat.csv")
write.csv(b_apod_ocorrencia_gbif, 
      file = "dados/ocorrencia/b_apodiensis_ocorrencia_gbif.csv") # Salvar daos do GBIF em um arquivo CSV

b_apod_ocorr_gbif <- read.csv(
  file = "dados/ocorrencia/b_apodiensis_ocorrencia_gbif.csv")

# Dados de coleções e artigos
b_apod_ocorr_lit <- read.csv(
  file = "dados/ocorrencia/b_apodiensis_literatura.csv")

# Dados do SpeciesLink
b_apod_ocorr_specieslink <- read.delim(
  file = "dados/ocorrencia/speciesLink-20240618192743-0010992.txt")

# Dados de coletas
b_apod_ocorr_coletas <- read.csv(
  file = "dados/ocorrencia/b_apodiensis_ocorrencia_coletas.csv")

# Tratamento dos dados
# GBIF
b_apod_tt_gbif <- b_apod_ocorr_gbif %>% 
  select(species, lat, lon) 

# Literatura
b_apod_tt_lit <- b_apod_ocorr_lit %>% 
  select(Latitude, Longitude) %>% 
  mutate(species = "Borreria apodiensis", .before = Latitude)

names(b_apod_tt_lit) <- c("species", "lat", "lon") # Padroniza os nomes das colunas.

# Species Link
b_apod_tt_specieslink <- b_apod_ocorr_specieslink %>% 
  select(scientificname, latitude, longitude) %>% 
  filter(latitude != 0 & longitude !=0) # Remove entradas do SP que não tinham coordenadas.

names(b_apod_tt_specieslink) <- c("species", "lat", "lon")

# Coletas
b_apod_tt_coletas <- b_apod_ocorr_coletas %>% 
  select(nome_cientifico, latitude, longitude) 

names(b_apod_tt_coletas) <- c("species", "lat", "lon")

# Compilado dos dados

b_apod <- bind_rows(
  b_apod_tt_gbif,
  b_apod_tt_lit,
  b_apod_tt_specieslink,
  b_apod_tt_coletas
) %>% 
  drop_na() %>%
  distinct() # combina todas as tabelas e exclui os NA e duplicados.

# b_apod_thin <- thin(
#   loc.data = b_apod,
#   lat.col = "lat",
#   long.col = "lon",
#   spec.col = "species",
#   thin.par = 0.5, # Especialização em km
#   reps = 1,
#   locs.thinned.list.return = TRUE,
#   write.files = FALSE,
#   write.log.file = FALSE
# )
# 
# b_apod_thin <- b_apod_thin[[1]]
# b_apod_thin <- b_apod_thin %>% select(Latitude, Longitude)
# names(b_apod_thin) <- c("lat", "lon")

# Excluíndo pontos discrepantes
b_apod_ok <- b_apod %>% filter(lon <= -37 & lon > -39)

write.csv(b_apod_ok, file = "dados/ocorrencia/b_apodiensis_distinct_coherent.csv", row.names = FALSE)
b_apod_ok <- read.csv("dados/ocorrencia/b_apodiensis_distinct_coherent.csv") # IMPORTANTE Coordenadas Distintas ----

pontos <- b_apod_ok  # p/ ggplot

b_apod_pts <- b_apod_ok %>% select(lon, lat)
b_apod_pts <- read.csv("dados/ocorrencia/b_apodiensis_pts.csv")# CHECKPOINT ~ b_apo_pts ---- 

#write.csv(b_apod_pts, "dados/ocorrencia/b_apodiensis_pts.csv", row.names = FALSE) 


#Mapa dos pontos para conferência com o ggplot
ggplot() +
  geom_sf(data = shape, fill = "lightgray") +
  geom_point(data = b_apod_pts, aes(x = lon, y = lat)) +
  coord_sf() +
  theme_minimal()

#ggsave(filename = "projecoes/plot_1_b_apodiensis.svg", plot = plot_1)

# 05. Dados de variáveis bioclimáticas ----

# Download das variáveis Bioclimáticas do WorldClim 
geodata::worldclim_global( # VAR BIO PASSO 0 ----
  var = "bio",
  res = 5, # resolução em graus/minutos
  path = "dados/bioclim/",
)

bioclim_var <- raster::stack( # CHECKPOINT ~ VAR BIO PASSO 1 ----
  list.files(
    path = "dados/bioclim/climate/wc2.1_5m/",
    pattern = ".tif",
    full.names = TRUE
  )
) 

bioclim_var <- subset( # CHECKPOINT ~ VAR BIO PASSO 2 ----
  bioclim_var, 
  c(
  "wc2.1_5m_bio_1",
  "wc2.1_5m_bio_2",
  "wc2.1_5m_bio_3",
  "wc2.1_5m_bio_4",
  "wc2.1_5m_bio_5",
  "wc2.1_5m_bio_6",
  "wc2.1_5m_bio_7",
  "wc2.1_5m_bio_8",
  "wc2.1_5m_bio_9",
  "wc2.1_5m_bio_10",
  "wc2.1_5m_bio_11",
  "wc2.1_5m_bio_12",
  "wc2.1_5m_bio_13",
  "wc2.1_5m_bio_14",
  "wc2.1_5m_bio_15",
  "wc2.1_5m_bio_16",
  "wc2.1_5m_bio_17",
  "wc2.1_5m_bio_18",
  "wc2.1_5m_bio_19"
  )
)

#Converter CRS para SIRGAS 2000
bioclim_var <- raster::projectRaster( # CHECKPOINT ~ VAR BIO PASSO 3 ----
  bioclim_var,
  crs = raster::crs(shape)
)

plot(bioclim_var[[1]]) #Teste

bio
writeRaster( # IMPORTANTE ~ Salvar rasterbrick ----
  bioclim_var, 
  filename = "dados/bioclim/climate/bio_var_raster", 
  overwrite=TRUE
)

# IMPORTANTE ~ Carregar direto o rasterbrick ----
bioclim_var <- raster::brick("dados/bioclim/climate/bio_var_raster") 

# Recorte Shape
bioclim_var <- crop(bioclim_var, shape) # CHECKPOINT ~ Recorte do SF das Variáveis Bio. ----
bioclim_var <- mask(bioclim_var, shape) # CHECKPOINT ~ Máscara do mesmo SF
plot(bioclim_var[[15]])

writeRaster( # CHECKPOINT ~ Salvar RasterBrick das variáveis ----
  bioclim_var, 
  filename = "dados/bioclim/climate/bio_NE", 
  overwrite=TRUE
) 




# 06. Tratamento das Variáveis Bioclimáticas ----

# Relaciona as varíaveis em função das coordenadas
b_apo <- raster::extract(bioclim_var, b_apod_pts) # IMPORTANTE Var(coord) ----


# Teste Fator de Inflação da Variância (VIF): análise de multicolinearidade
b_apod_vif <- usdm::vifcor(b_apo, th = 0.9) # CHECKPOINT ~ VIFSTEP/VIFCOR ----


b_apod_vif # IMPORTANTE Relatório do Teste VIF ----



# CHECKPOINT ~ Exclusão das variáveis colineares ----
b_apo_bioclim <- usdm::exclude(bioclim_var, b_apod_vif) 


writeRaster(# CHECKPOINT ~ Salvar RasterStack ----
  b_apo_bioclim, 
  "dados/bioclim/climate/b_apo_bioclim_variaveis_pos_vif",
  overwrite=TRUE
) 
# b_apo_bioclim <- raster::brick("dados/bioclim/climate/b_apo_bioclim_variaveis_pos_vif.tif") # SALVAR Tif




#Plot_2
plot(b_apo_bioclim) 
points(b_apod_pts, col="red")



# 07 Modelagem ----

#CHECKPOINT Corrigir nomes das VB OPICIONAL ----
names(b_apo_bioclim)
names(b_apo_bioclim) <- c("BIO10", "BIO11", "BIO14", "BIO19") #, "BIO13", "BIO14", "BIO18", "BIO19") 


# RasterBrick -> RasterStack
b_apo_bioclim <- raster::stack(b_apo_bioclim)


# Adiciona coluna de número 1.
b_apod_pts <- b_apod_pts %>% # CHECKPOINT Coluna de Presença ----
  mutate(borreria_apodiensis = 1) 

# Cria um SpatialPointsDataFrame com os dados de presença
coordinates(b_apod_pts) <- c("lon", "lat") # CHECKPOINT Georeferenciamento ----


model_data <- sdmData( # CHECKPOINT Model Data ----
  formula = borreria_apodiensis ~.,
  train = b_apod_pts,
  predictors = b_apo_bioclim,
  bg = list(
    n = 100,
    method = "eRandom",
    remove = TRUE
  )
)
model_data



model <- sdm( # CHECKPOINT Models ----
  formula = borreria_apodiensis ~.,
  data = model_data,
  methods = c("maxent","glm"),
  replication = "sub",
  n = 10,
  test.percent = 20,
  parallelSetting = list(
    ncore = 2,
    "parallel"
  )
) 
model
roc(model)
getVarImp(model)


dir <- "modelos/Maxent_GLM_eRandom/"
dir.create(
  dir,
  recursive = TRUE
)

sdm::write.sdm( # SALVAR Dados de Modelo
  model_data, 
  paste(dir, "model_data", sep = "")
)
model_data <- sdm::read.sdm(paste(dir, "model_data.sdd", sep = ""))


sdm::write.sdm( # SALVAR Modelo
  model, 
  paste(dir, "model", sep ="")
)
model <- sdm::read.sdm(paste(dir, "model.sdm", sep = ""))


writeRaster(
  b_apo_bioclim,
  paste(dir, "b_apo_bioclim", sep = "")
)
b_apo_bioclim <- raster::stack(paste(dir, "b_apo_bioclim.gri", sep = ""))


# 08 Pojeções ----
# a. Presente 
proj_b_apo <- predict(
  model,
  newdata = b_apo_bioclim,
  filename = "modelos/Maxent_GLM_eRandom/b_apo_maxent_GLM.grd",
  overwrite = TRUE,
  #nc = 2
)

raster::brick(b_apo_bioclim)
 raster::stack(b_apo_bioclim)
b_apo_ens <- ensemble(
  model,
  newdata = b_apo_bioclim,
  setting = list( 
  method = "weighted",
  stat = "TSS",
  opt = 2 #Max specificity
  )
)

# b. Futuro 




# 09. Mapas no ggplot ----

### Criar mapas finais -----



layer1 <- raster(proj_b_apo[[1]])

adeq_ambiental <- as.data.frame(
  raster(proj_b_apo[[1]]),
  xy = TRUE
)

# Paleta de cores LASEP :)
pal <- colorRampPalette(c("#fdf8b8", "#557755", "#006153"), bias = 1)
pal(5)

plot_1 <- ggplot() +
  geom_raster(data = adeq_ambiental, aes(x, y, fill = layer)) +
  geom_sf(data = shape, color = "#464646", fill = NA) +
  geom_point(
    data = pontos,
    aes(x = lon, y = lat), size = .4, color = "red"
  ) +
  scale_fill_gradientn(
    colours = pal(5),
    na.value = NA,
    limits = c(0.00, 1.00)
  ) + # parâmetros da escala de cor
  coord_sf() +
  annotation_north_arrow(
    location = "br",
    which_north = "true",
    style = north_arrow_fancy_orienteering,
    pad_x = unit(0.25, "in"), pad_y = unit(0.30, "in")
  ) +
  annotation_scale(location = "br", width_hint = 0.5) +
  labs(
    x = "Longitude", # texto do eixo x
    y = "Latitude", # texto do eixo y
    fill = "Adequabilidade"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_line(color = "lightgrey"),
    panel.grid.minor = element_line(color = "lightgrey"),
    legend.title = element_text(face = "bold", size = 12, vjust = 0.9)
  )



plot_1
names(adeq_ambiental) <- c("x", "y", "layer")
