## check the current working directory
# and make sure you are placed in Obiettivo2
# setwd("Obiettivo2")

library(tidyverse)
library(terra)
library(CoordinateCleaner)
library(readxl)
library(geodata)
library(patchwork)
library(raster)
library(sf)
library(dplyr)
library(ecospat)
library(rgbif)
library(devtools)
library(ade4)
library(ape)
library(biomod2)
library(maptools)
#library(rgeoboundaries)



### 1.IMPORTAZIONE E FILTRAGGIO OCCORRENZE NATIVE ###
# Importo il database e lo converto in un dataframe
occ <- read.delim("OccorrenzeAmerica2.csv",sep = "\t", header = TRUE, quote = "", stringsAsFactors = FALSE)
head(occ)

0# Creo 4 colonne chiave
occ_native <- occ |>
  dplyr::select(scientificName, decimalLatitude, decimalLongitude, year)

# Pulizia e filtraggio 
occ_native_clean <- occ_native |>
  dplyr::filter(year >= 2014 & year <= 2025) |>
  dplyr::filter(stringr::str_detect(tolower(scientificName), "^procyon lotor")) |>
  dplyr::mutate(scientificName = "Procyon lotor") |>                                   # Uniforma il nome
  dplyr::filter(!is.na(decimalLatitude) & !is.na(decimalLongitude)) |>
  dplyr::mutate(dplyr::across(dplyr::ends_with("tude"), ~as.numeric(.x))) |>  
  dplyr::rename(species = scientificName) |> 
  cc_dupl() |> 
  cc_zero() |> 
  cc_equ() |> 
  cc_val() |> 
  cc_sea() |> 
  cc_cap(buffer = 2000) |> 
  cc_cen(buffer = 2000) |> 
  cc_gbif(buffer = 2000) |> 
  cc_inst(buffer = 2000)

# Visualizzo il risultato
head(occ_native_clean)

# vector
on <- occ_native_clean |> 
  vect(geom=c("decimalLongitude", "decimalLatitude"), crs = "EPSG:4326", keepgeom=FALSE) |> 
  st_as_sf()

plot(on$geometry)

### 2.CARICAMENTO OCCORRENZE TOSCANA NON NATIVE ###
occ_no_native <- read_excel("OccorrenzeToscanaDef.xlsx")
occ_no_native <- subset(occ_no_native, select = -c(4:6))

# vector
onn <- occ_no_native |> 
  mutate(across(starts_with("Coord"), ~as.numeric(.x))) |>                                          #Converti le coordinate in numeri
  dplyr::select(lon = `Coordinata LONG (X)`, lat = `Coordinata LAT (Y)`, dplyr::everything()) |>    #Rinomina e riordina colonne
  vect(geom = c("lon", "lat"), crs = "EPSG:4326") |>                                                #Crea un oggetto spaziale (SpatVector)
  st_as_sf()    #Conversione in oggetto sf

plot(onn$geometry)

#Unione con dati nativi e conversione in data.frame con coordinate esplicite
on_df <- on |> 
  st_drop_geometry() |> 
  mutate(x = st_coordinates(on)[,1],
         y = st_coordinates(on)[,2])

onn_df <- onn |> 
  st_drop_geometry() |> 
  mutate(x = st_coordinates(onn)[,1],
         y = st_coordinates(onn)[,2])

# checkkkkkkk
head(on_df)
head(onn_df)

############### showbiz 
## Armonizza e unisci
soni_df <- bind_rows(
  on_df  |> mutate(ID = NA_real_,          ni = "Native"),
  onn_df |> mutate(species = "Procyon lotor", year = NA_integer_, ni = "Invasive")
)

## A SpatVector
soni <- vect(soni_df, geom = c("x","y"), crs = "EPSG:4326")

## Raster e extract senza colonna ID aggiuntiva
wclim <- rast("world_bio-001.tif")
env   <- extract(wclim, soni, ID = FALSE)
soni  <- cbind(soni, env)

## Tieni solo i record con valori ambientali completi
env_names <- names(wclim)                  
soni <- soni[complete.cases(as.data.frame(soni)[, env_names]), ]

plot(soni)

# crop 
itaExt <- ext(6, 19, 36, 47)
ameExt <- ext(-170, -50, 5, 72)

ocITA <- crop(soni, itaExt)
ocAME <- crop(soni, ameExt)

data(wrld_simpl)
par(mar = c(1, 0, 0, 0))
plot(wrld_simpl, border = "gray80")
points(ocAME, pch = 16, col = 2, cex = 0.3)
points(ocITA, pch = 16, col = 4, cex = 0.3)

dev.off()

# crop the environmental data to the native and invasive geographical ranges
ameEnvR <- crop(wclim, ameExt)
itaEnvR <- crop(wclim, itaExt)

ameEnvR_raster <- raster::stack(ameEnvR)
itaEnvR_raster <- raster::stack(itaEnvR)

ameEnvM <- getValues(ameEnvR_raster)
itaEnvM <- getValues(itaEnvR_raster)

# remove missing values
ameEnvM <- ameEnvM[complete.cases(ameEnvM), ]
ameEnvM

itaEnvM <- itaEnvM[complete.cases(itaEnvM), ]

# produce global environmental background data
globalEnvM <- rbind(ameEnvM, itaEnvM)



### 5. Niche quantification ###
pca.clim <- dudi.pca(globalEnvM, center = TRUE,
                     scale = TRUE, scannf = FALSE, nf = 2)

global.scores <- pca.clim$li

nativeLS.scores <- suprow(pca.clim,
                          data.frame(ocAME)[, colnames(globalEnvM)])$li   

invasiveLS.scores <-
  suprow(pca.clim,
         data.frame(ocITA)[, colnames(globalEnvM)])$li

nativeEnv.scores <- suprow(pca.clim, ameEnvM)$li
invasiveEnv.scores <- suprow(pca.clim, itaEnvM)$li


#Calcolo Occurrence Density Grid per nicchia invasiva e nicchia nativa --> Grafico Niche Overlap
nativeGrid <- ecospat.grid.clim.dyn(global.scores,
                                    nativeEnv.scores,
                                    nativeLS.scores)

invasiveGrid <- ecospat.grid.clim.dyn(global.scores,
                                      invasiveEnv.scores, 
                                      invasiveLS.scores)

ecospat.plot.niche.dyn(nativeGrid, invasiveGrid, quant = 0.1, interest = 2, title = "Niche Overlap", name.axis1 = "PC1", name.axis2 = "PC2")

# Plot Contributo delle variabili (frecce)
ecospat.plot.contrib(contrib=pca.clim$co, eigen=pca.clim$eig)


# Calcolo Expansion, Stability, Unfilling 
ecospat.niche.dyn.index(nativeGrid, invasiveGrid, intersection = 0.1)$dynamic.index.w
# expansion stability unfilling 
# 0.0000000 1.0000000 0.8807437

# Calcolo Niche Overlap (D Schoener)
ecospat.niche.overlap(nativeGrid, invasiveGrid, cor=T)
# $D
# [1] 0.103514

# $I
# [1] 0.2815228

#Niche Equivalency Test
eq.test <- ecospat.niche.equivalency.test(nativeGrid, invasiveGrid, rep = 100, ncores = 2)

#Niche Similarity Test
sim.test <- ecospat.niche.similarity.test(nativeGrid, invasiveGrid, rep = 100, rand.type = 2, ncores = 2)


set.seed(42)
eq.test1  <- ecospat.niche.equivalency.test(nativeGrid, invasiveGrid, rep = 10, ncores = 1)
sim.test1 <- ecospat.niche.similarity.test(nativeGrid, invasiveGrid, rep = 10, rand.type = 2, ncores = 1)


#Plot Equivalency e Similarity Test
par(mfrow=c(1,2))
ecospat.plot.overlap.test(eq.test, "D", "Equivalency") 
ecospat.plot.overlap.test(sim.test, "D", "Similarity")


#Dinamica di nicchia rispetto al biondicatore [1] "Annual_Mean_Temp"
#Provare a farlo anche con gli altri bioindicatori:
# [2] "Mean_Diurnal_Range"    
# [3] "Temp_Seasonality"      
# [4] "Max_Temp_Warmest_Month"
# [5] "Min_Temp_Coldest_Month"
# [6] "Annual_Precipitation"  

# Griglia della nicchia nativa
grid.clim.t.nat <- ecospat.grid.clim.dyn(glob = globalEnvM[,1],
                                         glob1 = data.frame(ameEnvM[,1]),
                                         sp = ocAME$Annual_Mean_Temp, R = 1000, th.sp = 0)
# Griglia della nicchia invasa
grid.clim.t.inv <- ecospat.grid.clim.dyn (glob = globalEnvM[,1], 
                                          glob1 = data.frame(itaEnvM[,1]), 
                                          sp = ocITA$Annual_Mean_Temp, R = 1000, th.sp = 0)

t.dyn <- ecospat.niche.dyn.index (grid.clim.t.nat, grid.clim.t.inv, intersection=0.1)

ecospat.plot.niche.dyn(grid.clim.t.nat, grid.clim.t.inv, quant=0.1, interest=2, title= "Niche Overlap", name.axis1="Annual_Mean_Temp")

# Plot Spostamento centroide lungo il gradiente biocliamtico [1] "Annual_Mean_Temp"
ecospat.shift.centroids(
  sp1  = ocAME$Annual_Mean_Temp,  
  sp2  = ocITA$Annual_Mean_Temp,   
  env1 = ameEnvM[,1],              
  env2 = itaEnvM[,1],      
)
