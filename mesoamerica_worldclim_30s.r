################################################################################
# Samantha M. Pacheco Gómez
# 13 marzo 2026
# mesoamerica_30s.r
# formato script para recortar el dataset de worldclim de 30 segundos
# a solo las coordenadas correspondientes a Mesoamérica (México a Costa Rica)
# data worldclim
# este código se corrió de manera interactiva en el cluster fenix
################################################################################

# module load r
# module load proj
# module load gdal
# module load geos
# 
# R       #inicializa r
# q()     #cierra r en el cluster

library(terra)

# 1. Cargar los 19 archivos globales
path_clima = "/mnt/data/sur/users/spacheco/data/worldclim/climdata_30s/"
listado_tifs = list.files(path_clima, pattern = ".tif$", full.names = TRUE)
listado_tifs = listado_tifs[order(as.numeric(gsub("\\D", "", basename(listado_tifs))))]
clima_global = rast(listado_tifs)

# 2. Definir el recorte (Oeste, Este, Sur, Norte)
ext_meso = ext(-118, -82, 8, 33)

# 3. Recortar
clima_meso = crop(clima_global, ext_meso)

# 4. GUARDAR como un solo archivo de 19 capas (Esto es lo que te ahorrará horas)
writeRaster(clima_meso, 
            filename = "/mnt/data/sur/users/spacheco/data/worldclim/mesoamerica_30s.tif", 
            overwrite = TRUE)

 # esto fue lo que sucedió una vez me asignaron un nodo de cómputo:

# > library(terra)
# terra 1.8.60
# > # 1. Cargar los 19 archivos globales
#   path_clima = "/mnt/data/sur/users/spacheco/data/worldclim/climdata_30s/"
# listado_tifs = list.files(path_clima, pattern = ".tif$", full.names = TRUE)
# listado_tifs = listado_tifs[order(as.numeric(gsub("\\D", "", basename(listado_tifs))))]
# clima_global = rast(listado_tifs)
# > # 2. Definir el recorte (Oeste, Este, Sur, Norte)
#   ext_meso = ext(-118, -82, 8, 33)
# > # 3. Recortar
#   clima_meso = crop(clima_global, ext_meso)
# > # 4. GUARDAR como un solo archivo de 19 capas (Esto es lo que te ahorrará horas)
#   writeRaster(clima_meso, 
#               filename = "/mnt/data/sur/users/spacheco/data/worldclim/mesoamerica_30s.tif", 
#               overwrite = TRUE)
# > clima_meso
# class       : SpatRaster 
# size        : 3000, 4320, 19  (nrow, ncol, nlyr)
# resolution  : 0.008333333, 0.008333333  (x, y)
# extent      : -118, -82, 8, 33  (xmin, xmax, ymin, ymax)
# coord. ref. : lon/lat WGS 84 (EPSG:4326) 
# source(s)   : memory
# names       : wc2.1~bio_1, wc2.1~bio_2, wc2.1~bio_3, wc2.1~bio_4, wc2.1~bio_5, wc2.1~bio_6, ... 
# min values  :     2.05000,    3.741667,    25.90498,    33.65096,         9.4,        -8.4, ... 
# max values  :    28.77917,   21.358334,    88.14815,   853.50201,        43.9,        23.5, ... 
# > 
  
