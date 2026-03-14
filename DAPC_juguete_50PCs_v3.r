################################################################################
# Samantha M. Pacheco Gómez
# 13 marzo 2026
# DAPC_juguete_50PCs_v3.r
# formato script para división por clusters k= 5
# data teosinte
# versión cluster 5 clusters, 50 PCs, data worldclim recortada
################################################################################

# cargar librerías

library(dplyr)
library(terra)
library(geodata)
library(adegenet)
library(SNPRelate)
library(ggplot2)

# --- 1. CONFIGURACIÓN DE RUTAS ---
# Definimos las rutas como variables para que sea fácil cambiarlas después
path_teosinte = "/mnt/data/sur/users/spacheco/data/teosinte/"
path_resultados = "/mnt/data/sur/users/spacheco/resultados/"
path_clima = "/mnt/data/sur/users/spacheco/data/worldclim/mesoamerica_30s.tif"

cat("Rutas configuradas\n")
# --- CARGA DE METADATA Y CLIMA ---
Data = read.csv2(paste0(path_teosinte, "data_teosinte.csv"), 
                 header = TRUE, sep=",", as.is = FALSE)

# Convertimos coordenadas a numérico para que terra no de problemas
Data$Latitude  = as.numeric(as.character(Data$Latitude))
Data$Longitude = as.numeric(as.character(Data$Longitude))

cat("Metadata cargada\n")

# Cargamos el archivo único (automáticamente detecta las 19 capas)

clima_raw = rast(path_clima)

nombres_chidos = c("Temp_Media_Anual", "Rango_Diurno", "Isotermalidad", 
                   "Estacionalidad_Temp", "Max_Temp_Mes_Calido", "Min_Temp_Mes_Frio", 
                   "Rango_Anual_Temp", "Temp_Cuartil_Humedo", "Temp_Cuartil_Seco", 
                   "Temp_Cuartil_Calido", "Temp_Cuartil_Frio", "Prec_Anual", 
                   "Prec_Mes_Humedo", "Prec_Mes_Seco", "Estacionalidad_Prec", 
                   "Prec_Cuartil_Humedo", "Prec_Cuartil_Seco", "Prec_Cuartil_Calido", 
                   "Prec_Cuartil_Frio")
names(clima_raw) = nombres_chidos

cat("Data de clima mesoamerica cargada\n")

# --- PROCESAMIENTO GEOGRÁFICO ---
# Primero calculamos los promedios por población

geo_pop = Data %>%
  group_by(POB_CODE) %>%
  summarise(lat = mean(Latitude, na.rm=TRUE),
            lon = mean(Longitude, na.rm=TRUE),
            .groups="drop") %>%
  filter(!is.na(lat) & !is.na(lon))

# Ahora extraemos los valores (esto será instantáneo ahora)
valores_climaticos = extract(clima_raw, geo_pop[, c("lon", "lat")])

# Unimos todo
geo_pop_clim = cbind(geo_pop, valores_climaticos)
Data_Final = left_join(Data, geo_pop_clim %>% select(-lat, -lon, -ID), by = "POB_CODE")

cat("procesamiento geográfico terminado\n")

# ---  CONVERSIÓN Y LECTURA GENÉTICA ---
# Solo crea el GDS si no existe todavía en la carpeta teosinte
archivo_gds = paste0(path_teosinte, "teosinte.gds")

if (!file.exists(archivo_gds)) {
  snpgdsBED2GDS(bed.fn = paste0(path_teosinte, "T3604_33929_all.bed"), 
                fam.fn = paste0(path_teosinte, "T3604_33929_all.fam"), 
                bim.fn = paste0(path_teosinte, "T3604_33929_all.bim"), 
                out.gdsfn = archivo_gds)
}

genofile = snpgdsOpen(archivo_gds)
geno_mat = snpgdsGetGeno(genofile)
samp_ids = read.gdsn(index.gdsn(genofile, "sample.id"))
snp_ids  = read.gdsn(index.gdsn(genofile, "snp.id"))
snpgdsClose(genofile)

gl_teosinte = new("genlight", geno_mat)
indNames(gl_teosinte) = samp_ids
locNames(gl_teosinte) = snp_ids

cat("conversión y lectura genética terminada\n")

# ---  DAPC DIAGNÓSTICO (PRUEBA DE 50 PCs) ---
set.seed(147) # para que sea replicable

# Buscamos 5 clusters usando solo 50 PCs para que sea rápido
grupos_test = find.clusters(gl_teosinte, n.pca = 50, n.clust = 5, choose.n.clust = FALSE)

# Corremos el DAPC también con 50 PCs
dapc_test = dapc(gl_teosinte, pop = grupos_test$grp, n.pca = 50, n.da = 4)

cat("DAPC terminado\n")

# ---  GUARDAR RESULTADOS ---
# pdf con nombre significativo segun yo jiji
pdf(paste0(path_resultados, "DAPC_juguete_50PCs_v3.pdf"), width = 10, height = 8)

# Gráfica del nicho térmico por taxón
print(ggplot(Data_Final, aes(x=Taxon, y=Temp_Media_Anual, fill=Taxon)) +
        geom_boxplot() + theme_minimal() +
        labs(title="Diferencias de Temperatura Media (30s res)", y="Temp (°C)") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)))

# Gráfica del DAPC (La constelación de puntos)
scatter(dapc_test, col = spectral(5), scree.da = FALSE, bg = "white", pch = 20, 
        main = "DAPC Diagnóstico (Mesoamérica 30s - 50 PCs - k=5")

dev.off()

cat("pdf guardado\n")

# Guardamos la tabla final que une Genética + Clima de 30s
pertenencia = as.data.frame(dapc_test$posterior)
pertenencia$Sample_name = rownames(pertenencia)
Data_Master = left_join(Data_Final, pertenencia, by = "Sample_name")

# Guardamos el .RData también en la carpeta de resultados
save(Data_Master, dapc_test, file = paste0(path_resultados, "DAPC_juguete_50PCs_v3.RData"))

cat("¡Proceso terminado con éxito!\n")


