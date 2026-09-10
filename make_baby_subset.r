# ============================================================
# subset bebé de pruebas
# 9 sept 2026

# Armar un subset chico (3 poblaciones, 1,000 SNPs) para probar
# que el pipeline de estructura corre bien y ver tiempos, antes
# de mandar la corrida completa al cluster.
# ============================================================

# ============================================================
# 1. RUTAS Y CONFIGURACIÓN
# ============================================================

RUTA_BFILE <- "/home/sam/Documents/sur_ecoevo_lab/data/teosinte/archivos/T3604_33929_all"   # sin extensión
RUTA_METADATA <- "/home/sam/Documents/sur_ecoevo_lab/data/teosinte/archivos/data_teosinte.csv"  # para asignar accession a cada individuo
RUTA_SALIDA_SUBSET <- "/home/sam/Documents/sur_ecoevo_lab/data/teosinte/archivos/T3604_baby_subset"

n_poblaciones <- 3   
n_snps <- 1000       
semilla <- 44         


# ============================================================
# 2. Elegir 3 Accessions con un tamaño de muestra razonable
# ============================================================

meta <- read.csv(RUTA_METADATA, stringsAsFactors = FALSE)
conteo_por_pop <- table(meta$Accession)

# solo considerar poblaciones con al menos 8 individuos, para que la prueba se parezca un poco más a una corrida real
candidatas <- names(conteo_por_pop[conteo_por_pop >= 8])

set.seed(semilla)
pops_elegidas <- sort(sample(candidatas, n_poblaciones))
cat("Poblaciones elegidas para la prueba:", paste(pops_elegidas, collapse = ", "), "\n")
cat("Individuos por población:\n")
print(conteo_por_pop[pops_elegidas])

individuos_elegidos <- meta[meta$Accession %in% pops_elegidas, ]


# ============================================================
# 3. Escribir el archivo --keep para plink2 (necesita FID + IID exactos 
# tal como aparecen en el .fam, por eso se sacan del .fam real, no se reconstruyen a mano)
# ============================================================

fam <- read.table(paste0(RUTA_BFILE, ".fam"), stringsAsFactors = FALSE,
                  col.names = c("FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE"))

fam_filtrado <- fam[fam$IID %in% individuos_elegidos$Sample_name, c("FID", "IID")]
cat("Individuos encontrados en el .fam para esas 3 poblaciones:", nrow(fam_filtrado), "\n")

write.table(fam_filtrado, "keep_individuos_prueba.txt",
            row.names = FALSE, col.names = FALSE, quote = FALSE)


# ============================================================
# 4. Elegir 1,000 SNPs al azar (de los ~33,929) y escribir la lista para --extract
# ============================================================

bim <- read.table(paste0(RUTA_BFILE, ".bim"), stringsAsFactors = FALSE,
                  col.names = c("CHR", "SNP_ID", "CM", "BP", "A1", "A2"))

set.seed(semilla)
snps_elegidos <- sample(bim$SNP_ID, min(n_snps, nrow(bim)))
writeLines(snps_elegidos, "extract_snps_prueba.txt")
cat("SNPs elegidos para la prueba:", length(snps_elegidos), "de", nrow(bim), "\n")


# ============================================================
# 5. Correr plink2 para generar el subset
# ============================================================

resultado <- system2("plink2", c(
  "--bfile", RUTA_BFILE,
  "--keep", "keep_individuos_prueba.txt",
  "--extract", "extract_snps_prueba.txt",
  "--make-bed",
  "--out", RUTA_SALIDA_SUBSET
))

if (resultado != 0) {
  stop("plink2 no corrió bien armando el subset -- revisar mensaje de arriba.")
}

cat("\nSubset de prueba listo en:", RUTA_SALIDA_SUBSET, ".bed/.bim/.fam\n")


