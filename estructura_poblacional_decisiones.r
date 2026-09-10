# ============================================================
# PCA + DAPC ANTES de cualquier filtro (línea base para comparar contra los
# estructura_poblacional_decisiones.r
# decisiones porque apenas decidiré si vale la pena usar los 3300 componentes con los
# que se quedan ellos o si hay mejores parámteros
# ============================================================

library(adegenet)

reportar_tiempo <- function(etiqueta, t_referencia) {
  transcurrido <- as.numeric(difftime(Sys.time(), t_referencia, units = "mins"))
  cat(sprintf("[tiempo] %s: %.2f min\n", etiqueta, transcurrido))
  Sys.time()
}
t0 <- Sys.time()
t_inicio_total <- t0

# ============================================================
# 1. RUTAS Y CONFIGURACIÓN
# ============================================================
# cluster

RUTA_BFILE <- "/mnt/data/sur/users/spacheco/data/teosinte/T3604_33929_all"   # sin extensión
RUTA_METADATA <- "/mnt/data/sur/users/spacheco/data/teosinte/data_teosinte.csv"  # para asignar accession a cada individuo
RUTA_CARPETA_SALIDA <- "/mnt/data/sur/users/spacheco/results/sep_2026/teo/9_sep/"

# Rango de K a probar en find.clusters (el paper probó 1 a 40)
K_MIN <- 1
K_MAX <- 40

# Número de PCs "grande" para replicar el enfoque del paper
# (ahí usaron 3,300 -- ver la nota de la sección 4 sobre por qué
# ese número específico no se explica en el texto del paper)
N_PCS_GRANDE <- 3300  # <-- ajustable

# Configuración de xvalDapc (validación cruzada para elegir PCs)
XVAL_N_PCA_MAX <- 200   # <-- rango de PCs a explorar; ajustar si hace falta
XVAL_N_REP <- 30        # <-- default de adegenet; bajar si tarda demasiado


# ============================================================
# 2. Convertir .bed/.bim/.fam a .raw con plink2, y arreglar el
#    delimitador (plink2 usa tabs; adegenet::read.PLINK necesita
#    espacios -- por eso NO se usa "plink1.9 --recodeA" aquí, para
#    no depender de un programa aparte que quizás no tengas
#    instalado; plink2 ya lo tienes)
# ============================================================

dir.create(RUTA_CARPETA_SALIDA, recursive = TRUE, showWarnings = FALSE)
ruta_raw_original <- paste0(RUTA_BFILE, ".raw")
ruta_raw_convertido <- paste0(RUTA_BFILE, "_convertido.raw")

if (!file.exists(ruta_raw_convertido)) {
  message("Generando .raw con plink2 --export A...")
  resultado <- system2("plink2", c("--bfile", RUTA_BFILE, "--export", "A", "--out", RUTA_BFILE))
  
  if (resultado != 0 || !file.exists(ruta_raw_original)) {
    stop("plink2 no corrió correctamente (¿está instalado? probá 'which plink2' en la ",
         "terminal). No se puede continuar sin el .raw generado.")
  }
  
  message("Convirtiendo delimitador (tabs -> espacios) para que adegenet lo lea...")
  lineas <- readLines(ruta_raw_original)
  lineas <- gsub("\t", " ", lineas)
  writeLines(lineas, ruta_raw_convertido)
} else {
  message(ruta_raw_convertido, " ya existe, no se regenera. Bórralo a mano si cambiaste el .bed/.bim/.fam.")
}
t0 <- reportar_tiempo("conversión a .raw", t0)


# ============================================================
# 3. Leer como genlight y asignar población (Accession)
# ============================================================

gl <- read.PLINK(ruta_raw_convertido, quiet = TRUE)
cat("Individuos leídos:", nInd(gl), "| Loci:", nLoc(gl), "\n")

# Asignar Accession como población, cruzando por Sample_name.
# OJO: esto asume que los IID del .fam coinciden con Sample_name
# de data_teosinte.csv -- confírmalo si algo no cuadra (debería
# ser así porque ese archivo fue la fuente para armar el pipeline
# de genotipos desde el principio).
meta <- read.csv(RUTA_METADATA, stringsAsFactors = FALSE)
orden_gl <- data.frame(IID = indNames(gl))
cruce <- merge(orden_gl, meta[, c("Sample_name", "Accession")],
               by.x = "IID", by.y = "Sample_name", all.x = TRUE, sort = FALSE)
# el merge no garantiza el orden -- reordenar según indNames(gl)
cruce <- cruce[match(indNames(gl), cruce$IID), ]

faltantes <- sum(is.na(cruce$Accession))
if (faltantes > 0) {
  message("[AVISO] ", faltantes, " individuo(s) del genlight no encontraron ",
          "Accession en ", RUTA_METADATA, " -- revisar antes de continuar.")
}

pop(gl) <- cruce$Accession
cat("Poblaciones (Accession) distintas asignadas:", length(unique(pop(gl))), "\n")
t0 <- reportar_tiempo("leer genlight + asignar Accession", t0)


# ============================================================
# 4. PCA (independiente del DAPC, para la gráfica de PCA sola)
# ============================================================

n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "1"))
cat("Usando", n_cores, "núcleo(s) para glPca (detectado de SLURM_CPUS_PER_TASK)\n")
pca <- glPca(gl, nf = 10, parallel = (n_cores > 1), n.cores = n_cores)

varianza_pca <- pca$eig / sum(pca$eig) * 100
write.csv(
  data.frame(eje = paste0("PC", seq_along(varianza_pca)),
             varianza_pct = varianza_pca)[1:10, ],
  paste0(RUTA_CARPETA_SALIDA, "pca_varianza.csv"), row.names = FALSE
)
write.csv(
  data.frame(Accession = pop(gl), pca$scores),
  paste0(RUTA_CARPETA_SALIDA, "pca_scores.csv"), row.names = FALSE
)
cat("Varianza explicada, primeros 3 ejes del PCA:", round(varianza_pca[1:3], 2), "\n")
t0 <- reportar_tiempo("PCA (glPca)", t0)


# ============================================================
# 5a. DAPC replicando el enfoque del paper (N_PCS_GRANDE fijo)
# ============================================================

buscar_K_con_salvaguarda <- function(gl, n_pca, k_max, etiqueta, carpeta_salida) {
  grupos <- find.clusters(gl, n.pca = n_pca, max.n.clust = k_max, choose.n.clust = FALSE)
  
  # guardar la curva completa de BIC vs K -- el equivalente a la Fig S3
  # del paper, para poder inspeccionarla visualmente, no solo confiar
  # en la elección automática
  bic_df <- data.frame(K = seq_along(grupos$Kstat), BIC = as.numeric(grupos$Kstat))
  write.csv(bic_df, paste0(carpeta_salida, "bic_vs_K_", etiqueta, ".csv"), row.names = FALSE)
  png(paste0(carpeta_salida, "bic_vs_K_", etiqueta, ".png"), width = 800, height = 500)
  plot(bic_df$K, bic_df$BIC, type = "b", xlab = "K", ylab = "BIC",
       main = paste("BIC vs K --", etiqueta))
  dev.off()
  
  k_elegido <- length(unique(grupos$grp))
  if (k_elegido == k_max) {
    message("[AVISO SERIO] Para '", etiqueta, "', el K elegido automáticamente (",
            k_elegido, ") choca exactamente contra el techo (max.n.clust=", k_max,
            "). Esto normalmente significa que el BIC no encontró un mínimo real ",
            "dentro del rango probado -- NO confiar en este K sin revisar el .png ",
            "de la curva a mano. Puede que necesites subir max.n.clust, o elegir K ",
            "manualmente viendo dónde el codo/mínimo de la curva ocurre de verdad.")
  }
  grupos
}

max_pcs_posible <- min(nInd(gl), nLoc(gl)) - 1
if (N_PCS_GRANDE > max_pcs_posible) {
  message("[AVISO] N_PCS_GRANDE (", N_PCS_GRANDE, ") es mayor al máximo posible (",
          max_pcs_posible, "). Se usará el máximo posible en su lugar.")
  N_PCS_GRANDE <- max_pcs_posible
}

grupos <- buscar_K_con_salvaguarda(gl, N_PCS_GRANDE, K_MAX, "pcs_grande", RUTA_CARPETA_SALIDA)
t0 <- reportar_tiempo(paste0("find.clusters (", N_PCS_GRANDE, " PCs, K=1 a ", K_MAX, ")"), t0)

if (length(unique(grupos$grp)) < 2) {
  stop("K encontrado fue 1 (ningún grupo distinto) con N_PCS_GRANDE=", N_PCS_GRANDE,
       " -- el DAPC no puede correr con un solo grupo. Revisa bic_vs_K_pcs_grande.png; ",
       "puede que este número de PCs no esté capturando estructura real, o que de ",
       "verdad no haya estructura distinguible en este punto del pipeline.")
}

dapc_grande <- dapc(gl, pop = grupos$grp, n.pca = N_PCS_GRANDE, n.da = length(unique(grupos$grp)) - 1)

saveRDS(dapc_grande, paste0(RUTA_CARPETA_SALIDA, "dapc_grande_", N_PCS_GRANDE, "pcs.rds"))
cat("\nDAPC con", N_PCS_GRANDE, "PCs -- K encontrado:", length(unique(grupos$grp)), "\n")
cat("Proporción de reasignación correcta:", round(summary(dapc_grande)$assign.prop, 4), "\n")
t0 <- reportar_tiempo(paste0("dapc (", N_PCS_GRANDE, " PCs)"), t0)


# ============================================================
# 5b. DAPC con el número de PCs sugerido por xvalDapc
# ============================================================

cat("\n--- Corriendo xvalDapc (puede tardar) ---\n")
mat <- as.matrix(gl)
mat[is.na(mat)] <- 0  # xvalDapc no acepta NA -- imputación simple; revisar cuántos NA hay antes

xval <- xvalDapc(mat, pop(gl), n.pca.max = XVAL_N_PCA_MAX, n.rep = XVAL_N_REP,
                 xval.plot = FALSE)

n_pcs_xval <- as.numeric(xval$`Number of PCs Achieving Highest Mean Success`)
cat("Número de PCs sugerido por xvalDapc:", n_pcs_xval, "\n")
t0 <- reportar_tiempo(paste0("xvalDapc (hasta ", XVAL_N_PCA_MAX, " PCs, ", XVAL_N_REP, " repeticiones)"), t0)

grupos_xval <- buscar_K_con_salvaguarda(gl, n_pcs_xval, K_MAX, "pcs_xval", RUTA_CARPETA_SALIDA)

if (length(unique(grupos_xval$grp)) < 2) {
  stop("K encontrado fue 1 (ningún grupo distinto) con n_pcs_xval=", n_pcs_xval,
       " -- el DAPC no puede correr con un solo grupo. Revisa bic_vs_K_pcs_xval.png.")
}

dapc_xval <- dapc(gl, pop = grupos_xval$grp, n.pca = n_pcs_xval,
                  n.da = length(unique(grupos_xval$grp)) - 1)

saveRDS(dapc_xval, paste0(RUTA_CARPETA_SALIDA, "dapc_xval_", n_pcs_xval, "pcs.rds"))
cat("DAPC con", n_pcs_xval, "PCs (xvalDapc) -- K encontrado:", length(unique(grupos_xval$grp)), "\n")
cat("Proporción de reasignación correcta:", round(summary(dapc_xval)$assign.prop, 4), "\n")
t0 <- reportar_tiempo(paste0("find.clusters + dapc (", n_pcs_xval, " PCs, xvalDapc)"), t0)
cat(sprintf("\n[tiempo] TOTAL del script: %.2f min\n", as.numeric(difftime(Sys.time(), t_inicio_total, units = "mins"))))


# ============================================================
# 6. Comparación rápida para decidir qué número de PCs usar
#    en los siguientes checkpoints (paso 3, 4, 5)
# ============================================================

cat("\n=== COMPARACIÓN ===\n")
cat("PCs estilo paper (", N_PCS_GRANDE, "): K =", length(unique(grupos$grp)),
    ", reasignación =", round(summary(dapc_grande)$assign.prop, 4), "\n")
cat("PCs por xvalDapc (", n_pcs_xval, "): K =", length(unique(grupos_xval$grp)),
    ", reasignación =", round(summary(dapc_xval)$assign.prop, 4), "\n")
cat("\nSi la reasignación con PCs grande está muy cerca de 1.0 (100%) y la de\n")
cat("xvalDapc es notablemente más baja, es señal de sobreajuste con PCs grande\n")
cat("(como en la prueba con datos sintéticos). Decide aquí cuál número de PCs\n")
cat("usar de forma FIJA en los checkpoints 2 a 5, para no volver a correr\n")
cat("xvalDapc cada vez.\n")