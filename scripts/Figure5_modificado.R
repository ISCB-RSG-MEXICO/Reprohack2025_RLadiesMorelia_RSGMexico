# Título del manuscrito: Paisaje de elementos genéticos móviles y su carga de resistencia
# a antibióticos en genomas procariotas

# El siguiente código puede utilizarse para crear la Figura 5 del manuscrito

# Fecha de Creacion: 19-02-2021
# Autor: supriya.khedkar@embl.de

# Modificado por: Evelia Coss
# Fechad de modificacion: 6 de abril 2025

#---- Cargar paquetes -----

library(tidyverse)
library(cowplot)
library(viridis)
library(scales)

#---- PASO 1: Importar datos ----- 

arg_mge <- read_tsv("data/raw_data/mge_arg_final.txt.gz", col_names = T)
arg_genome <- read_tsv("data/raw_data/genome_arg.txt.gz", col_names = F)
glist <- read_tsv("data/raw_data/genome_status_supplementary_tableS2.txt.gz", col_names = T)

m_arg_rm <- read_tsv("data/raw_data/mge_arg_resistance_mechanism_final.txt.gz", col_names = T)
g_arg_rm <- read_tsv("data/raw_data/genome_arg_resistance_mechanism.txt.gz", col_names = T)
mge_solitary <- read_tsv("data/processed_data/solitary_mge_bins_final.txt.gz", col_names = T)
rm_list <- read_tsv("data/raw_data/resistance_mechanism.list.gz", col_names = T)

#---- PASO 2: Manipulación y Limpieza de los datos  ----- 

# Obtener los genomas con la mas alta calidad
glist_high <- glist %>% 
  filter(genome_quality == "high")

# ---- PASO 3: Procesamiento del dataframe arg_mge ----------

# 1. Filtrado de datos y obtencion de genomas
arg_mge <- arg_mge %>% 
dplyr::rename(abr_count = X3, total = X4) %>% 
mutate(mge = if_else(mge == "Cellular", "Genome", mge)) %>% 
mutate(island1 = island) %>% 
separate(island1, c("g1","g2","g3")) %>% 
mutate(genome = paste(g1, g2, sep = ".")) %>% 
select(-g1,-g2,-g3) %>% 
filter(genome %in% glist_high$genome) %>% 
select(-genome)

# 2. Procesamiento del dataframe arg_genome
arg_genome <- arg_genome %>% 
  dplyr::rename(mge = X1, island = X2, abr_count = X3, total = X4) %>% 
  mutate(abr_count = ifelse(abr_count > 0, abr_count, 0)) 

# 3. Combina ambos datasets
arg_combined <- rbind (arg_mge,arg_genome)

# 4. Obtiene el nombre del genoma de cada island
arg_combined_gNames <-  arg_combined %>%
  mutate(island1 = island) %>%
  separate(island1, c("one","two","three"), sep = "\\.", extra = "merge") %>%
  unite(.,genome, c("one","two"), sep = ".") %>%
  select(-three)

# ------ PASO 4: Filtrado y agrupación de ARGs en MGEs -----------------------
# 1. Filtrado de datos y obtencion de ARGs
arg_mge_filtered <- arg_mge %>% 
  mutate(island1 = island) %>% 
  separate(island1, c("g1","g2","g3")) %>% 
  mutate(genome = paste(g1, g2, sep = ".")) %>% 
  select(-g1,-g2,-g3) %>% 
  filter(genome %in% glist_high$genome) %>% 
  select(-genome) %>% 
  filter(.,total > 2) %>% 
  group_by(mge) %>% 
  summarise(abr_total = sum(abr_count), all = sum(total) - sum(abr_count))

# 2. Aplica lo mismo pero a arg_combined_gNames
arg_combined_filtered <- arg_combined_gNames %>% 
  filter(genomeID %in% glist_high$genome) %>% 
  select(-genomeID) %>% 
  filter(.,total > 2) %>% 
  group_by(mge) %>% 
  summarise(abr_total = sum(abr_count), all = sum(total) - sum(abr_count))

# 3. Preparar la tabla para análisis estadístico (pre_arg_combined_filtered_m)
pre_arg_combined_filtered_m <- arg_combined_filtered %>% 
  dplyr::rename(abr = abr_total, non_abr = all) %>% 
  reshape2::melt() %>% 
  spread(mge,value)

# ------ PASO 5: Ajustes de los valores de pvalue ---------------------------
#  1. Preparar nombres
allrec <- subset(colnames(pre_arg_combined_filtered_m), colnames(pre_arg_combined_filtered_m)!="variable") 
alldomains <- as.vector(pre_arg_combined_filtered_m$variable)
rownames(pre_arg_combined_filtered_m) <- as.vector(pre_arg_combined_filtered_m$variable)
arg_combined_filtered_m <- pre_arg_combined_filtered_m %>% select(-variable)

#  2. Test de Fisher
FE <- matrix(NA, nrow = length(alldomains), ncol = length(allrec),dimnames = list(alldomains, allrec))
OD <- matrix(NA, nrow = length(alldomains), ncol = length(allrec),dimnames = list(alldomains, allrec))
for (i in alldomains) {
  for(j in allrec) {
    print(paste0(i, "_", j))
    incog <- arg_combined_filtered_m [i,j]
    outcog <- sum(arg_combined_filtered_m [,j]) - incog
    notmge <- sum(arg_combined_filtered_m [i,]) - incog
    outcogoutmge <- sum(arg_combined_filtered_m ) - incog - outcog - notmge
    m <- matrix(data = c(incog, outcog, notmge, outcogoutmge), ncol = 2)
    test <- fisher.test(m,alternative = "greater")
    FE[i,j] <- test$p.value
    OD[i,j] <- test$estimate
  }
}

# 3. Corrección de Bonferroni
FE_adj <- matrix(p.adjust(p = FE, method = "bonferroni"), ncol = ncol(FE), nrow = nrow(FE), dimnames = list(rownames(FE), colnames(FE)))

# ------ PASO 6: Preparación de datos para graficar -----

# Se transforman las matrices de odds ratios y p-values ajustados a formato largo.
OD_mat <- matrix(OD, ncol = ncol(OD), nrow = nrow(OD), dimnames = list(rownames(OD), colnames(OD)))

od_melt <- OD_mat %>% 
  reshape2::melt() %>% 
  dplyr::rename(resistance = Var1, mge = Var2, score = value)

fe_melt <- FE_adj %>% 
  reshape2::melt() %>% 
  dplyr::rename(resistance = Var1, mge = Var2, adj_pval = value)

# Esto selecciona solo las filas para resistance == "abr" y clasifica MGEs según la significancia del test de enriquecimiento.
# Después se hace un left_join con los datos de abundancia de ARGs para obtener:
# abr_total = conteo total de genes de resistencia por MGE.
# signP = categorización de significancia.

merge_od_pval <- left_join(od_melt,fe_melt, by = c("resistance", "mge")) %>% 
  filter(.,resistance == "abr") %>% 
  mutate(signP=ifelse(adj_pval>0.05,"p value > 0.05","p value < 0.05"))

combine <- left_join(merge_od_pval,arg_combined_filtered, by = "mge")
abr_nplot <- ggplot(combine, aes(x = reorder(mge,-abr_total,mean), y = log10(abr_total), fill = as.factor(signP))) + geom_bar(stat = "identity") + scale_fill_manual(values=c('#E69F00','#999999')) + geom_text(data = combine, aes(x = mge, y = log10(abr_total), label = abr_total, angle = 90)) 
abr_nplot + theme_cowplot(font_size = 20) + theme(axis.text.x=element_text(angle = 90, hjust = 1, vjust = 0.5)) + labs (x = "", y = "#ABR genes" ,fill = "") 


### Figura 5B: Asociación de mecanismos de resistencia ARG

# 1. Renombrar columnas
mge_solitary <- mge_solitary %>% 
  dplyr::rename(IS_Tn = X1, Phage = X2, Phage_like = X3, CE = X4, Integron = X5, MI = X6,	Hotspot	= X7, UC = X8,	Cellular = X9, island = X10, island_size = X11, prot_count = X12,	phage_count = X13, CONJ_T4SS = X14,	mgeR = X15)

# 2. Filtrar y anotar m_arg_rm y g_arg_rm
m_arg_rm <- m_arg_rm %>% 
  dplyr::rename(island1 = ptn, prot = island) %>% 
  dplyr::rename(island = island1) %>% 
  mutate(island1 = island) %>% 
  separate(island1, c("g1","g2","g3")) %>% 
  mutate(genome = paste(g1, g2, sep = ".")) %>% 
  select(-g1,-g2,-g3) %>% 
  filter(genome %in% glist_high$genome) %>% 
  select(-genome)

#  3. Crear mge_non_nested
m_arg_rm %>% pull(prot) %>%  n_distinct()
g_arg_rm <- g_arg_rm %>% 
  dplyr::rename(island1 = ptn, prot = island) %>% 
  dplyr::rename(island = island1) %>% 
  mutate(mge = "Genome") %>% 
  mutate(island1 = island) %>% 
  separate(island1, c("g1","g2","g3")) %>% 
  mutate(genome = paste(g1, g2, sep = ".")) %>% 
  select(-g1,-g2,-g3) %>% 
  filter(genome %in% glist_high$genome) %>% 
  select(-genome)

mge_non_nested <- mge_solitary %>% 
  select(1:11,15) %>% 
  gather(mge, mge_pa, 1:9) %>% 
  mutate(island1 = island) %>% 
  separate(island1, c("g1","g2","g3")) %>% 
  mutate(genome = paste(g1, g2, sep = ".")) %>% 
  select(-g1,-g2,-g3) %>% 
  filter(genome %in% glist_high$genome) %>% 
  select(-genome) %>% 
  filter(mge_pa ==1) %>% 
  filter(.,!grepl("UC", mge)) %>% 
  filter(.,!grepl("Hotspot", mge)) %>% 
  mutate(mgeRn = str_replace_all(mgeR,"_","")) %>% 
  filter(!str_detect(mgeRn, '[:alnum:] &{1,}')) %>% 
  select(-mgeRn) %>% 
  mutate(mgeR = str_replace_all(mgeR,"&","")) %>% 
  select(-mgeR,-mge_pa)

# 4. Unir proteínas con sus MGEs
m_arg_rm_final <- left_join(m_arg_rm,mge_non_nested, by="island") %>% 
drop_na(mge)

## Figure 5B – Association of ARG Resistance Mechanisms with MGEs
# Paso 1: Limpiar nombres de categorías de mecanismos de resistencia en `rm_list`
rm_list <- rm_list %>% 
  mutate(abr_categories = str_replace_all(abr_categories,"\\s","_"))

# Paso 2: Limpiar y seleccionar columnas relevantes de ARGs en MGEs (`m_arg_rm_final`) y en genoma (`g_arg_rm`)
m_arg_rm_final <- m_arg_rm_final %>% 
  select(-confers_resistance_to,-island_size)
g_arg_rm_final <- g_arg_rm %>% 
  select(-confers_resistance_to)  

# Paso 3: Combinar ARGs en MGEs y en el genoma en un solo objeto (`together_rm`)
together_rm <- rbind(m_arg_rm_final,g_arg_rm_final) %>% 
  mutate(mge = ifelse(mge == "Cellular", "Genome", mge)) %>% 
  mutate(participates_in_1 = str_replace_all(participates_in,"^ARO:(\\d+)\\s!\\s",""), is_a_1 = str_replace_all(is_a,"^ARO:(\\d+)\\s!\\s","")) %>% 
  mutate(participates_in_1 = str_replace_all(participates_in_1,"\\s","_"), is_a_1 = str_replace_all(is_a_1,"\\s","_")) %>% 
  mutate(comb =if_else(participates_in_1 %in% rm_list$abr_categories, participates_in_1, is_a_1)) %>% 
  filter(!is.na(comb)) %>%
  filter(comb != "reduced_permeability_to_antibiotic" | comb != "GO:0015291_!_secondary_active_transmembrane_transporter_activity") %>% 
  select(island,mge,comb) %>%
#  group_by(mge,island,comb) %>%
  unique() %>%
#  ungroup() %>%
  group_by(mge,comb) %>%
  summarise(count = n()) %>% 
  filter(comb != "GO:0015291_!_secondary_active_transmembrane_transporter_activity")

# Paso 4: Preparar tabla final de mecanismos vs MGEs
# Renombrar columna para claridad
pre_arg_rm <- together_rm %>%
  filter(comb != "hydrolysis_of_antibiotic_conferring_resistance") %>%
  dplyr::rename(abr_categories = comb) %>% 
  spread(mge, count,fill = 0)

# Paso 5: Transformar a formato largo para visualizaciones posteriores
m_pre_arg <- pre_arg_rm %>% 
reshape2::melt()

# ---- Figure 5C – Enrichment of Resistance Mechanisms per MGE via Fisher's Exact Test----------
# Paso 1: Preparar matriz de conteo de mecanismos de resistencia vs tipos de MGE
arg_rm <- as.data.frame(pre_arg_rm)
allMges <- subset(colnames(arg_rm), colnames(arg_rm)!="abr_categories") 
allClasses <- as.vector(arg_rm$abr_categories)
rownames(arg_rm) <- as.vector(arg_rm$abr_categories)
arg_rm <- arg_rm %>% select(-abr_categories)

# Paso 2: Pruebas de Fisher para evaluar enriquecimiento de cada mecanismo en cada MGE
FE_arg <- matrix(NA, nrow = length(allClasses), ncol = length(allMges),dimnames = list(allClasses, allMges))
for (i in allClasses) {
  for(j in allMges) {
    print(paste0(i, "_", j))
    inarg <- arg_rm[i,j]
    outarg <- sum(arg_rm[,j]) - inarg
    notmge <- sum(arg_rm[i,]) - inarg
    outargoutmge <- sum(arg_rm) - inarg - outarg - notmge
    m <- matrix(data = c(inarg, notmge,outarg, outargoutmge), ncol = 2)
    test <- fisher.test(m, alternative = "greater")
    FE_arg[i,j] <- test$p.value
  }
}

# Paso 3: Corrección de valores p con Bonferroni
FE_arg_adj <- matrix(p.adjust(p = FE_arg, method = "bonferroni"), ncol = ncol(FE_arg), nrow = nrow(FE_arg), dimnames = list(rownames(FE_arg), colnames(FE_arg)))
# Paso 4: Preparar tabla larga con significancia y anotaciones
FE_arg_adj_t <- as.tibble(FE_arg_adj) %>% 
  mutate(abr = rownames(FE_arg_adj))

# Paso 5: Integrar con la matriz de conteo original para graficar
m_FE_arg_adj_t <- FE_arg_adj_t %>% 
  gather(mge,p_adj,1:7) %>% 
  mutate(signP=ifelse(p_adj>0.05,"","*"))


m_pre_arg_n <- m_pre_arg %>% 
  dplyr::rename(abr = abr_categories,mge = variable)
m_FE_arg_adj_t_rs <- left_join(m_FE_arg_adj_t, m_pre_arg_n, by = c("abr", "mge")) 

# Paso 6: Generar heatmap de mecanismos de resistencia vs MGEs
arg_hmap <- ggplot(m_FE_arg_adj_t_rs,aes(x = reorder(mge,value, sum),y = reorder(abr,value,sum), fill = log10(value))) + geom_tile(colour = "black") + geom_text(aes(label = signP), colour = "white", size = 8) + scale_fill_gradientn(colours = viridis(10), values = rescale(c(0,3,4,6)), na.value = "grey80")
# Paso 7: Personalizar tema de la figura
arg_hmap + theme_cowplot(font_size=20) + theme(axis.text.x=element_text(angle = 90, hjust = 1))  + labs(x = "", y = "", fill = "#ABR genes")
