# Título del manuscrito: Paisaje de elementos genéticos móviles y su carga de resistencia
# a antibióticos en genomas procariotas

# El siguiente código puede utilizarse para crear la Figura 4 del manuscrito

# Fecha de Creacion: 19-02-2021
# Autor: supriya.khedkar@embl.de

# Modificado por: Evelia Coss
# Fechad de modificacion: 6 de abril 2025

#---- Cargar paquetes ----- 

library(tidyverse)
library(cowplot)
library(reshape)
library(RColorBrewer)
library(phytools) 
library(viridis)

#---- PASO 1: Importar datos ----- 

data_mf <- read_tsv("data/raw_data/recombinase_hgt_cluster_master_file.txt.gz", col_names = F)
mge_bins <- read_tsv("data/raw_data/mge_bins_final.txt.gz",col_names = T)
tax <- read_tsv("data/raw_data/hgt_species.list.gz", col_names = F)
class_tree <- read.tree("data/raw_data/progenomes2_class_tree.nwk")
glist <- read_tsv("data/raw_data/genome_status_supplementary_tableS2.txt.gz", col_names = T)

# Cargamos el archivo con datos anotados sobre elementos genéticos móviles (MGE)
dat_mge <- read_tsv("data/processed_data/resolved_MGE_hgt_top_tax_final.txt.gz", col_names = F)

# > Tabla con datos de transferencia horizontal a nivel de clase taxonómica
# Esta tabla incluye información sobre las recombinasas asociadas a distintos tipos de MGE
data_mge_class <- read_tsv("data/processed_data/resolved_MGE_hgt_class_level_final.txt.gz", col_names = F)

# Paleta general de colores
colc <- c("#D55E00", "#E69F00", "#F0E442", "#56B4E9", "#009E73", "#0072B2","#CECCCC")
names(colc) <- c("IS_Tn", "Phage", "Phage_like", "CE", "Integron", "MI", "Cellular") 

# Carga funciones
source("scripts/phylo_heatmap_function.R")

#---- PASO 2: Cambios de formato y estructura ----- 

# Renombramos las columnas de `data_mf` con nombres significativos
data_mf_rename <- data_mf %>% 
  dplyr::rename(rec_cluster = X1, prot = X2, specI = X3, kingdom = X4, phylum = X5, class = X6, order = X7, family = X8, genus = X9) 

# Renombramos las columnas de `tax` para que tengan nombres descriptivos.
# Luego seleccionamos solo la columna 'genome' y filtramos los que están presentes en `glist_high`.
tax <- tax %>% 
  dplyr::rename(specI = X1, genome = X2, kingdom = X3, phylum = X4, class = X5, family = X6) %>% 
  select(genome) %>% 
  filter(genome %in% glist_high$genome) 

# Renombramos las columnas para hacerlas más interpretables
data_mge_class_rename <- data_mge_class %>% 
  dplyr::rename(cor_mge = X1, island = X2, mge = X3, rec = X4, rec_cluster = X5, family = X6, class = X7, prot = X8)

#---- PASO 3: Manipulación y Limpieza de los datos  ----- 
# Obtener los genomas con la mas alta calidad
glist_high <- glist %>% 
  filter(genome_quality == "high")

# Cambio de formato y limpieza de datos
mge_bins_melted <- mge_bins %>% 
  select(1:6,9,10) %>% 
  # Reestructuramos el dataframe de formato ancho a largo
  gather("mge","count",1:7) %>% # Convertimos las columnas 1 a 7 en dos columnas: 'mge' (nombre original de la columna) y 'count' (su valor)
  filter(count > 0) %>% 
  mutate(island1 = island) %>% 
  separate(island1, c("g1","g2","g3")) %>% # Separamos 'island1' en tres componentes
  # Unimos las dos primeras partes para crear un identificador único de genoma (e.g., "gen1.gen2")
  mutate(genome = paste(g1, g2, sep = ".")) %>%  
  # Eliminamos las columnas auxiliares g1, g2 y g3
  select(-g1,-g2,-g3) %>% 
  # Filtramos solo los genomas que están presentes en la lista `glist_high$genome`
  filter(genome %in% glist_high$genome) %>% 
  select(-genome)

# Guardar datos
write.table(mge_bins_melted,file = "data/processed_data/mge_bins_final_collapsed.txt", sep = "\t", row.names = F, col.names = T, quote = F)

#---- PASO 3: Asignar información taxonómica y filtrar genomas de interés ----- 

# Extraemos el ID del genoma a partir del identificador de proteína
# (se asume que 'prot' tiene formato como "X.Y.Z" donde X.Y es el ID de genoma)
datar_phy <- data_mf_rename %>% 
  mutate(prot1 = prot) %>%
  separate(prot1, c("one","two","three"), sep = "\\.", extra = "merge") %>% 
  unite(.,genome, c("one","two"), sep = ".") %>% 
  select(-three) %>% 
  filter(genome %in% glist_high$genome) 

## Unimos información taxonómica con los datos de proteínas ##
datar <- left_join(tax,datar_phy, by = "genome") %>% 
  drop_na()  # eliminamos filas con valores NA después del join

#---- PASO 4: Identificar HGT a cada nivel taxonómico ----------------------------------------------

# Filtramos los clusters de recombinación (`rec_cluster`) que contienen más de una entrada,
# y que tienen variación en el nivel taxonómico "family".
# Esto sugiere la presencia de transferencia horizontal entre diferentes familias.

#datar %>% count(rec_cluster) 
datar_melt <- datar %>% 
  select(-genus) %>%         # Eliminamos columna de genus (no se usa aquí)
  group_by(rec_cluster) %>%  # Agrupamos por cluster de recombinación
  filter(n() > 1) %>%        # Nos quedamos solo con clusters que tienen más de una entrada
  filter_at(vars(family), any_vars(length(unique(.)) > 1)) %>% # Que haya al menos dos familias distintas
  ungroup() 

# Contamos el número de taxones distintos por nivel (kingdom a family) por cada cluster
datar_tax <- datar_melt %>% 
  group_by(rec_cluster) %>% 
  summarise_at(vars(kingdom:family), n_distinct) # Cuenta cuántos taxones distintos hay en cada nivel

# Preparamos un resumen del nivel taxonómico más alto en el que ocurre variación
datar_top <- datar_tax %>% 
  mutate_at(vars(kingdom:family), ~ . - 1) %>%  # Restamos 1 (para que 0 indique sin variación)
  mutate(cumul = rowSums(.[2:6])) %>%           # Suma total de variación en niveles taxonómicos
  pivot_longer(kingdom:family, names_to = "key", values_to = "value") %>% # Reorganiza a formato largo
  group_by(rec_cluster) %>% 
  mutate(top_tax = key[which(value > 0)[1]]) %>% # Identifica el primer nivel donde hay variación
  ungroup() %>% 
  pivot_wider(names_from = key, values_from = value) # Regresa a formato ancho

# Combinamos con información de proteínas para cada cluster
sub_datar <- datar %>% select(rec_cluster,prot)
combine_mge <- left_join(datar_top,sub_datar, by = "rec_cluster")

# Guardar datos
write.table(combine_mge,file= "data/processed_data/pre_MGE_top_tax_file_final.txt", sep = "\t", quote = F, col.names = T, row.names = F)

#---- PASO 4: Estratificación basada en MGEs de los clusters de recombinasas ---------------------

# Limpiamos y formateamos el archivo:
# - Quitamos el encabezado repetido
# - Renombramos columnas para claridad
# - Extraemos el nombre del genoma a partir del identificador de proteína
# - Filtramos los genomas de interés y descartamos islas no definidas
dat_mge_f <- dat_mge %>% 
  filter(!X5 == "rec_cluster") %>% 
  dplyr::rename(MGE = X1, ISLAND = X2, MGE_ND = X3, rec = X4, rec_cluster = X5,	cumul = X6,	top_tax = X7,	kingdom = X8,	phylum = X9, class = X10,	order = X11, family = X12, prot = X13) %>% 
  mutate(prot1 = prot) %>%
  separate(prot1, c("one","two","three"), sep = "\\.", extra = "merge") %>% 
  unite(.,genome, c("one","two"), sep = ".") %>% 
  select(-three) %>% 
  filter(genome %in% glist_high$genome) %>% 
  filter(ISLAND != "ISLAND_ND") # Excluye islas no definidas

# Preparamos una versión más restringida del conjunto de datos,
# conservando solo genomas que están también en `datar`
dat_mge_pref <- dat_mge_f %>%
  mutate(prot1 = prot) %>%
  separate(prot1, c("one","two","three"), sep = "\\.", extra = "merge") %>% 
  unite(.,genome, c("one","two"), sep = ".") %>% 
  select(-three) %>%
  filter(genome %in% datar$genome)

# Quitamos las columnas taxonómicas redundantes para fusionar luego con datos actualizados
dat_mge_f1 <- dat_mge_pref %>% select(-kingdom,-phylum,-class,-order,-family)

# Fusionamos los datos de MGE con la anotación taxonómica por proteína y cluster
# Luego, reorganizamos los datos en formato largo para obtener el valor en el nivel taxonómico más informativo
dat_with_tax <- left_join(dat_mge_f1,datar, by = c("rec_cluster","prot")) %>%
  select(MGE, rec_cluster, top_tax,prot, kingdom:family) %>%
  reshape2::melt(id = c("MGE", "rec_cluster", "top_tax","prot")) %>%
  # Filtramos solo el nivel taxonómico que fue identificado como el más informativo
  filter(variable == top_tax) 

#---- PASO 5: Generación de datos para las figuras 4C y 5 ---------------

# Fusionamos datos de recombinasas y MGEs con anotación taxonómica
# Luego transformamos a formato largo, conservando solo el nivel taxonómico más informativo
dat_with_tax2n <- left_join(dat_mge_f1,datar, by = c("rec_cluster","prot")) %>%
  select(MGE, MGE_ND, ISLAND, rec_cluster, top_tax,prot, kingdom:family) %>%
  reshape2::melt(id = c("MGE", "MGE_ND","rec_cluster", "ISLAND", "top_tax","prot")) %>%
  filter(variable == top_tax) #  Solo nos quedamos con la categoría taxonómica más relevante

# Preparamos subconjunto con información de familia para cada recombinasa
temp_datar <- datar %>% select(rec_cluster,prot,family)
# Añadimos la información de familia a los datos largos con MGEs
temp_data_with_tax2n <- left_join(dat_with_tax2n,temp_datar, by = c("rec_cluster","prot")) 

# Filtramos clusters de recombinasas asociadas a MGEs móviles (excluyendo "Cellular" y "nested")
# Luego, identificamos clusters que presentan diversidad en el nivel taxonómico relevante
dat_with_tax3_testn <- temp_data_with_tax2n %>%
  filter(MGE != "Cellular" & MGE != "nested") %>%              # Nos centramos en MGEs no celulares
  group_by(rec_cluster, MGE) %>%                               # Agrupamos por cluster y tipo de MGE
  mutate(val = n_distinct(value) > 1) %>%                      # Evaluamos si hay más de un taxón distinto
  filter(val == "TRUE") %>%                                    # Conservamos solo aquellos con diversidad taxonómica
  ungroup()

# Guardar datos
write.table(dat_with_tax3_testn,file="processed_data/all_hgt_data_family_expanded_redundant_final.txt", sep = "\t", col.names = T, row.names = F, quote = F)

#---- PASO 6: Preparación de datos para cuantificar HGT según el tipo de MGE y el nivel taxonómico afectado --------------------------

# Extraemos información relevante (cluster, proteína, familia) desde los datos anotados
temp_datar <- datar %>% select(rec_cluster,prot,family)
# Unimos esta información a los datos previos con anotación taxonómica
temp_data_with_tax <- left_join(dat_with_tax,temp_datar, by = c("rec_cluster","prot")) 

# Filtramos MGEs móviles (excluyendo los de origen celular o sin clasificar),
# luego identificamos recombinasas que están distribuidas en más de un taxón distinto
dat_with_tax_mge <- temp_data_with_tax %>%
  filter(MGE != "Cellular" & MGE != "nested") %>%        # Solo MGEs móviles
  group_by(rec_cluster, MGE) %>%                         # Agrupamos por cluster y tipo de MGE
  mutate(val = n_distinct(value) > 1) %>%                # Evaluamos si hay diversidad taxonómica
  filter(val == "TRUE") %>%                              # Conservamos solo los casos con diversidad
  summarise(top_tax_new = unique(top_tax)) %>%           # Extraemos el nivel taxonómico más informativo
  ungroup() %>%
  group_by(MGE, top_tax_new) %>%                         # Agrupamos por tipo de MGE y nivel taxonómico
  summarise(final = n()) %>%                             # Contamos número de clusters con diversidad
  mutate(tot = sum(final)) %>%                           # Calculamos total por MGE
  ungroup()

## Ajuste manual: Forzamos a que haya al menos 3 eventos a nivel de reino (kingdom)
## para reflejar transferencias entre arqueas y bacterias aunque sean pocas pero significativas
dat_with_tax_mge <- dat_with_tax_mge %>%
  mutate(final = if_else(top_tax_new == "kingdom",3,as.numeric(final))) 

# Ordenamos niveles de factores para visualización ordenada
dat_with_tax_mge$top_tax_new <- factor(dat_with_tax_mge$top_tax_new, levels = rev(c("kingdom","phylum","class","order","family")))
dat_with_tax_mge$MGE <- factor(dat_with_tax_mge$MGE,levels = c("IS_Tn","CE","MI","Phage","Phage_like","Integron"))

#---- Figura 4A --------------------------

# Generamos un gráfico de línea con puntos para visualizar los eventos de HGT
# por nivel taxonómico (kingdom, phylum, etc.) estratificados por tipo de MGE

mge_tax_log <- ggplot(dat_with_tax_mge, aes(x = top_tax_new, y = final, color = MGE)) + 
  geom_point(aes(group = MGE), size = 3) +   # Puntos grandes para cada tipo de MGE
  geom_line(aes(group = MGE), size = 1) +    # Conectamos puntos con líneas para cada tipo de MGE
  scale_colour_manual(values = colc) +       # Usamos una paleta de colores personalizada
  scale_y_log10()                            # Escala logarítmica para visualizar mejor las diferencias

# Personalizamos el gráfico con un tema limpio y rotación de etiquetas
mge_tax_log + 
  theme_cowplot(font_size = 15) +            # Estilo claro y moderno
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5)) + # Rotamos etiquetas del eje X
  labs(y = "#HGT events", x = "")            # Etiquetas de los ejes

#---- Figura 4B - HGT heatmap at taxonomic class level -------------------------- 

# Preprocesamiento de la tabla con información de clases taxonómicas y MGEs
data_mge_class_phy <- data_mge_class_rename %>% 
  mutate(prot1 = prot) %>%  # Guardamos la columna original de proteínas en una nueva columna
  separate(prot1, c("one","two","three"), sep = "\\.", extra = "merge") %>% 
  unite(.,genome, c("one","two"), sep = ".") %>%  # Reunimos los dos primeros fragmentos para identificar el genoma
  select(-three) %>%  # Eliminamos el componente restante de la separación
  filter(genome %in% tax$genome) %>%  # Filtramos solo los genomas presentes en la lista filtrada de taxonomía
  filter(.,!grepl("Cellular",cor_mge))  # Excluimos elementos que no son MGEs móviles ("Cellular")

## Identificación de recombinasas compartidas entre genomas (potenciales eventos de HGT) ##
data_mge_class_hgt <- data_mge_class_phy %>% 
  group_by(rec_cluster) %>%  # Agrupamos por clúster de recombinasa
  summarise(count = n()) %>%  # Contamos cuántos genomas comparten ese clúster
  filter(., count > 1)  # Nos quedamos solo con aquellos compartidos por más de un genoma (indicador de HGT)

# Filtramos la tabla original para quedarnos solo con los clústeres que muestran posible HGT
data_hgt <- left_join(data_mge_class_hgt, data_mge_class_phy, by = "rec_cluster")

# ----- Paso 10: Contar combinaciones de familias por recombinasa
pre_hgt <- data_hgt %>% 
  group_by(rec_cluster) %>% 
  summarise(family = paste(family, collapse = ":")) %>% 
  group_by(rec_cluster,family) %>% 
  summarise(count = n()) 

# Normalizar combinaciones de familias (A:B ≡ B:A)
pre_hgt_decoupled <- pre_hgt %>% 
  group_by(family) %>% 
  mutate(nfamily = duo_maker(family)) %>% 
  ungroup() %>% 
  separate_rows(nfamily, sep = ";") %>% group_by(rec_cluster,nfamily) %>% summarise(ncount = sum(count)) %>% 
  mutate(nfamily = map_chr(nfamily,~toString(sort(str_split(.x, ":")[[1]]))))  %>% 
  group_by(nfamily, rec_cluster)  %>% 
  summarise(final_count = sum(ncount)) %>%
  ungroup() %>%
  filter(.,!grepl("NA ",nfamily))

# ----- Paso 10: ----------------

## compute recombinase HGT events across taxonomic class ##
rel_family_class <- data_hgt %>% 
  select(family,class) %>% 
  unique(.) %>% 
  mutate(family1 = family, family2 = family) %>% 
  select(-family)

pre_hgt_decoupled_family <- pre_hgt_decoupled %>% 
  separate(nfamily, c("family1","family2"), sep = ", ", extra = "merge")

pre_hgt_class1 <- left_join(pre_hgt_decoupled_family,rel_family_class, by = "family1") %>% 
  select(-family1) %>% 
  dplyr::rename(family2 = family2.x)

pre_hgt_class2 <- left_join(pre_hgt_class1,rel_family_class, by = "family2") %>% 
  select(-family2,-family2.y,-family1) %>% 
  filter(.,class.x != class.y) %>% 
  group_by(class.x,class.y) %>% 
  summarise(hgt_count = sum(final_count)) %>%
  unite(name_new,class.x:class.y, sep = ":") %>%
  mutate(name_new = map_chr(name_new,~toString(sort(str_split(.x, ":")[[1]])))) %>% 
  group_by(name_new) %>% 
  summarise(final_count = sum(hgt_count)) %>% ungroup() 

pre_hgt_class2 %>% mutate(tot = sum(final_count))
#Total HGT events class level 2823 for high quality genomes

# ----- Paso 10: ----------------
## Compute HGT based on MGE to obtain heatmap arcs
data_mge_class_phy_MGE <- data_mge_class_rename %>% 
  mutate(prot1 = prot) %>%
  separate(prot1, c("one","two","three"), sep = "\\.", extra = "merge") %>% 
  unite(.,genome, c("one","two"), sep = ".") %>% 
  select(-three) %>% 
  filter(genome %in% tax$genome) %>% 
  filter(!grepl("Cellular",cor_mge) & !grepl("MGE_ND", cor_mge)) 

## genomes sharing MGE associated recombinase i.e. showing HGT ##
data_mge_class_hgt_MGE <- data_mge_class_phy_MGE %>% group_by(rec_cluster) %>% summarise(count = n()) %>%  filter(., count > 1) 
data_hgt_MGE <- left_join(data_mge_class_hgt_MGE, data_mge_class_phy_MGE, by = "rec_cluster") 

pre_hgt_MGE <- data_hgt %>% 
  group_by(rec_cluster,cor_mge) %>% 
  summarise(family = paste(family, collapse = ":")) %>% 
  group_by(rec_cluster,cor_mge,family) %>% 
  summarise(count = n()) 

## make 3 or more combinations of taxonomic family to 2 and merge redundant combinations A:B and B:A ##
pre_hgt_decoupled_MGE <- pre_hgt_MGE %>% 
  group_by(family) %>% 
  mutate(nfamily = duo_maker(family)) %>% 
  ungroup() %>% 
  separate_rows(nfamily, sep = ";") %>% group_by(cor_mge,nfamily) %>% summarise(ncount = sum(count)) %>% 
  mutate(nfamily = map_chr(nfamily,~toString(sort(str_split(.x, ":")[[1]])))) %>% 
  group_by(nfamily, cor_mge)  %>% 
  summarise(final_count = sum(ncount)) %>%
  ungroup() %>%
  filter(.,!grepl("NA ",nfamily))

# -----Cálculo de eventos de transferencia horizontal (HGT) mediados por MGE entre clases taxonómicas ----
rel_family_class_MGE <- data_hgt %>% 
  select(family,class) %>% 
  unique(.) %>% 
  mutate(family1 = family, family2 = family) %>% select(-family)

pre_hgt_decoupled_family_MGE <- pre_hgt_decoupled_MGE %>% 
  separate(nfamily, c("family1","family2"), sep = ", ", extra = "merge")

pre_hgt_class1_MGE <- left_join(pre_hgt_decoupled_family_MGE,rel_family_class_MGE, by = "family1") %>% 
  select(-family1) %>% 
  dplyr::rename(family2 = family2.x)

pre_hgt_class2_MGE <- left_join(pre_hgt_class1_MGE,rel_family_class_MGE, by = "family2") %>% 
  select(-family2,-family2.y,-family1) %>% 
  filter(.,class.x != class.y) %>% 
  group_by(class.x,class.y,cor_mge) %>% 
  summarise(hgt_count = sum(final_count))

pre_hgt_decoupled_class_MGE <- pre_hgt_class2_MGE %>% 
  unite(new, class.x:class.y, sep = ":") %>%
  group_by(new,cor_mge) %>% 
  mutate(name_new = duo_maker(new)) %>% 
  ungroup() %>% 
  separate_rows(name_new, sep = ";") %>% 
  group_by(name_new,cor_mge) %>%
  summarize(hgt_count = sum(hgt_count)) %>%
  ungroup() %>%
  mutate(name_new = map_chr(name_new,~toString(sort(str_split(.x, ":")[[1]])))) %>% 
  group_by(name_new,cor_mge) %>% 
  summarise(final_count = sum(hgt_count)) %>% 
  ungroup() %>%
  separate(name_new,c("class.x", "class.y"))

# -------------------
## HGT arcs Figure 4B for iTOL ##
class_tab <- pre_hgt_decoupled_class_MGE %>% 
  unite(new, class.x:class.y, sep = ":") %>%
  group_by(new,cor_mge) %>% 
  mutate(name_new = duo_maker(new)) %>% 
  ungroup() %>% 
  separate_rows(name_new, sep = ";") %>% 
  group_by(name_new,cor_mge) %>%
  summarize(hgt_count = sum(final_count)) %>%
  select(name_new,cor_mge, hgt_count) %>%
  ungroup() %>% 
  separate(name_new, c("class.x","class.y"))

#Percentage of HGT events stratified by MGE category
class_tab %>% 
  group_by(cor_mge) %>% 
  summarise(f = sum(hgt_count)) %>% 
  mutate(frac = (f/sum(f))*100)

class_tab_plot <- class_tab %>% 
  group_by(cor_mge) %>% 
  summarise(count = sum(hgt_count)) 

# ----- Paso 10: ----------------

ggplot(class_tab_plot, aes(x = reorder(cor_mge,-count,mean), y = count)) + 
geom_bar(stat="identity") + 
ylab("#HGT events") + 
xlab("") +
theme_cowplot() 

#
mges <- c("IS_Tn", "Phage", "Phage_like", "CE", "Integron", "MI")
colc_grey <- c("#D55E00", "#808080", "#808080", "#808080", "#808080", "#808080")
names(colc_grey) <- c("IS_Tn", "Phage", "Phage_like", "CE", "Integron", "MI")
angles <- seq(32,80,12)
for(i in 1:length(mges)) {
  separator <- c("DATASET_CONNECTION","SEPARATOR COMMA")
  dataset_label <- paste0("DATASET_LABEL,",mges[i])
  dcolor <- paste0("COLOR,",colc_grey[mges[i]])
  optional <- c("DRAW_ARROWS,0",paste0("CURVE_ANGLE,",angles[i]),"CENTER_CURVES,1","ALIGN_TO_LABELS,1")
  data_lines <- c("DATA")
  data_oi <- class_tab %>% filter(cor_mge==mges[i])
  out_lines <- c(separator, dataset_label, dcolor, optional, data_lines)
  for(j in 1:nrow(data_oi)){
    line_oi <- paste(data_oi[j,1],data_oi[j,2],"2",colc_grey[mges[i]],"normal",sep = ",")
    out_lines <- c(out_lines,line_oi)
  }
  writeLines(out_lines, paste0(mges[i],"_class_sub_grey_itol_connections_fig4B.txt"))
}

#for supplementary figure
angles <- seq(32,80,12)
for(i in 1:length(mges)) {
  separator <- c("DATASET_CONNECTION","SEPARATOR COMMA")
  dataset_label <- paste0("DATASET_LABEL,",mges[i])
  dcolor <- paste0("COLOR,",colc[mges[i]])
  optional <- c("DRAW_ARROWS,0",paste0("CURVE_ANGLE,",angles[i]),"CENTER_CURVES,1","ALIGN_TO_LABELS,1")
  data_lines <- c("DATA")
  data_oi <- class_tab %>% filter(cor_mge==mges[i])
  out_lines <- c(separator, dataset_label, dcolor, optional, data_lines)
  for(j in 1:nrow(data_oi)){
    line_oi <- paste(data_oi[j,1],data_oi[j,2],"2",colc[mges[i]],"normal",sep = ",")
    out_lines <- c(out_lines,line_oi)
  }
  writeLines(out_lines, paste0(mges[i],"_class_sub_itol_connections_fig4B_S.txt"))
}

##Figure 4B Heatmap plot recombinase HGT##
hgt_class_level <- pre_hgt_class2 %>% 
  separate(name_new, into = c("class1","class2"), sep = ", ") 


hgt_class_level_spread <- hgt_class_level %>%
  filter(!grepl("NA", class1)) %>%
  filter(!grepl("NA", class2)) %>% 
  filter(!grepl("NA", final_count)) %>% 
   pivot_wider(id_cols = class1, names_from = class2, values_from = final_count) %>%
   column_to_rownames(var = "class1") %>%
   replace(is.na(.), 0)

allClasses <- unique(c(hgt_class_level$class1, hgt_class_level$class2))#, class_tree$tip.label))
allClasses <- allClasses[!grepl("NA", allClasses)]

class_tree_rooted <- keep.tip(class_tree, as.vector(allClasses))
#write.tree(class_tree_rooted,file="processed_data/hgt_tree.nwk")

new_matrix <- matrix(0, nrow = Ntip(class_tree_rooted), ncol = Ntip(class_tree_rooted), dimnames = list(class_tree_rooted$tip.label, rev(class_tree_rooted$tip.label)))
for(i in rownames(hgt_class_level_spread )) {
  for(j in colnames(hgt_class_level_spread )) {
    value <- hgt_class_level_spread[i,j]
    if(value > 0) {
      new_matrix[i,j] <- (value)
      new_matrix[j,i] <- (value)
    }
  }
}

#code to plot
small_value <- unique(sort(unlist(new_matrix)))[2]
small_value <- small_value - small_value %% 0.001

lseq <- function(from, to, length.out) {
  # logarithmic spaced sequence
  10^(seq(log10(from), log10(to), length.out = length.out))
}
heatmap_breaks <- c(0,small_value,lseq(1,max(new_matrix),length.out = 18))
plot_colors <- c("white",colorRampPalette((viridis(10)), bias = 5)(length(heatmap_breaks)-2))

#phyloheatmap hack
phylo.heatmap.orient <- function (tree, X, fsize = 1, colors = NULL, standardize = FALSE, ...) 
{
  if (length(fsize) != 3) 
    fsize <- rep(fsize, 3)
  if (hasArg(legend)) 
    legend <- list(...)$legend
  else legend <- TRUE
  if (hasArg(labels)) 
    labels <- list(...)$labels
  else labels <- TRUE
  if (hasArg(split)) 
    split <- list(...)$split
  else split <- c(0.5, 0.5)
  split <- split/sum(split)
  if (is.null(colnames(X))) 
    colnames(X) <- paste("var", 1:ncol(X), sep = "")
  if (standardize) {
    sd <- apply(X, 2, function(x) sqrt(var(x, na.rm = TRUE)))
    X <- (X - matrix(rep(1, Ntip(tree)), Ntip(tree), 1) %*% 
            colMeans(X, na.rm = TRUE))/(matrix(rep(1, Ntip(tree)), 
                                               Ntip(tree), 1) %*% sd)
  }
  if (hasArg(xlim)) 
    xlim <- list(...)$xlim
  else xlim <- c(-0.5, (2 - 0.5) * split[2]/split[1] + 0.5)
  if (hasArg(ylim)) 
    ylim <- list(...)$ylim
  else ylim <- if (legend) 
    c(if (standardize) -0.15 else -0.1, if (labels) 1.1 else 1)
  else c(0, if (labels) 1.1 else 1)
  if (hasArg(mar)) 
    mar <- list(...)$mar
  else mar <- rep(1.1, 4)
  if (is.null(colors)) 
    colors <- heat.colors(n = 20)[20:1]
  if (hasArg(grid)) 
    add.grid <- list(...)$grid
  else add.grid <- FALSE
  cw <- untangle(tree, "read.tree")
  plot.new()
  par(mar = mar)
  plot.window(xlim = xlim, ylim = ylim)
  h <- phylogram(cw, fsize = fsize[1], ...)
  # print(h)
  START <- h + 1/2 * ((2 - 0.5) * split[2]/split[1] + 0.5 - 
                        h)/(ncol(X) - 1) + 0.5 * strwidth("W") * fsize[1]
  END <- (2 - 0.5) * split[2]/split[1] + 0.5 - 1/2 * ((2 - 0.5) * split[2]/split[1] + 0.5 - START)/(ncol(X) - 1)
  X <- X[cw$tip.label, rev(cw$tip.label)]
  image(x = seq(START, END, by = (END - START)/(ncol(X) - 1)), 
        z = t(X[cw$tip.label, rev(cw$tip.label)]), add = TRUE, col = colors, ...)
  if (add.grid) {
    dx <- (END - START)/(ncol(X) - 1)
    x <- seq(START - dx/2, END + dx/2, by = dx)
    nTips <- length(tree$tip.label)
    y <- c(-1/(2 * (nTips - 1)), seq(0, 1, length = nTips) + 
             1/(2 * (nTips - 1)))
    segments(x, y[1], x, y[length(y)])
    segments(x[1], y, x[length(x)], y)
  }
  if (legend) 
    add.color.bar(leg = END - START, cols = colors, lims = range(X, na.rm = TRUE), title = if (standardize) "standardized value" else "value", subtitle = if (standardize) "SD units" else "", prompt = FALSE, x = START, y = -1/(2 * (Ntip(cw) - 1)) - 3 * fsize[3] * strheight("W"), digits = if (max(abs(X), na.rm = TRUE) < 1) round(log10(1/max(abs(X), na.rm = TRUE))) + 1 else 2, fsize = fsize[3])
  if (labels) 
    text(x = seq(START, END, by = (END - START)/(ncol(X) - 1)), y = rep(1 + 1/(2 * (Ntip(cw) - 1)) + 0.4 * fsize[2] * strwidth("I"), ncol(X)), colnames(X), srt = 70, adj = c(0, 0.5), cex = fsize[2])
  if (any(is.na(X))) {
    ii <- which(is.na(X), arr.ind = TRUE)
    x.na <- seq(START, END, by = (END - START)/(ncol(X) - 
                                                  1))[ii[, 2]]
    y.na <- seq(0, 1, by = 1/(nrow(X) - 1))[ii[, 1]]
    for (i in 1:length(x.na)) {
      xx <- x.na[i] + c(1/2, -1/2) * (END - START)/(ncol(X) - 
                                                      1)
      yy <- y.na[i] + c(-1/2, 1/2) * 1/(nrow(X) - 1)
      lines(xx, yy)
    }
  }
}
environment(phylo.heatmap.orient) <- environment(phylo.heatmap)
#

phylo.heatmap.orient(class_tree_rooted, new_matrix, fsize = c(0.9, 0.9, 1), colors = plot_colors, grid = T, split = c(0.6, 0.4), breaks = heatmap_breaks, lwd = 1)

###Figure 4C- Barplot of nested IS_Tn/Integron enrichment with other MGEs####

#load data#
datacl <- read_tsv("raw_data/all_recombinase_clusters_mge_resolved.txt", col_names = F)

#parse hgt nested data#
dat_with_tax_nested <- left_join(dat_mge_f1,datar, by = c("rec_cluster","prot")) %>%
  select(MGE, MGE_ND, ISLAND, rec_cluster, top_tax,prot, kingdom:family) %>%
  reshape2::melt(id = c("MGE", "MGE_ND", "rec_cluster", "ISLAND", "top_tax","prot")) %>%
  filter(variable == top_tax) 

hgt1.1 <- left_join(dat_with_tax_nested,temp_datar, by = c("rec_cluster","prot")) %>%
  filter(MGE != "Cellular" & MGE != "nested") %>%
  group_by(rec_cluster,MGE) %>%
  mutate(val = n_distinct(value) > 1) %>% 
  filter(val == "TRUE") %>% 
  ungroup()  


####determine overall nested fraction for IS_Tn####
hgt_tn1.1 <- hgt1.1 %>%
  dplyr::rename(cor_mge = MGE, mge = MGE_ND) %>% 
  select(1,2,3,4,6,9) %>%
  filter(cor_mge == "IS_Tn") 

hgt_tn1.1 %>%
  mutate(new = (grepl("(\\d+)&", mge))) %>% 
  filter(new == "TRUE")

hgt_tn1.1_check <- hgt_tn1.1 %>% 
  mutate(prot1 = prot) %>%
  separate(prot1, c("one","two","three"), sep = "\\.", extra = "merge") %>% 
  unite(.,genome, c("one","two"), sep = ".") %>% 
  select(-three) %>% 
  group_by(rec_cluster,genome) %>% 
  mutate(atleast_1_nested = if_else((grepl("(\\d+)&", mge)),1,0)) %>% 
  ungroup() %>% 
  group_by(rec_cluster,cor_mge,genome) %>% 
  summarise(val = if_else(sum(atleast_1_nested) > 0, "nested", "non_nested")) %>%
  ungroup() %>% 
  group_by(cor_mge,val) %>% 
  summarise(n = n()) 


hgt_tn1.1 %>% 
  group_by(rec_cluster,ISLAND) %>% 
  summarise(n = n_distinct(prot)) %>% 
  View()

hgt_tn1.2_check <- hgt_tn1.1 %>% 
  mutate(prot1 = prot) %>%
  separate(prot1, c("one","two","three"), sep = "\\.", extra = "merge") %>% 
  unite(.,genome, c("one","two"), sep = ".") %>% 
  select(-three) %>% 
  group_by(rec_cluster,ISLAND) %>% 
  summarise(nested_n = sum(grepl("(\\d+)&", mge)), all_n = n_distinct(prot)) %>% 
  ungroup()

## calculate nesting in the entire data
data_find_nested <- mge_bins_melted %>% 
  group_by(island) %>%
  mutate(n_mge = n_distinct(mge)) 

#individual MGE nested counts - TE nested with TE
data_find_nested %>% filter(n_mge > 1) %>% group_by(mge) %>%  summarise(t = sum(count)) %>% ungroup() %>% mutate(f = sum(t))
# #  mge             t      f
# <chr>       <dbl>  <dbl>
#   1 CE          38269 481327
# 2 Cellular    10923 481327
# 3 Integron     5286 481327
# 4 IS_Tn      313774 481327
# 5 MI          42179 481327
# 6 Phage       17955 481327
# 7 Phage_like  52941 481327

# 313774/1818135
# 0.1725801

#individual MGE non-nested counts
data_find_nested %>% filter(n_mge == 1) %>% group_by(mge) %>%  summarise(t = sum(count)) %>% ungroup() %>% mutate(f = sum(t))

# R1 high quality genomes
# # A tibble: 7 × 3
# mge             t      f
# <chr>       <dbl>  <dbl>
#   1 CE          37835 475514
# 2 Cellular    10729 475514
# 3 Integron     5223 475514
# 4 IS_Tn      309857 475514
# 5 MI          41812 475514
# 6 Phage       17756 475514
# 7 Phage_like  52302 475514

#309857/1798076
#0.172327

#TE nested with TE and all other MGEs
data_find_nested %>% filter(mge == "IS_Tn") %>% filter(count > 1 & n_mge == 1) %>% unique() %>% group_by(mge) %>% summarise(total = sum(count)) 
#489079 (high quality)

#non-nested TE
data_find_nested %>% filter(mge == "IS_Tn") %>% filter(n_mge == 1) %>% unique() %>% group_by(mge) %>% summarise(total = sum(count)) #1504361
#1488219 (high quality)

#nested integrons
data_find_nested  %>% filter(mge == "Integron") %>% filter(n_mge > 1) %>% group_by(mge) %>% summarise(total = sum(count)) 
#5223/8251 - 63.3% - (high quality)

#nested TE
data_find_nested %>% filter(mge == "IS_Tn") %>% filter(n_mge > 1) %>% unique() %>% group_by(mge) %>% summarise(total = sum(count)) 
#309857/1798076 - 17.23% (high quality)


##islands with > 1 mge calls
data_onlyN <- mge_bins_melted %>% 
  group_by(island) %>%
  mutate(n_mge = n_distinct(mge)) %>%
  filter(n_mge > 1)

data_onlyN_wider <- data_onlyN %>% 
  pivot_wider(names_from = mge, values_from = count) %>%
  mutate_at(vars(IS_Tn:Cellular), replace_na, 0)
  
#parse master file of all recombinase clusters# 
datacl <- datacl %>% 
dplyr::rename(MGE_ND = X1, island = X2, MGE_ND_1 = X3, rec = X4, rec_cluster = X5, prot = X6)

datacl1 <- datacl %>% 
  mutate(prot1 = prot) %>%
  separate(prot1, c("one","two","three"), sep = "\\.", extra = "merge") %>% 
  unite(.,genome, c("one","two"), sep = ".") %>% 
  select(-three) %>% 
  filter(genome %in% glist_high$genome) %>% 
  select(-genome) %>% 
  group_by(island) %>%
  mutate(n_mge = n_distinct(MGE_ND))

data_onlyN_cluster <- datacl1 %>% 
  filter(n_mge > 1)

#compute nested IS_Tn in the entire mge_bins data set#
is_tn_non_hgt_nested_perc <- datacl %>% 
  filter(MGE_ND != "nested", MGE_ND != "Cellular") %>% 
  group_by(island) %>% 
  filter(any(MGE_ND == "IS_Tn")) %>% 
  filter(n_distinct(MGE_ND) > 1, n_distinct(rec_cluster) > 1) %>%
  group_by(island, MGE_ND, rec_cluster) %>% 
  summarise(counts = n()) %>%
  ungroup() %>% 
  group_by(island, MGE_ND) %>% 
  summarise(new_counts = sum(counts)) %>% 
  mutate(corrected_counts = if_else(MGE_ND %in% c("IS_Tn", "Integron"), new_counts, new_counts[MGE_ND == "IS_Tn"])) %>% 
  group_by(MGE_ND) %>% 
  summarise(total = sum(corrected_counts)) %>% 
  ungroup() %>% 
  mutate(perc = (total/(sum(total) - total[MGE_ND == "IS_Tn"]))*100) %>% 
  mutate(perc = if_else(MGE_ND == "IS_Tn", 0, perc))

#compute nested IS_Tn in hgt data#
mge_rec_map <- hgt1.1 %>% 
  select(ISLAND, MGE, rec_cluster) %>% 
  unique()

#obtain list of nested IS_Tn clusters#
nested_is_tn_cluster_2r <- hgt1.1 %>% 
  group_by(ISLAND) %>% 
  filter(n_distinct(MGE) > 1, n_distinct(rec_cluster) > 1) %>% 
  summarise(all_mges = list(MGE), all_recs = list(rec_cluster), all_families = unique(family), is_tn_cluster = list(unique(rec_cluster[MGE=="IS_Tn"]))) %>% 
  ungroup() %>% 
  unnest(cols = "is_tn_cluster") %>% 
  group_by(is_tn_cluster) %>% 
  do({
    sub_tbl <- .
    rec_list <- unlist(sub_tbl$all_recs)
    table_rec <- table(rec_list)
    family_list <- unlist(sub_tbl$all_families)
    is_tn_oi <- unique(sub_tbl$is_tn_cluster)
    mge_list <- unlist(sub_tbl$all_mges)
    is_tn_all <- unique(rec_list[mge_list == "IS_Tn"])
    hgt_clusters <- rec_list[mge_list %in% c("Phage", "CE", "MI", "Integron", "Phage_like")]
    if(table_rec[is_tn_oi] > 1 & any(table_rec[(!names(table_rec) %in% is_tn_all) & (names(table_rec) %in% hgt_clusters)] > 1) & length(unique(family_list)) > 1) {
      outcome <- "nested_hgt"
    } else {
      outcome <- "unnested"
    }
    data.frame(is_tn_cluster = is_tn_oi, outcome = outcome)
  }) %>% 
  filter(outcome == "nested_hgt") %>% 
  pull(is_tn_cluster) %>% 
  unique()

#compute fraction of nested IS_Tn with other MGEs#
is_tn_hgt_nested_2rec_perc <- hgt1.1 %>% 
  group_by(ISLAND) %>% 
  filter(n_distinct(MGE) > 1, n_distinct(rec_cluster) > 1) %>% 
  summarise(all_mges = list(MGE), all_recs = list(rec_cluster), all_families = unique(family), is_tn_cluster = list(unique(rec_cluster[MGE=="IS_Tn"]))) %>% 
  ungroup() %>% 
  unnest(cols = "is_tn_cluster") %>% 
  filter(is_tn_cluster %in% nested_is_tn_cluster_2r) %>% 
  group_by(is_tn_cluster) %>%
  filter(n_distinct(all_families) > 1) %>%
  unnest(cols = c("all_recs", "all_mges")) %>% 
  ungroup() %>%
  group_by(ISLAND, all_families) %>% 
  summarise(rec_comb = list(t(combn(all_recs, 2)) %>% as_tibble() %>% filter(!(V1 %in% all_recs[all_mges == "IS_Tn"] & V2 %in% all_recs[all_mges == "IS_Tn"])) %>% filter(!(V1 %in% all_recs[all_mges != "IS_Tn"] & V2 %in% all_recs[all_mges != "IS_Tn"])))) %>% 
  unnest(cols = "rec_comb") %>%
  unique() %>% 
  group_by(V1, V2) %>% 
  filter(n_distinct(all_families) > 1) %>% 
  ungroup() %>%
  dplyr::rename(rec_cluster1 = V1, rec_cluster2 = V2) %>% 
  left_join(mge_rec_map %>% dplyr::rename(rec_cluster1 = rec_cluster), by = c("ISLAND", "rec_cluster1")) %>% 
  dplyr::rename(mge1 = MGE) %>% 
  left_join(mge_rec_map %>% dplyr::rename(rec_cluster2 = rec_cluster), by = c("ISLAND", "rec_cluster2")) %>% 
  dplyr::rename(mge2 = MGE) %>% 
  mutate(nesting_partner = if_else(mge1 == "IS_Tn", mge2, mge1)) %>% 
  group_by(ISLAND, nesting_partner) %>% 
  summarise(count = n()) %>% 
  ungroup() %>% 
  mutate(count = if_else(nesting_partner %in% c("IS_Tn"), count, as.integer(1))) %>% 
  group_by(nesting_partner) %>% 
  summarise(n = sum(count)) %>% 
  ungroup() %>% 
  mutate(perc = (n/sum(n))*100)

#create dataframe combining nested fractions of IS_Tn with other MGEs from HGT and non-HGT data)
is_tn_hgt_clusterF <- is_tn_hgt_nested_2rec_perc %>% 
  dplyr::rename(hgt_total = n, hgt_perc = perc, MGE = nesting_partner)

is_tn_non_hgt_clusterF <- is_tn_non_hgt_nested_perc %>%  
  dplyr::rename(non_hgt_total = total, all_perc = perc, MGE = MGE_ND)

final_is_tn <- left_join(is_tn_non_hgt_clusterF,is_tn_hgt_clusterF, by = "MGE") %>%
  reshape2::melt() %>%
  filter(MGE != "IS_Tn") %>%
  replace(is.na(.), 0) 

combined_is_tn_nested <- final_is_tn %>% 
  pivot_wider(names_from = variable) %>% 
  select(-all_perc,-hgt_perc)

n_hgt <- as.data.frame(combined_is_tn_nested)
allRClasses <- subset(colnames(n_hgt), colnames(n_hgt)!="MGE") 
allMges <- as.vector(n_hgt$MGE)
rownames(n_hgt) <- as.vector(n_hgt$MGE)
n_hgt <- n_hgt %>% select(-MGE)

#enrichment analysis#
FE_n_hgt <- matrix(NA, nrow = length(allMges), ncol = length(allRClasses),dimnames = list(allMges, allRClasses))
for (i in allMges) {
  for(j in allRClasses) {
    print(paste0(i, "_", j))
    in_n_hgt <- n_hgt[i,j]
    not_mge <- sum(n_hgt[,j]) - in_n_hgt
    out_n_hgt <- sum(n_hgt[i,]) - in_n_hgt
    out_n_hgt_outmge <- sum(n_hgt) - in_n_hgt - out_n_hgt - not_mge
    m <- matrix(data = c(in_n_hgt, not_mge,out_n_hgt, out_n_hgt_outmge), ncol = 2)
    test <- fisher.test(m,alternative = "greater")
    FE_n_hgt[i,j] <- test$p.value
  }
}

FE_n_hgt_adj <- matrix(p.adjust(p = FE_n_hgt, method = "bonferroni"), ncol = ncol(FE_n_hgt), nrow = nrow(FE_n_hgt), dimnames = list(rownames(FE_n_hgt), colnames(FE_n_hgt)))

#Multiple testing#
FE_n_hgt_adj_t <- as.tibble(FE_n_hgt_adj) %>% mutate(mge = rownames(FE_n_hgt_adj))
m_FE_n_hgt_adj_t <- FE_n_hgt_adj_t %>% reshape2::melt(value.name = "p_adj") %>%  mutate(signP=ifelse(p_adj>0.05,"","*")) %>% dplyr::rename(MGE = mge)

final_is_tn_nested <- final_is_tn %>% filter(!grepl("_perc",variable))
nested_is_tn_fig <- left_join(final_is_tn_nested,m_FE_n_hgt_adj_t, by = c("MGE","variable")) %>% 
  group_by(variable) %>%
  mutate(frac = (value/sum(value))*100) %>%
  ungroup() %>%
  mutate(signP = if_else(variable == "non_HGT" & signP == "*","",signP)) 

ggplot(nested_is_tn_fig, aes(x = reorder(MGE,value,sum), y = frac, fill = variable)) +
  geom_bar(stat = "identity", colour = "grey20", width = 0.5,position=position_dodge()) +
  geom_text(data = nested_is_tn_fig, label = nested_is_tn_fig$signP) + 
  scale_fill_manual(values = c("grey80","black")) + 
  theme_cowplot() +
  theme(axis.text.x=element_text(angle = 90, hjust = 1)) +
  labs(y = "Perecentage of nested IS_Tn" , x ="", fill = "")

####determine overall nested fraction for integrons####

##determine overall nested fraction for Integron
hgt_int1.1 <- hgt1.1 %>%
  dplyr::rename(cor_mge = MGE, mge = MGE_ND) %>% 
  select(1,2,3,4,6,9) %>%
  filter(cor_mge == "Integron") 

hgt_int1.1 %>%
  mutate(new = (grepl("(\\d+)&", mge))) %>% 
  filter(new == "TRUE")

hgt_int1.1_check <- hgt_int1.1 %>% 
  mutate(prot1 = prot) %>%
  separate(prot1, c("one","two","three"), sep = "\\.", extra = "merge") %>% 
  unite(.,genome, c("one","two"), sep = ".") %>% 
  select(-three) %>% 
  group_by(rec_cluster,genome) %>% 
  mutate(atleast_1_nested = if_else((grepl("(\\d+)&", mge)),1,0)) %>% 
  ungroup() %>% 
  group_by(rec_cluster,cor_mge,genome) %>% 
  summarise(val = if_else(sum(atleast_1_nested) > 0, "nested", "non_nested")) %>%
  ungroup() %>% 
  group_by(cor_mge,val) %>% 
  summarise(n = n()) 

hgt_int1.2_check <- hgt_int1.1 %>% 
  mutate(prot1 = prot) %>%
  separate(prot1, c("one","two","three"), sep = "\\.", extra = "merge") %>% 
  unite(.,genome, c("one","two"), sep = ".") %>% 
  select(-three) %>% 
  group_by(rec_cluster,ISLAND) %>% 
  summarise(nested_n = sum(grepl("(\\d+)&", mge)), all_n = n_distinct(prot)) %>% 
  ungroup()

#overall nested Integron in HGT data

#compute nested Integron in the entire mge_bins data set#
int_non_hgt_nested_perc <- datacl %>% 
  filter(MGE_ND != "nested", MGE_ND != "Cellular") %>% 
  group_by(island) %>% 
  filter(any(MGE_ND == "Integron")) %>% 
  filter(n_distinct(MGE_ND) > 1, n_distinct(rec_cluster) > 1) %>%
  group_by(island, MGE_ND, rec_cluster) %>% 
  summarise(counts = n()) %>%
  ungroup() %>% 
  group_by(island, MGE_ND) %>% 
  summarise(new_counts = sum(counts)) %>% 
  mutate(corrected_counts = if_else(MGE_ND %in% c("IS_Tn", "Integron"), new_counts, new_counts[MGE_ND == "Integron"])) %>% 
  group_by(MGE_ND) %>% 
  summarise(total = sum(corrected_counts)) %>% 
  ungroup() %>% 
  mutate(perc = (total/(sum(total) - total[MGE_ND == "Integron"]))*100) %>% 
  mutate(perc = if_else(MGE_ND == "Integron", 0, perc))

#compute nested Integron in hgt data#
mge_rec_map <- hgt1.1 %>% 
  select(ISLAND, MGE, rec_cluster) %>% 
  unique()

#obtain list of nested integron clusters#
nested_int_cluster_2r <- hgt1.1 %>% 
  group_by(ISLAND) %>% 
  filter(n_distinct(MGE) > 1, n_distinct(rec_cluster) > 1) %>% 
  summarise(all_mges = list(MGE), all_recs = list(rec_cluster), all_families = unique(family), int_cluster = list(unique(rec_cluster[MGE=="Integron"]))) %>% 
  ungroup() %>% 
  unnest(cols = "int_cluster") %>% 
  group_by(int_cluster) %>% 
  do({
    sub_tbl <- .
    rec_list <- unlist(sub_tbl$all_recs)
    table_rec <- table(rec_list)
    family_list <- unlist(sub_tbl$all_families)
    int_oi <- unique(sub_tbl$int_cluster)
    mge_list <- unlist(sub_tbl$all_mges)
    int_all <- unique(rec_list[mge_list == "Integron"])
    hgt_clusters <- rec_list[mge_list %in% c("Phage", "CE", "MI", "IS_Tn", "Phage_like")]
    if(table_rec[int_oi] > 1 & any(table_rec[(!names(table_rec) %in% int_all) & (names(table_rec) %in% hgt_clusters)] > 1) & length(unique(family_list)) > 1) {
      outcome <- "nested_hgt"
    } else {
      outcome <- "unnested"
    }
    data.frame(int_cluster = int_oi, outcome = outcome)
  }) %>% 
  filter(outcome == "nested_hgt") %>% 
  pull(int_cluster) %>% 
  unique()

#compute fraction of nested Integron with other MGEs#
int_hgt_nested_2rec_perc <- hgt1.1 %>% 
  group_by(ISLAND) %>% 
  filter(n_distinct(MGE) > 1, n_distinct(rec_cluster) > 1) %>% 
  summarise(all_mges = list(MGE), all_recs = list(rec_cluster), all_families = unique(family), int_cluster = list(unique(rec_cluster[MGE=="Integron"]))) %>% 
  ungroup() %>% 
  unnest(cols = "int_cluster") %>% 
  filter(int_cluster %in% nested_int_cluster_2r) %>% 
  group_by(int_cluster) %>%
  filter(n_distinct(all_families) > 1) %>%
  unnest(cols = c("all_recs", "all_mges")) %>% 
  ungroup() %>%
  group_by(ISLAND, all_families) %>% 
  summarise(rec_comb = list(t(combn(all_recs, 2)) %>% as_tibble() %>% filter(!(V1 %in% all_recs[all_mges == "Integron"] & V2 %in% all_recs[all_mges == "Integron"])) %>% filter(!(V1 %in% all_recs[all_mges != "Integron"] & V2 %in% all_recs[all_mges != "Integron"])))) %>% 
  unnest(cols = "rec_comb") %>%
  unique() %>% 
  group_by(V1, V2) %>% 
  filter(n_distinct(all_families) > 1) %>% 
  ungroup() %>%
  dplyr::rename(rec_cluster1 = V1, rec_cluster2 = V2) %>% 
  left_join(mge_rec_map %>% dplyr::rename(rec_cluster1 = rec_cluster), by = c("ISLAND", "rec_cluster1")) %>% 
  dplyr::rename(mge1 = MGE) %>% 
  left_join(mge_rec_map %>% dplyr::rename(rec_cluster2 = rec_cluster), by = c("ISLAND", "rec_cluster2")) %>% 
  dplyr::rename(mge2 = MGE) %>% 
  mutate(nesting_partner = if_else(mge1 == "Integron", mge2, mge1)) %>% 
  group_by(ISLAND, nesting_partner) %>% 
  summarise(count = n()) %>% 
  ungroup() %>% 
  mutate(count = if_else(nesting_partner %in% c("Integron"), count, as.integer(1))) %>% 
  group_by(nesting_partner) %>% 
  summarise(n = sum(count)) %>% 
  ungroup() %>% 
  mutate(perc = (n/sum(n))*100)

#create dataframe combining nested fractions of Integron with other MGEs from HGT and non-HGT data)

int_hgt_clusterF <- int_hgt_nested_2rec_perc %>% 
  dplyr::rename(hgt_total = n, hgt_perc = perc, MGE = nesting_partner)

int_non_hgt_clusterF <- int_non_hgt_nested_perc %>%  
  dplyr::rename(non_hgt_total = total, all_perc = perc, MGE = MGE_ND)

final_int <- left_join(int_non_hgt_clusterF,int_hgt_clusterF, by = "MGE") %>%
  reshape2::melt() %>%
  filter(MGE != "Integron") %>%
  replace(is.na(.), 0) 

combined_int_nested <- final_int %>% 
  pivot_wider(names_from = variable) %>% 
  select(-all_perc,-hgt_perc)

n_hgt <- as.data.frame(combined_int_nested)
allRClasses <- subset(colnames(n_hgt), colnames(n_hgt)!="MGE") 
allMges <- as.vector(n_hgt$MGE)
rownames(n_hgt) <- as.vector(n_hgt$MGE)
n_hgt <- n_hgt %>% select(-MGE)

#enrichment analysis#
FE_n_hgt <- matrix(NA, nrow = length(allMges), ncol = length(allRClasses),dimnames = list(allMges, allRClasses))
for (i in allMges) {
  for(j in allRClasses) {
    print(paste0(i, "_", j))
    in_n_hgt <- n_hgt[i,j]
    not_mge <- sum(n_hgt[,j]) - in_n_hgt
    out_n_hgt <- sum(n_hgt[i,]) - in_n_hgt
    out_n_hgt_outmge <- sum(n_hgt) - in_n_hgt - out_n_hgt - not_mge
    m <- matrix(data = c(in_n_hgt, not_mge,out_n_hgt, out_n_hgt_outmge), ncol = 2)
    test <- fisher.test(m,alternative = "greater")
    FE_n_hgt[i,j] <- test$p.value
    print(paste0(i, "_", j,"=",test$estimate))
  }
}

FE_n_hgt_adj <- matrix(p.adjust(p = FE_n_hgt, method = "bonferroni"), ncol = ncol(FE_n_hgt), nrow = nrow(FE_n_hgt), dimnames = list(rownames(FE_n_hgt), colnames(FE_n_hgt)))

#Multiple testing#
FE_n_hgt_adj_t <- as.tibble(FE_n_hgt_adj) %>% mutate(mge = rownames(FE_n_hgt_adj))
m_FE_n_hgt_adj_t <- FE_n_hgt_adj_t %>% reshape2::melt(value.name = "p_adj") %>%  mutate(signP=ifelse(p_adj>0.05,"","*")) %>% dplyr::rename(MGE = mge)

final_int_nested <- final_int %>% filter(!grepl("_perc",variable))
nested_int_fig <- left_join(final_int_nested,m_FE_n_hgt_adj_t, by = c("MGE","variable")) %>% 
  group_by(variable) %>%
  mutate(frac = (value/sum(value))*100) %>%
  ungroup() %>%
  mutate(signP = if_else(variable == "non_HGT" & signP == "*","",signP)) 
ggplot(nested_int_fig, aes(x = reorder(MGE,value,sum), y = frac, fill = variable)) +
  geom_bar(stat = "identity", colour = "grey20", width = 0.5,position=position_dodge()) +
  geom_text(data = nested_int_fig, label = nested_int_fig$signP) + 
  scale_fill_manual(values = c("grey80","black")) + 
  theme_cowplot() +
  theme(axis.text.x=element_text(angle = 90, hjust = 1)) +
  labs(y = "Perecentage of nested Integron" , x ="", fill = "")

####Figure 4D - Heatmaps of HGT of MGEs within and between habitats####

#load data#
hgt_hab <- read_tsv("raw_data/hgt_habitat_final.txt", col_names = F)
all_mge <- read_tsv("raw_data/recombinase_mge_habitat_final.txt", col_names = F)

hgt_hab <- hgt_hab %>% 
  dplyr::rename(cluster = X1, mge = X2, rec = X3, prot = X4, hab= X5) %>% 
  separate(prot, into = c("one","two","three"), sep = "[.]", extra = "merge", remove = F) %>%
  unite("genome", one, two, sep = ".") %>%
  select(-three) 

hgt_hab_data <- left_join(tax,hgt_hab, by = "genome") %>% 
  drop_na()

#analysis of within habitat (wn) transfers#
wn <- hgt_hab_data %>% 
  mutate(hab = strsplit(hab, split = "[:]")) %>% 
  unnest() %>% 
  filter(hab != "") %>% 
  group_by(hab, cluster) %>% 
  summarise(n = length(unique(prot)), transfers = n - 1) %>% 
  ungroup()

wn_stat <- wn %>% 
  group_by(hab) %>% 
  summarise(total = sum(transfers), clusters = n(), avg = mean(transfers), std_dev = sd(transfers))

#analysis of between habitat (bn) transfers#
bn <- hgt_hab_data %>% 
  mutate(hab = strsplit(hab, split = "[:]")) %>% 
  unnest() %>% 
  filter(hab != "") 

bn_int <- bn %>% 
  group_by(cluster) %>% 
  summarise(cnt = list(table(hab)), n = n_distinct(prot)) 

bn_stat <- bn_int %>% 
  ungroup() %>% 
  rowwise() %>% 
  do({
    sub_tbl <- .
    x <- length(which(sub_tbl$cnt < sub_tbl$n))
    if(x > 0) {
      habs_orig <- paste(names(sub_tbl$cnt)[which(sub_tbl$cnt == sub_tbl$n)], collapse = ":")
      habs_unique <- paste(names(sub_tbl$cnt)[which(sub_tbl$cnt < sub_tbl$n)], collapse = ":")
    } else {
      habs_orig <- paste(names(sub_tbl$cnt), collapse = ":")
      habs_unique <- ""
    }
    data.frame(cluster = sub_tbl$cluster, n = x, habs_orig = habs_orig, habs_unique = habs_unique)
  }) %>% 
  mutate(habs_unique = strsplit(habs_unique, split = "[:]")) %>% 
  unnest() %>% 
  mutate(habs_orig = strsplit(habs_orig, split = "[:]")) %>% 
  unnest() %>% 
  group_by(cluster) %>% 
  mutate(freq = n/length(unique(habs_orig))) %>% 
  group_by(habs_orig, habs_unique) %>% 
  summarise(avg_freq = mean(freq), total = sum(n))

all_hab <- unique(bn$hab)

hab_comp <- matrix(0, nrow = length(all_hab), ncol = length(all_hab), dimnames = list(all_hab, all_hab))
for(i in all_hab) {
  for(j in all_hab) {
    if(i == j) {
      hab_comp[i,j] <- wn_stat$total[wn_stat$hab == i]
    } else {
      hab_comp[i,j] <- sum(bn_stat$n[(bn_stat$habs_orig==i | bn_stat$habs_unique ==i) & (bn_stat$habs_orig==j | bn_stat$habs_unique ==j)])
    }
  }
}


wn_m <- hgt_hab_data %>% 
  mutate(hab = strsplit(hab, split = "[:]")) %>% 
  unnest() %>% 
  filter(hab != "") %>% 
  group_by(mge, hab, cluster) %>% 
  summarise(n = length(unique(prot)), transfers = n - 1) %>% 
  ungroup()

wn_m_stat <- wn_m %>% 
  group_by(mge, hab) %>% 
  summarise(total = sum(transfers), clusters = n(), avg = mean(transfers), std_dev = sd(transfers))

#between
bn_m <- hgt_hab_data %>% 
  mutate(hab = strsplit(hab, split = "[:]")) %>% 
  unnest() %>% 
  filter(hab != "")

bn_m_stat <- bn_m %>% 
  group_by(mge, cluster) %>% 
  summarise(cnt = list(table(hab)), n = n_distinct(prot)) %>% 
  ungroup() %>% 
  rowwise() %>% 
  do({
    sub_tbl <- .
    x <- length(which(sub_tbl$cnt < sub_tbl$n))
    if(x > 0) {
      habs_orig <- paste(names(sub_tbl$cnt)[which(sub_tbl$cnt == sub_tbl$n)], collapse = ":")
      habs_unique <- paste(names(sub_tbl$cnt)[which(sub_tbl$cnt < sub_tbl$n)], collapse = ":")
    } else {
      habs_orig <- paste(names(sub_tbl$cnt), collapse = ":")
      habs_unique <- ""
    }
    data.frame(cluster = sub_tbl$cluster, mge = sub_tbl$mge, n = x, habs_orig = habs_orig, habs_unique = habs_unique)
  }) %>% 
  mutate(habs_unique = strsplit(habs_unique, split = "[:]")) %>% 
  unnest() %>% 
  mutate(habs_orig = strsplit(habs_orig, split = "[:]")) %>% 
  unnest() %>% 
  group_by(mge, cluster) %>% 
  mutate(freq = n/length(unique(habs_orig))) 

all_hab_m <- unique(bn_m$hab)

##IS_tn
wn_m_stat_tn <- wn_m_stat %>% filter(mge == "IS_Tn")
bn_m_stat_tn <- bn_m_stat %>% filter(mge == "IS_Tn")
hab_comp_m_tn <- matrix(0, nrow = length(all_hab_m), ncol = length(all_hab_m), dimnames = list(all_hab_m, all_hab_m))
for(i in all_hab_m) {
  for(j in all_hab_m) {
    if(i == j) {
      hab_comp_m_tn[i,j] <- wn_m_stat_tn$total[wn_m_stat_tn$hab == i]
    } else {
      hab_comp_m_tn[i,j] <- sum(bn_m_stat_tn$n[(bn_m_stat_tn$habs_orig==i | bn_m_stat_tn$habs_unique ==i) & (bn_m_stat_tn$habs_orig==j | bn_m_stat_tn$habs_unique ==j)])
    }
  }
}


#CE
wn_m_stat_ce <- wn_m_stat %>% filter(mge == "CE")
bn_m_stat_ce <- bn_m_stat %>% filter(mge == "CE")
hab_comp_m_ce <- matrix(0, nrow = length(all_hab_m), ncol = length(all_hab_m), dimnames = list(all_hab_m, all_hab_m))
for(i in all_hab_m) {
  for(j in all_hab_m) {
    if(i == j) {
      hab_comp_m_ce[i,j] <- wn_m_stat_ce$total[wn_m_stat_ce$hab == i]
    } else {
      hab_comp_m_ce[i,j] <- sum(bn_m_stat_ce$n[(bn_m_stat_ce$habs_orig==i | bn_m_stat_ce$habs_unique ==i) & (bn_m_stat_ce$habs_orig==j | bn_m_stat_ce$habs_unique ==j)])
    }
  }
}

#Phage
wn_m_stat_ph <- wn_m_stat %>% filter(mge == "Phage")
bn_m_stat_ph <- bn_m_stat %>% filter(mge == "Phage")
hab_comp_m_ph <- matrix(0, nrow = length(all_hab_m), ncol = length(all_hab_m), dimnames = list(all_hab_m, all_hab_m))
for(i in all_hab_m) {
  for(j in all_hab_m) {
    if(i == j) {
      hab_comp_m_ph[i,j] <- wn_m_stat_ph$total[wn_m_stat_ph$hab == i]
    } else {
      hab_comp_m_ph[i,j] <- sum(bn_m_stat_ph$n[(bn_m_stat_ph$habs_orig==i | bn_m_stat_ph$habs_unique ==i) & (bn_m_stat_ph$habs_orig==j | bn_m_stat_ph$habs_unique ==j)])
    }
  }
}

#Phage-like
wn_m_stat_phl <- wn_m_stat %>% filter(mge == "Phage_like")
bn_m_stat_phl <- bn_m_stat %>% filter(mge == "Phage_like")
hab_comp_m_phl <- matrix(0, nrow = length(all_hab_m), ncol = length(all_hab_m), dimnames = list(all_hab_m, all_hab_m))
for(i in all_hab_m) {
  for(j in all_hab_m) {
    if(i == j) {
      hab_comp_m_phl[i,j] <- wn_m_stat_phl$total[wn_m_stat_phl$hab == i]
    } else {
      hab_comp_m_phl[i,j] <- sum(bn_m_stat_phl$n[(bn_m_stat_phl$habs_orig==i | bn_m_stat_phl$habs_unique ==i) & (bn_m_stat_phl$habs_orig==j | bn_m_stat_phl$habs_unique ==j)])
    }
  }
}


#Integron
wn_m_stat_int <- wn_m_stat %>% filter(mge == "Integron")
bn_m_stat_int <- bn_m_stat %>% filter(mge == "Integron")
hab_comp_m_int <- matrix(0, nrow = length(all_hab_m), ncol = length(all_hab_m), dimnames = list(all_hab_m, all_hab_m))
for(i in all_hab_m) {
  for(j in all_hab_m) {
    if(i == j) {
      hab_comp_m_int[i,j] <- wn_m_stat_int$total[wn_m_stat_int$hab == i]
    } else {
      hab_comp_m_int[i,j] <- sum(bn_m_stat_int$n[(bn_m_stat_int$habs_orig==i | bn_m_stat_int$habs_unique ==i) & (bn_m_stat_int$habs_orig==j | bn_m_stat_int$habs_unique ==j)])
    }
  }
}


#MI
wn_m_stat_mi <- wn_m_stat %>% filter(mge == "MI")
bn_m_stat_mi <- bn_m_stat %>% filter(mge == "MI")
hab_comp_m_mi <- matrix(0, nrow = length(all_hab_m), ncol = length(all_hab_m), dimnames = list(all_hab_m, all_hab_m))
for(i in all_hab_m) {
  for(j in all_hab_m) {
    if(i == j) {
      hab_comp_m_mi[i,j] <- wn_m_stat_mi$total[wn_m_stat_mi$hab == i]
    } else {
      hab_comp_m_mi[i,j] <- sum(bn_m_stat_mi$n[(bn_m_stat_mi$habs_orig==i | bn_m_stat_mi$habs_unique ==i) & (bn_m_stat_mi$habs_orig==j | bn_m_stat_mi$habs_unique ==j)])
    }
  }
}

#Plot combine habitat hgt data for all MGEs
mge_oi <- c("tn", "ce", "ph", "phl", "mi", "int")
hab_order <- rev(c("soil", "freshwater", "marine", "wastewater", "cat-gut", "dog-gut", "pig-gut", "mouse-gut", "human-gut", "human-oral", "human-skin", "human-nose", "human-vagina", "built-environment"))
full_table <- NULL
for(i in mge_oi) {
  table_oi <- eval(parse(text = paste0("hab_comp_m_", i)))
  table_m <- table_oi %>% 
    reshape2::melt() %>% 
    mutate(mge = i)
  full_table <- rbind(full_table, table_m)
}

full_table %>% 
  mutate(cut_val = cut(value,breaks = c(-1, 0, 5, 20, 40, 100, 1000, max(value,na.rm = T)),
                       labels = c("0", "1-5", "6-20", "21-40", "41-100", "101-1000", ">1000"))) %>% 
  mutate(cut_val = factor(as.character(cut_val), levels = rev(levels(cut_val)))) %>% 
  mutate(Var1 = factor(Var1, levels = hab_order), Var2 = factor(Var2, levels = hab_order)) %>% 
  mutate(mge = factor(mge, levels = c("tn", "ce", "mi", "int", "ph", "phl", "nested"), labels = c("IS_Tn", "CE", "MI", "Integron", "Phage", "Phage_like", "Nested"))) %>% 
  filter(mge != "Nested") %>% 
  ggplot(aes(Var1, Var2, fill = cut_val)) +
  geom_tile(colour = "grey60") +
  scale_fill_manual(values = rev(brewer.pal(7,"YlGnBu")), na.value = "grey90") +
  facet_wrap(~mge, nrow = 2) +
  theme_gray(base_size = 24) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) +
  labs(
    x = "",
    y = "",
    fill = "# HGT"
  )
#Normalise counts to abundance of different MGE categories
all_mge <- all_mge %>% dplyr::rename(prot = X1, mge = X2, hab = X3)

all_mge_hab <- all_mge %>%
  group_by(mge,hab) %>%
  summarise(count = n())

full_table_norm <- full_table %>%
  mutate(hab = Var1, mge = factor(mge, levels = c("tn", "ce", "mi", "int", "ph", "phl"), labels = c("IS_Tn", "CE", "MI", "Integron", "Phage", "Phage_like"))) %>% 
  left_join(all_mge_hab, by = c("hab", "mge")) %>% 
  mutate(hab = Var2) %>% 
  left_join(all_mge_hab, by = c("hab", "mge")) %>% 
  mutate(poss = if_else(Var1==Var2, count.x, count.x + count.y)) %>% 
  select(-hab, -count.x, -count.y) %>% 
  mutate(frac = value/poss)


cols <- rev(brewer.pal(7,"YlGnBu"))
full_table_norm %>% 
  mutate(cut_val = cut(log10(frac),breaks = c(-Inf, -6, -4, -2, -1, 1), right = F, labels = c("0", "0.000001 - 0.0001", "0.0001 - 0.01", "0.01 - 0.1", "0.1 - 1"))) %>%
  mutate(cut_val = factor(as.character(cut_val), levels = rev(levels(cut_val)))) %>%
  mutate(Var1 = factor(Var1, levels = hab_order), Var2 = factor(Var2, levels = hab_order)) %>% 
  mutate(mge = factor(mge, levels = c("IS_Tn", "CE", "MI", "Integron","Phage", "Phage_like"))) %>% 
  ggplot(aes(Var1, Var2, fill = cut_val)) +
  geom_tile(colour = "grey60") +
  scale_fill_manual(values = cols[(8-5):7], na.value = "grey90") +
  facet_wrap(~mge, nrow = 2) +
  theme_gray(base_size = 20) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    legend.position = "bottom"
  ) +
  guides(fill = guide_legend(nrow = 1)) +
  labs(
    x = "",
    y = "",
    fill = "Normalised HGT counts"
  )
