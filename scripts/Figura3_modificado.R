# Título del manuscrito: Paisaje de elementos genéticos móviles y su carga de resistencia
# a antibióticos en genomas procariotas

# El siguiente código se puede utilizar para crear la Figura 3 en el manuscrito

# Fecha de Creacion: 19-02-2021
# Autor: supriya.khedkar@embl.de

# Modificado por: Evelia Coss
# Fechad de modificacion: 6 de abril 2025

#---- Cargar paquetes ----- 
library(tidyverse)
library(cowplot)
library(phytools)
library(viridisLite)
library(viridis)
library(scales)

#---- PASO 1: Importar datos ----- 
tax <- read_tsv("data/raw_data/species_with_atleast_2genomes.list.gz", col_names=F)
db <-read_tsv("data/processed_data/mge_bins_per_genome_final.txt.gz", col_names = T)
gs <-read_tsv("data/raw_data/genome_size.txt.gz", col_names = T)
class_tree <- read.tree("data/raw_data/progenomes2_class_tree.nwk")
glist <- read_tsv("data/raw_data/genome_status_supplementary_tableS2.txt.gz", col_names = T)

#---- PASO 2: Manipulación y Limpieza de los datos  ----- 
# Obtener los genomas con la mas alta calidad
glist_high <- glist %>% 
  filter(genome_quality == "high")

# Filtra los genomas de interés y calcula el número de especies y genomas distintos
tax %>% 
  # Filtra los genomas presentes en la lista de genomas de alta calidad (glist_high$genome)
  filter(X2 %in% glist_high$genome) %>% 
  # Extrae la columna X1 (especies)
  pull(X1) %>% 
  # Calcula el número de especies distintas
  n_distinct()
# Esto devuelve: 76,228 genomas y 3207 especies.

#---- Figura 3B:----- 

# Renombra las columnas del dataframe tax para un nombre más legible
colnames(tax) <- c("specI", "genomeID", "kingdom", "phylum", "class", "genus")

# Reestructura la base de datos 'db' de formato ancho a largo usando reshape2::melt
mdb <- db %>% reshape2::melt() 

# Renombra las columnas de la base de datos reestructurada para claridad
mdb <- mdb %>% 
  dplyr::rename(mge = variable, count = value, genomeID = 'Genome') %>% 
  # Filtra para incluir solo los genomas de alta calidad (genomas en glist_high$genome)
  filter(genomeID %in% glist_high$genome)

##Format data 
#create a data.frame containing for each specI for each MGE: genome counts (genomeCnt), 
#total counts (cnt_tot), 
#average counts across genomes (avg_cnt), 
#number of genomes where the MGE was present (pa), 
#the fraction of genomes with the MGE present (frac) plus taxonomic information

# Selecciona las columnas necesarias para trabajar con los datos
mdb_cnt <- mdb %>% select(1, 7, 8)
# Obtiene los valores únicos de MGE
all_mge <- unique(mdb$mge)

# Inicializa un data.frame vacío donde se almacenarán las combinaciones
db_tax <- NULL       

# Itera sobre todos los valores únicos de MGE
for(i in 1:length(all_mge)){
  # Para cada MGE, se agrega una nueva columna 'mge' con el valor del MGE actual
  db_tax_sing <- tax %>% add_column(mge = all_mge[i])
  # Une el data.frame actual con el anterior (rbind), creando una tabla con todas las combinaciones
  db_tax <- rbind(db_tax, db_tax_sing)
}

# Muestra las primeras filas del data.frame resultante
head(db_tax)

# Filtra los genomas que están en 'glist_high$genome' y los combina con los recuentos de MGE
db_cnt_all <- db_tax %>% 
  filter(genomeID %in% glist_high$genome) %>%   # Filtra para obtener solo los genomas que están en 'glist_high$genome'
  left_join(., mdb_cnt, by = c("genomeID", "mge"))  # Realiza una unión de izquierda con 'mdb_cnt' usando 'genomeID' y 'mge'

# Reemplaza todos los valores NA en el data.frame por 0
db_cnt_all[is.na(db_cnt_all)] <- 0

# Contar el número de genomas por especie y clase taxonómica
db_cnt_all %>% 
  select(specI, genomeID, class) %>%        # Selecciona las columnas relevantes: 'specI', 'genomeID', y 'class'
  unique() %>%                              # Elimina duplicados, dejando solo combinaciones únicas de especie, genoma y clase
  group_by(class) %>%                       # Agrupa los datos por clase taxonómica ('class')
  summarise(count = n()) %>%                # Cuenta el número de genomas por cada clase
  filter(count > 9) %>%                     # Filtra para solo mostrar las clases que tienen más de 9 genomas
  View()                                    # Muestra el resultado en una ventana de visualización en RStudio

# Crear un data.frame de presencia y ausencia
## Agregar información de presencia y ausencia
db_pa_all <- db_cnt_all %>% mutate(presAbs = ifelse(count > 0, 1, 0))

## Agrupar todos los genomas por especie (specI) 
## Resumir: avg_cnt = promedio de recuento por MGE por especie (specI), frac = fracción de genomas por especie con cada MGE
db_specI <- db_pa_all %>% 
  group_by(specI, mge) %>%                  # Agrupar los datos por especie (specI) y MGE
  summarise(genomeCnt = n(),                # Contar el número de genomas por combinación de 'specI' y 'mge'
            avg_cnt = mean(count),          # Calcular el promedio de recuento por MGE por especie
            pa = sum(presAbs),              # Calcular la presencia (1) total para cada MGE por especie
            cnt_tot = sum(count)) %>%       # Calcular el total de recuentos para cada MGE por especie
  mutate(frac = pa / genomeCnt) %>%         # Calcular la fracción de genomas por especie que tienen el MGE
  left_join(., tax, by = "specI") %>%       # Unir los datos con la información taxonómica de 'tax' utilizando 'specI'
  select(-genomeID) %>%                     # Eliminar la columna 'genomeID' que no es necesaria
  unique(.)                                 # Eliminar filas duplicadas



## Resumir: avg_cnt = promedio de recuento por MGE por genoma, frac = fracción de genomas con cada MGE
db_genome <- db_pa_all %>% 
  group_by(genomeID, mge) %>%                # Agrupar los datos por genoma (genomeID) y MGE
  summarise(genomeCnt = n(),                 # Contar el número de combinaciones de genoma y MGE
            avg_cnt = mean(count),           # Calcular el promedio de recuento de cada MGE por genoma
            pa = sum(presAbs),               # Calcular la presencia total de MGE en cada genoma
            cnt_tot = sum(count)) %>%        # Calcular el total de recuentos de MGE por genoma
  mutate(frac = pa / genomeCnt) %>%          # Calcular la fracción de genomas con el MGE
  left_join(., tax, by = "genomeID") %>%     # Unir los datos con la información taxonómica utilizando 'genomeID'
  unique(.)                                  # Eliminar filas duplicadas


## Crear data frames separados para los recuentos de MGE y las fracciones de MGE
# Promedio de recuentos de cada MGE por especie (specI), donde los MGEs son las columnas y las especies (specI) las filas
all_cnt_specI <- db_specI %>% 
  select(specI, mge, avg_cnt)   # Selecciona las columnas 'specI' (especie), 'mge' (MGE), y 'avg_cnt' (promedio de recuentos)

all_frac_specI <- db_specI %>% 
  select(specI, mge, frac)      # Selecciona las columnas 'specI' (especie), 'mge' (MGE), y 'frac' (fracción de genomas con el MGE)

all_cnt_genome <- db_genome %>% 
  select(genomeID, mge, avg_cnt) # Selecciona las columnas 'genomeID' (ID de genoma), 'mge' (MGE), y 'avg_cnt' (promedio de recuentos)

## Análisis por especie (specI) a nivel de clase
# Prueba de Mann-Whitney
# Se utiliza la prueba de Mann-Whitney para analizar asociaciones de los MGEs a nivel de clase  
# El análisis se basa en la fracción de genomas por especiación (SpecI) que contiene el MGE para evitar sesgos de muestreo

# entrada de datos
all_avg_cnt_tax_class <-  db_specI %>% 
  select(specI, mge, avg_cnt, class) %>%   # Selecciona las columnas: 'specI', 'mge', 'avg_cnt' y 'class'
  filter(., !grepl("Hotspot", mge)) %>%    # Filtra para excluir los MGEs que contienen "Hotspot" en su nombre
  filter(., !grepl("Cellular", mge))       # Filtra para excluir los MGEs que contienen "Cellular" en su nombre

allMges <- unique(all_avg_cnt_tax_class$mge)  # Extrae los valores únicos de MGEs

# Filtrar para clases con al menos 10 genomas
allClasses <- all_avg_cnt_tax_class %>% 
  filter(mge == all_avg_cnt_tax_class$mge[1]) %>%  # Filtra para mantener solo el primer MGE
  group_by(class) %>%   # Agrupa los datos por clase
  summarise(cntClass = n()) %>%  # Cuenta cuántos genomas existen por clase
  filter(cntClass > 9)   # Filtra las clases con más de 9 genomas

# Selección final de datos para el análisis
all_avg_cnt_tax_class_sel_pre <- all_avg_cnt_tax_class %>% 
  filter(class %in% allClasses$class)   # Filtra las especies para incluir solo aquellas pertenecientes a clases con al menos 10 genomas


### Normalización por tamaño de genoma ###
# Se calcula el tamaño promedio del genoma por cada especie (SpecI) y se normaliza el recuento de MGEs por el tamaño del genoma

# Promedio de tamaño de genoma por especie (SpecI)
gs_int <- gs %>% 
  group_by(SpecI_id_v3) %>% 
  summarise(avg_gs = mean(ProteinGeneCounts)) %>% 
  dplyr::rename(specI = SpecI_id_v3)

# Unir la información de tamaño de genoma con el promedio de recuentos de MGEs por especie y normalizar los recuentos por tamaño de genoma
all_avg_cnt_tax_class_sel <- left_join(all_avg_cnt_tax_class_sel_pre, gs_int, by = "specI") %>% 
  mutate(norm_count = avg_cnt / avg_gs) %>%  # Normalización: se divide el recuento promedio por el tamaño del genoma
  dplyr::rename(count = avg_cnt, avg_cnt = norm_count)  # Renombrar las columnas para mayor claridad

# Realizar la prueba de Mann-Whitney para comparar recuentos de MGE entre clases
MW_all <- NULL
for(i in 1:length(allMges)){  # Iterar sobre todos los MGEs
  mgeX <- allMges[i]  # Seleccionar el MGE actual
  mgeDat <- all_avg_cnt_tax_class_sel %>% filter(mge == mgeX)  # Filtrar los datos para el MGE seleccionado
  
  MW <- sapply(seq_along(allClasses$class), function(j){
    class_sel <- allClasses$class[j]  # Seleccionar la clase actual
    ingroup <- mgeDat %>% filter(class == class_sel)  # Filtrar datos de la clase seleccionada (ingroup)
    outgroup <- mgeDat %>% filter(class != class_sel)  # Filtrar datos de las clases no seleccionadas (outgroup)
    
    # Realizar la prueba de Mann-Whitney para comparar los recuentos promedio de MGE entre las clases
    MWout <- wilcox.test(ingroup$avg_cnt, outgroup$avg_cnt, alternative = "greater")
    out <- c(class_sel, MWout$p.value)  # Guardar el resultado de la prueba (clase y p-value)
    return(out)
  })
  
  MWOut <- as.data.frame(t(MW[2,]))  # Convertir los resultados de la prueba en un data frame
  colnames(MWOut) <- MW[1,]  # Asignar los nombres de las clases como nombres de columna
  rownames(MWOut) <- mgeX  # Asignar el nombre del MGE como nombre de fila
  MW_all <- rbind(MW_all, MWOut)  # Unir los resultados de la prueba para todos los MGEs
}

# Transponer los resultados finales de la prueba
MW_all_mod <- as.data.frame(t(MW_all))


## Ajuste de valores p para múltiples comparaciones
MW_all_adj <- NULL  # Inicializar un objeto vacío para almacenar los resultados ajustados

# Iterar sobre cada columna de los resultados de la prueba de Mann-Whitney
for(k in 1:ncol(MW_all_mod)){
  MW_all_mod2 <- as.matrix(MW_all_mod)  # Convertir el data frame en una matriz para facilitar el acceso a las columnas
  pVec <- c(MW_all_mod2[,k])  # Extraer la columna k, que contiene los valores p de la prueba de Mann-Whitney
  adjP <- p.adjust(pVec, "BH")  # Ajustar los valores p utilizando el método de Benjamini-Hochberg ("BH")
  MW_all_adj <- cbind(MW_all_adj, adjP)  # Unir los valores p ajustados a la matriz final
}

# Asignar los nombres de las columnas originales a los resultados ajustados
colnames(MW_all_adj) <- colnames(MW_all_mod)

# Convertir el resultado final en un data frame
MW_all_fin <- as.data.frame(MW_all_adj)

## Crear una tabla binaria con valores p < 0.1 --> 1, de lo contrario 0
MW_all_stat <- NULL  # Inicializar un objeto vacío para almacenar los resultados binarios

# Iterar sobre cada columna de los resultados ajustados (MW_all_fin)
for(l in 1:ncol(MW_all_fin)){
  # Crear una columna "signP" con valor 1 si el valor p es menor a 0.1, y 0 en caso contrario
  MW_test <- MW_all_fin %>% mutate(signP=ifelse(MW_all_fin[,l]>0.1,0,1))  
  MW_all_stat <- cbind(MW_all_stat, MW_test$signP)  # Unir los resultados binarios
}

# Convertir la matriz a un data frame
MW_all_stat <- as.data.frame(MW_all_stat)

# Asignar los nombres de las columnas y filas a la tabla binaria
colnames(MW_all_stat) <- colnames(MW_all_fin)
rownames(MW_all_stat) <- rownames(MW_all_fin)

# Sumar los valores de "signP" por fila para obtener el número total de pruebas significativas por MGE
MW_all_stat_sum <- MW_all_stat %>% mutate(sumSign=rowSums(.))

# Agregar los nombres de las filas a la tabla para obtener un formato adecuado
MW_all_stat2 <- cbind(rownames(MW_all_stat), MW_all_stat)
colnames(MW_all_stat2)[1] <- "class"  # Renombrar la primera columna como "class"

# Crear una versión de los resultados ajustados con nombres de clase y valores p ajustados
MW_all_fin2


## Crear un mapa de calor con los promedios de cuentas por especie y clase (con al menos 10 genomas)
specI_class_hmap_w_g9 <- ggplot(plot_table_w_g9, aes(y = reorder(variable, avg_class, sum), x = reorder(class, avg_class, sum), fill = avg_class)) + 
  geom_tile(colour = "black") +  # Añadir los azulejos con bordes negros
  geom_text(aes(label = signP), colour = "white") +  # Añadir etiquetas con el valor de significancia (asterisco) en color blanco
  scale_fill_gradientn(colours = alpha(viridis(10), 0.75), values = rescale(c(0, 10, 100, 500)), na.value = "grey80")  # Escala de colores para el mapa de calor

# Aplicar un tema personalizado (theme_cowplot) y ajustar la orientación de las etiquetas del eje x
specI_class_hmap_w_g9 + 
  theme_cowplot() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
  labs(x = "", y = "", fill = "Average counts per species")  # Etiquetas de los ejes y leyenda del mapa de calor

## Crear el mapa de calor filogenético con normalización de los valores
normalize <- function(x) {  # Función para normalizar los valores entre 0 y 1
  x/max(x)
}

# Crear la matriz de datos para el mapa de calor filogenético
phylo_heatmap_mat <- all_avg_cnt_tax_class_sel1 %>% 
  filter(!grepl("NA ", class)) %>%  # Filtrar clases que no contienen "NA"
  spread(variable, avg_class) %>%  # Convertir los datos en formato ancho (cada variable como columna)
  as_tibble(column_to_rownames(var = "class")) %>%  # Convertir en tibble y establecer las clases como nombres de fila
  select(class, IS_Tn, Phage, Phage_like, CE, Integron, MI) %>%  # Seleccionar las variables de interés
  column_to_rownames(var = "class")  # Establecer las clases como nombres de fila

# Crear el árbol filogenético manteniendo solo las especies que están en la matriz de datos del mapa de calor
class_tree_w_g9 <- keep.tip(class_tree, tip = rownames(phylo_heatmap_mat))

## Ajuste de colores para el mapa de calor filogenético para datos sesgados
small_value <- unique(sort(unlist(phylo_heatmap_mat)))[2]  # Encuentra el segundo valor más pequeño en la matriz de datos para establecer un umbral
small_value <- small_value - small_value %% 0.001  # Redondea el valor a 3 decimales para crear un valor pequeño ajustado

# Histograma para visualizar la distribución de los datos en la matriz
hist(unlist(phylo_heatmap_mat), nclass = 50)

# Definir los puntos de corte para el mapa de calor según la distribución de los datos
heatmap_breaks <- c(0, small_value, 1, 2, 3, 4, seq(5, 55, 10))

# Histograma con los puntos de corte definidos para el mapa de calor
hist(unlist(phylo_heatmap_mat), breaks = heatmap_breaks, freq = T)

# Definir una paleta de colores basada en viridis, con un sesgo hacia valores más bajos
plot_colors <- c("white", colorRampPalette((viridis(10)), bias = 5)(length(heatmap_breaks) - 2))

# Crear una matriz de significancia (signP) para etiquetar el mapa de calor
sig_mat <- plot_table_w_g9 %>% 
  select(class, variable, signP)  # Selecciona las columnas relevantes de los datos

# Filtra las filas que contienen "NA" en la columna "class" y transforma los datos a formato ancho
sig_mat <- sig_mat %>% 
  filter(!grepl("NA ", class)) %>%
  spread(variable, signP) %>%  # Transforma los datos para que las variables sean columnas
  column_to_rownames(var = "class")  # Establece la columna "class" como nombres de fila

# Alinea la matriz de significancia con el árbol filogenético
sig_mat <- sig_mat[class_tree_w_g9$tip.label, colnames(phylo_heatmap_mat)]

# ---------- phylo.heatmap.coords() --------------------
#' Extract X and Y coordinates from phylogenetic heatmap
#'
#' This function computes the x and y coordinates for plotting a phylogenetic heatmap. It also allows the customization of font size, color schemes, and standardization of the data.
#' The coordinates are returned as a list for use in further plotting or data analysis.
#'
#' @param tree A phylogenetic tree of class `phylo` from the `ape` package.
#' @param X A matrix of data (variables x samples) to be visualized in the heatmap. 
#' @param fsize A numeric vector of length 3 specifying the font size for different components of the plot.
#' @param colors A character vector of colors to be used in the heatmap. Defaults to `NULL`, in which case the `heat.colors` palette is used.
#' @param standardize A logical indicating whether to standardize the data (default is `FALSE`).
#' @param legend A logical indicating whether to display a legend (default is `TRUE`).
#' @param labels A logical indicating whether to display labels (default is `TRUE`).
#' @param split A numeric vector of length 2 indicating how to split the plot. Defaults to `c(0.5, 0.5)`.
#' @param xlim A numeric vector of length 2 specifying the x-axis limits (default is `NULL`).
#' @param ylim A numeric vector of length 2 specifying the y-axis limits (default is `NULL`).
#' @param mar A numeric vector of length 4 specifying the plot margins (default is `rep(1.1, 4)`).
#' @param grid A logical indicating whether to add a grid to the plot (default is `FALSE`).
#' @param ... Additional parameters passed to the `phylogram` function or other plotting functions.
#'
#' @return A list containing two elements:
#' \item{xx}{A numeric vector of x coordinates for the heatmap.}
#' \item{yy}{A numeric vector of y coordinates for the heatmap.}
#'
#' @examples
#' # Example usage:
#' # Load a phylogenetic tree (example from the ape package)
#' library(ape)
#' data(bird.orders)
#' 
#' # Create some random data for the heatmap
#' set.seed(123)
#' X <- matrix(rnorm(100), ncol = 10)
#' 
#' # Get the coordinates for the heatmap
#' coords <- phylo.heatmap.coords(tree = bird.orders, X = X, fsize = 1.5)
#' 
#' # Plot the heatmap (not shown in this example)
#' plot(coords$xx, coords$yy, type = "n")  # Placeholder plot
#' 
#' @export

# Función para obtener las coordenadas x e y del mapa de calor filogenético
phylo.heatmap.coords <- function (tree, X, fsize = 1, colors = NULL, standardize = FALSE, ...) 
{
  # Verifica que el tamaño de la fuente tenga 3 valores
  if (length(fsize) != 3) 
    fsize <- rep(fsize, 3)
  # Si se pasa el parámetro "legend", se extrae
  if (hasArg(legend)) 
    legend <- list(...)$legend
  else legend <- TRUE # La leyenda por defecto es TRUE
  # Si se pasa el parámetro "labels", se extrae
  if (hasArg(labels)) 
    labels <- list(...)$labels
  else labels <- TRUE
  # Si se pasa el parámetro "split", se extrae
  if (hasArg(split)) 
    split <- list(...)$split 
  else split <- c(0.5, 0.5) # Por defecto, el gráfico se divide 50-50
  split <- split/sum(split) # Normaliza el valor de la división
  # Si la matriz X no tiene nombres de columnas, asigna valores predeterminados
  if (is.null(colnames(X))) 
    colnames(X) <- paste("var", 1:ncol(X), sep = "")
  # Si la normalización es requerida, estandariza los datos
  if (standardize) {
    sd <- apply(X, 2, function(x) sqrt(var(x, na.rm = TRUE))) # Calcula la desviación estándar
    X <- (X - matrix(rep(1, Ntip(tree)), Ntip(tree), 1) %*% 
            colMeans(X, na.rm = TRUE))/(matrix(rep(1, Ntip(tree)), 
                                               Ntip(tree), 1) %*% sd)
  }
  # Si se pasa un límite para el eje X, se utiliza; de lo contrario, se define uno predeterminado
  if (hasArg(xlim)) 
    xlim <- list(...)$xlim
  else xlim <- c(-0.5, (2 - 0.5) * split[2]/split[1] + 0.5)
  # Si se pasa un límite para el eje Y, se utiliza; de lo contrario, se define uno predeterminado
  if (hasArg(ylim)) 
    ylim <- list(...)$ylim
  else ylim <- if (legend) 
    c(if (standardize) -0.15 else -0.1, if (labels) 1.1 else 1)
  else c(0, if (labels) 1.1 else 1)
  # Si se pasa un margen para el gráfico, se utiliza; de lo contrario, se define uno predeterminado
  if (hasArg(mar)) 
    mar <- list(...)$mar
  else mar <- rep(1.1, 4)
  # Si no se define una paleta de colores, se utiliza la paleta predeterminada
  if (is.null(colors)) 
    colors <- heat.colors(n = 20)[20:1]
  # Si se pasa la opción para añadir una cuadrícula, se extrae
  if (hasArg(grid)) 
    add.grid <- list(...)$grid
  else add.grid <- FALSE
  
  # Desenreda el árbol filogenético
  cw <- untangle(tree, "read.tree")
  
  # Inicializa el gráfico
  plot.new()
  par(mar = mar)
  plot.window(xlim = xlim, ylim = ylim)
  
  # Dibuja el árbol filogenético
  h <- phylogram(cw, fsize = fsize[1], ...)
  
  # Calcula las coordenadas de inicio y fin del mapa de calor
  START <- h + 1/2 * ((2 - 0.5) * split[2]/split[1] + 0.5 - 
                        h)/(ncol(X) - 1) + 0.5 * strwidth("W") * fsize[1]
  END <- (2 - 0.5) * split[2]/split[1] + 0.5 - 1/2 * ((2 - 0.5) * split[2]/split[1] + 0.5 - START)/(ncol(X) - 1)
  
  # Calcula las coordenadas y para las especies/clases
  nTips <- length(tree$tip.label)
  y <- c(-1/(2 * (nTips - 1)), seq(0, 1, length = nTips) + 2/(2 * (nTips - 1)))
  # Calcula las coordenadas x para las variables
  x <- seq(START, END, by = (END - START)/(ncol(X) - 1))
  # Devuelve las coordenadas
  list(xx = x, yy = y)
}
# Establece el entorno de la función
environment(phylo.heatmap.coords) <- environment(phylo.heatmap)


# ---------- phylo.heatmap.legendmod() --------------------
#' Modified Phylogenetic Heatmap with Custom Legend
#'
#' This function generates a phylogenetic heatmap with an optional legend and additional customization options. It allows for standardization of data, color selection, and the addition of a grid to the plot. The function also includes the ability to customize the font size and labels.
#'
#' @param tree A phylogenetic tree of class `phylo` from the `ape` package.
#' @param X A matrix of data (variables x samples) to be visualized in the heatmap. 
#' @param fsize A numeric vector of length 3 specifying the font size for different components of the plot.
#' @param colors A character vector of colors to be used in the heatmap. Defaults to `NULL`, in which case the `heat.colors` palette is used.
#' @param standardize A logical indicating whether to standardize the data (default is `FALSE`).
#' @param legend A logical indicating whether to display a legend (default is `TRUE`).
#' @param labels A logical indicating whether to display labels (default is `TRUE`).
#' @param split A numeric vector of length 2 indicating how to split the plot. Defaults to `c(0.5, 0.5)`.
#' @param xlim A numeric vector of length 2 specifying the x-axis limits (default is `NULL`).
#' @param ylim A numeric vector of length 2 specifying the y-axis limits (default is `NULL`).
#' @param mar A numeric vector of length 4 specifying the plot margins (default is `rep(1.1, 4)`).
#' @param grid A logical indicating whether to add a grid to the plot (default is `FALSE`).
#' @param ... Additional parameters passed to the `phylogram` function or other plotting functions.
#'
#' @return A plot with a phylogenetic heatmap, a modified legend, and optional customizations.
#'
#' @examples
#' # Example usage:
#' # Load a phylogenetic tree (example from the ape package)
#' library(ape)
#' data(bird.orders)
#' 
#' # Create some random data for the heatmap
#' set.seed(123)
#' X <- matrix(rnorm(100), ncol = 10)
#' 
#' # Plot the modified phylogenetic heatmap with custom legend
#' phylo.heatmap.legendmod(tree = bird.orders, X = X, fsize = 1.5)
#'
#' @export

# Función para crear un mapa de calor filogenético con leyenda interactiva
phylo.heatmap.legendmod <- function (tree, X, fsize = 1, colors = NULL, standardize = FALSE, ...) 
{
  # Verifica que el tamaño de la fuente tenga 3 valores
  if (length(fsize) != 3) 
    fsize <- rep(fsize, 3)
  # Si se pasa el parámetro "legend", se extrae
  if (hasArg(legend)) 
    legend <- list(...)$legend
  else legend <- TRUE
  # Si se pasa el parámetro "labels", se extrae
  if (hasArg(labels)) 
    labels <- list(...)$labels
  else labels <- TRUE
  # Si se pasa el parámetro "split", se extrae
  if (hasArg(split)) 
    split <- list(...)$split
  else split <- c(0.5, 0.5) # Por defecto, el gráfico se divide 50-50
  split <- split/sum(split)# Normaliza el valor de la división
  # Si la matriz X no tiene nombres de columnas, asigna valores predeterminados
  if (is.null(colnames(X))) 
    colnames(X) <- paste("var", 1:ncol(X), sep = "")
  # Si la normalización es requerida, estandariza los datos
  if (standardize) {
    sd <- apply(X, 2, function(x) sqrt(var(x, na.rm = TRUE))) # Calcula la desviación estándar
    X <- (X - matrix(rep(1, Ntip(tree)), Ntip(tree), 1) %*% 
            colMeans(X, na.rm = TRUE))/(matrix(rep(1, Ntip(tree)), Ntip(tree), 1) %*% sd)
  }
  # Si se pasa un límite para el eje X, se utiliza; de lo contrario, se define uno predeterminado
  if (hasArg(xlim)) 
    xlim <- list(...)$xlim
  else xlim <- c(-0.5, (2 - 0.5) * split[2]/split[1] + 0.5)
  # Si se pasa un límite para el eje Y, se utiliza; de lo contrario, se define uno predeterminado
  if (hasArg(ylim)) 
    ylim <- list(...)$ylim
  else ylim <- if (legend) 
    c(if (standardize) -0.15 else -0.1, if (labels) 1.1 else 1)
  else c(0, if (labels) 1.1 else 1)
  # Si se pasa un margen para el gráfico, se utiliza; de lo contrario, se define uno predeterminado
  if (hasArg(mar)) 
    mar <- list(...)$mar
  else mar <- rep(1.1, 4)
  # Si no se define una paleta de colores, se utiliza la paleta predeterminada
  if (is.null(colors)) 
    colors <- heat.colors(n = 20)[20:1]
  # Si se pasa la opción para añadir una cuadrícula, se extrae
  if (hasArg(grid)) 
    add.grid <- list(...)$grid
  else add.grid <- FALSE
  
  # Desenreda el árbol filogenético
  cw <- untangle(tree, "read.tree")
  
  # Inicializa el gráfico
  plot.new()
  par(mar = mar)
  plot.window(xlim = xlim, ylim = ylim)
  
  # Dibuja el árbol filogenético
  h <- phylogram(cw, fsize = fsize[1], ...)
  
  # Calcula las coordenadas de inicio y fin del mapa de calor
  START <- h + 1/2 * ((2 - 0.5) * split[2]/split[1] + 0.5 - h)/(ncol(X) - 1) + 0.5 * strwidth("W") * fsize[1]
  END <- (2 - 0.5) * split[2]/split[1] + 0.5 - 1/2 * ((2 - 0.5) * split[2]/split[1] + 0.5 - START)/(ncol(X) - 1)
  
  # Ajusta los valores de X a las etiquetas del árbol
  X <- X[cw$tip.label, ]
  
  # Dibuja el mapa de calor
  image(x = seq(START, END, by = (END - START)/(ncol(X) - 1)), 
        z = t(X[cw$tip.label, ]), add = TRUE, col = colors, right = T, ...)
  # Si se desea, añade una cuadrícula
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
    add.color.bar(leg = END - START, cols = colors, lims = range(X, na.rm = TRUE), title = if (standardize) "standardized value" else "value", subtitle = if (standardize) "SD units" else "", prompt = F, digits = if (max(abs(X),na.rm = TRUE) < 1) round(log10(1/max(abs(X), na.rm = TRUE))) + 1 else 2, fsize = fsize[3], lwd = 20, outline = T, x = START, y = -1/(2 * (Ntip(cw) - 1)) - 3 * fsize[3] * strheight("W"))
  if (labels) 
    text(x = seq(START, END, by = (END - START)/(ncol(X) - 1)), y = rep(1 + 1/(2 * (Ntip(cw) - 1)) + 0.4 * fsize[2] * strwidth("I"), ncol(X)), colnames(X), srt = 70, adj = c(0,0.5), cex = fsize[2])
  if (any(is.na(X))) {
    ii <- which(is.na(X), arr.ind = TRUE)
    x.na <- seq(START, END, by = (END - START)/(ncol(X) - 1))[ii[, 2]]
    y.na <- seq(0, 1, by = 1/(nrow(X) - 1))[ii[, 1]]
    for (i in 1:length(x.na)) {
      xx <- x.na[i] + c(1/2, -1/2) * (END - START)/(ncol(X) - 1)
      yy <- y.na[i] + c(-1/2, 1/2) * 1/(nrow(X) - 1)
      lines(xx, yy)
    }
  }
}
environment(phylo.heatmap.legendmod) <- environment(phylo.heatmap)


# ------ PASO 5: Analizar datos -----

# Extraer las coordenadas xx y yy de la función personalizada 'phylo.heatmap'
# Esta línea de código llama a la función 'phylo.heatmap.coords' para extraer las coordenadas 'xx' y 'yy',
# necesarias para agregar los puntos de significancia en el mapa de calor filogenético.
# Los parámetros como el tamaño de fuente ('fsize'), los colores ('colors'), la cuadrícula ('grid'), 
# la separación del gráfico ('split'), el grosor de las líneas ('lwd'), los puntos de corte para los colores 
# ('breaks') y los márgenes del gráfico ('mar') son personalizados para adaptarse al estilo deseado.
xx_yy <- phylo.heatmap.coords(
  class_tree_w_g9,                   # Árvore filogenética con grupos de clase
  phylo_heatmap_mat,                  # Matriz de datos para el mapa de calor
  fsize = c(0.8, 0.9, 0.7),          # Tamaño de la fuente para diferentes elementos del gráfico
  colors = plot_colors,              # Colores a usar en el mapa de calor
  grid = T,                          # Agregar cuadrícula al gráfico
  split = c(0.7, 0.3),               # División del gráfico en dos partes
  lwd = 1,                           # Grosor de las líneas
  breaks = heatmap_breaks,            # Valores de corte para los colores del mapa de calor
  mar = c(1.2, 1.2, 1.2, 1.2)        # Márgenes del gráfico
)

# Graficar el mapa de calor filogenético con la leyenda personalizada
# La función 'phylo.heatmap.legendmod' crea el mapa de calor filogenético con las opciones de personalización.
# Los parámetros de la función permiten ajustar el tamaño de fuente, los colores, la cuadrícula, 
# la separación de las columnas y los márgenes.
phylo.heatmap.legendmod(
  class_tree_w_g9,                   # Árvore filogenética con grupos de clase
  phylo_heatmap_mat,                  # Matriz de datos para el mapa de calor
  fsize = c(0.8, 0.9, 0.7),          # Tamaño de la fuente
  colors = plot_colors,              # Colores a usar en el mapa de calor
  grid = T,                          # Agregar cuadrícula al gráfico
  split = c(0.7, 0.3),               # División del gráfico en dos partes
  lwd = 1,                           # Grosor de las líneas
  breaks = heatmap_breaks,            # Puntos de corte para los colores del mapa de calor
  mar = c(1.2, 1.2, 1.2, 1.2)        # Márgenes del gráfico
)

# Bucle para agregar puntos de significancia en el mapa de calor
# Este bucle recorre cada elemento de la matriz de significancia 'sig_mat' para buscar aquellos 
# valores que tienen el símbolo "*" (indicando significancia). Luego, agrega puntos en el gráfico
# en las posiciones correspondientes usando las coordenadas 'xx' y 'yy' extraídas previamente.
# Se asignan colores en función del valor de la matriz de datos del mapa de calor.
for(i in 1:nrow(sig_mat)) {          # Recorrer las filas de la matriz de significancia
  for(j in 1:ncol(sig_mat)) {        # Recorrer las columnas de la matriz de significancia
    if(sig_mat[i,j] == "*") {        # Si el valor en la matriz es "*" (significativo)
      rnm <- rownames(sig_mat)[i]    # Obtener el nombre de la fila (muestra)
      cnm <- colnames(sig_mat)[j]    # Obtener el nombre de la columna (variable)
      
      # Agregar un punto en la posición correspondiente en el gráfico
      points(
        xx_yy$xx[j],                 # Coordenada x extraída previamente
        xx_yy$yy[i],                 # Coordenada y extraída previamente
        cex = 1.5,                   # Tamaño del símbolo
        col = if(phylo_heatmap_mat[rnm, cnm] < heatmap_breaks[3]) "white" else "black", # Color del símbolo dependiendo del valor
        pch = "*"                    # Usar "*" como símbolo para significancia
      )
    }
  }
}


# ------ Figura 3A - barplot -------
# Filtrar clases que tienen más de 9 observaciones de MGE
# Se filtra el dataframe 'all_avg_cnt_tax_class' para mantener solo aquellas clases que aparecen más de 9 veces.
allclass_all <- all_avg_cnt_tax_class %>% 
  filter(mge == all_avg_cnt_tax_class$mge[1]) %>%  # Filtrar por el primer MGE
  group_by(class) %>%                               # Agrupar por clase
  summarise(cntclass = n()) %>%                     # Contar el número de elementos por clase
  filter(., cntclass > 9)                           # Mantener solo clases con más de 9 observaciones

# Filtrar y calcular el promedio de las clases seleccionadas, excluyendo "Cellular" y "Hotspot" de "mge"
# Se filtran las clases para excluir aquellos elementos relacionados con "Cellular" o "Hotspot" en 'mge'.
# Luego, se agrupan por 'mge' y 'class' y se calcula el promedio de 'avg_cnt' para cada clase.
all_avg_cnt_tax_class_sel_all <- all_avg_cnt_tax_class %>% 
  filter(., !grepl("Cellular", mge)) %>%             # Excluir "Cellular" de 'mge'
  filter(., !grepl("Hotspot", mge)) %>%              # Excluir "Hotspot" de 'mge'
  filter(class %in% allclass_all$class) %>%          # Mantener solo clases que aparecen más de 9 veces
  group_by(mge, class) %>%                           # Agrupar por 'mge' y 'class'
  summarise(avg_class = mean(avg_cnt))               # Calcular el promedio de 'avg_cnt' para cada grupo

# Sumar los valores de 'avg_class' por clase para obtener el total de MGE
all_avg_cnt_tax_class_sel_all_mge_all <- all_avg_cnt_tax_class_sel_all %>% 
  group_by(class) %>% 
  summarise(total_mge = sum(avg_class))             # Sumar 'avg_class' por clase para obtener el total de MGE

# Calcular la proporción relativa por clase
# Se combina la información sobre las clases y las proporciones calculadas, dividiendo el valor promedio 
# de 'avg_class' entre el total de MGE para cada clase.
all_avg_cnt_tax_class_sel_4bar_all <- left_join(all_avg_cnt_tax_class_sel_all, all_avg_cnt_tax_class_sel_all_mge_all) %>% 
  group_by(class, mge, total_mge) %>% 
  summarise(frac = avg_class / total_mge)           # Calcular la fracción

# Filtrar los datos para la gráfica de barras, excluyendo las clases con "NA "
barplot_mat <- all_avg_cnt_tax_class_sel_4bar_all %>% 
  filter(., !grepl("NA ", class))                    # Excluir clases con "NA "
barplot_mat$class <- factor(barplot_mat$class, levels = class_tree_w_g9$tip.label)  # Ordenar las clases según el árbol filogenético

# Filtrar los datos para agregar el número total de MGE por clase en la gráfica
barplot_mat_for_n <- all_avg_cnt_tax_class_sel_all_mge_all %>% 
  filter(., !grepl("NA ", class))                    # Excluir clases con "NA "
barplot_mat_for_n$class <- factor(barplot_mat_for_n$class, levels = class_tree_w_g9$tip.label)  # Ordenar las clases

## Graficar el barplot de proporciones relativas por clase ##
# Usamos 'ggplot' para crear una gráfica de barras horizontales ('geom_barh'), mostrando la proporción relativa 
# de cada clase y coloreando según el tipo de MGE.
barplot_specI_class_count_all <- ggplot(barplot_mat, aes(y = class, x = frac, fill = mge)) +  
  geom_barh(stat = "identity", color = "grey60") +  # Graficar barras horizontales con borde gris
  scale_fill_manual("MGE", values = colc)           # Asignar colores personalizados a cada tipo de MGE

# Personalizar el gráfico: definir límites y etiquetas para el eje x, agregar etiquetas con los totales de MGE por clase,
# ajustar la apariencia general, y configurar la leyenda en la parte inferior.
barplot_specI_class_count_all +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1.0), limits = c(0, 1.1)) +  # Definir los valores y límites del eje x
  geom_text(data = barplot_mat_for_n, mapping = aes(y = class, x = 1.05, label = paste("", round(total_mge))), inherit.aes = F, col = "black", size = 6) +  # Etiquetas con los totales de MGE
  theme_cowplot(font_size = 20) +  # Usar tema 'cowplot' con tamaño de fuente 20
  theme(axis.text.x = element_text(angle = 0), legend.position = "bottom") +  # Configurar el texto del eje x y la leyenda en la parte inferior
  theme(axis.text.y = element_text(hjust = 0, vjust = 0.5)) +  # Ajustar la orientación de las etiquetas del eje y
  labs(y = "", x = "relative proportion")   # Etiquetas de los ejes
