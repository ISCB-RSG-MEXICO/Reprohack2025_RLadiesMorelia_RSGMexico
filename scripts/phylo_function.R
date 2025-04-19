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
# Establece el entorno de la función
environment(phylo.heatmap.legendmod) <- environment(phylo.heatmap)