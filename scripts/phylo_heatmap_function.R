# ---------------------- duo_maker() -------------------------------
#' Genera combinaciones por pares de elementos separados por un delimitador
#'
#' Esta función toma una cadena de texto con elementos separados por un delimitador
#' (por defecto `":"`) y devuelve todas las combinaciones posibles por pares
#' de esos elementos, concatenadas con el delimitador original. Si la cadena contiene
#' dos elementos o menos, se devuelve sin cambios.
#'
#' @param x Cadena de texto que contiene elementos separados por un delimitador.
#' @param sep Carácter delimitador usado para separar los elementos dentro de la cadena. Por defecto es `":"`.
#'
#' @return Una cadena con las combinaciones por pares separadas por `";"`, o la cadena original si contiene dos elementos o menos.
#'
#' @examples
#' duo_maker("A:B:C")
#' # Devuelve: "A:B;A:C;B:C"
#'
#' duo_maker("A:B")
#' # Devuelve: "A:B"
#'
#' duo_maker("A:B:C:D", sep = ":")
#' # Devuelve: "A:B;A:C;A:D;B:C;B:D;C:D"
#'
#' @export

duo_maker <- function(x, sep = ":") {
  a <- strsplit(x, split = sep)[[1]]
  if(length(a) > 2) {
    aa <- combn(a, 2, paste, collapse = ":")
    paste(aa, collapse = ";")
  } else {
    x
  }
}


# ---------------------- lseq() -------------------------------

#' Genera una secuencia con espaciado logarítmico
#'
#' Esta función genera una secuencia de números con espaciado logarítmico base 10
#' entre los valores `from` y `to`, utilizando `length.out` puntos.
#'
#' @param from Valor inicial de la secuencia (debe ser mayor que 0).
#' @param to Valor final de la secuencia (debe ser mayor que 0).
#' @param length.out Número de valores a generar en la secuencia.
#'
#' @return Un vector numérico con valores espaciados logarítmicamente.
#'
#' @examples
#' lseq(1, 1000, length.out = 4)
#' # Devuelve: 1, 10, 100, 1000
#'
#' lseq(0.1, 10, length.out = 5)
#' # Devuelve: 0.1, 0.316, 1, 3.16, 10
#'
#' @export

lseq <- function(from, to, length.out) {
  # logarithmic spaced sequence
  10^(seq(log10(from), log10(to), length.out = length.out))
}

# ---------------------- phylo.heatmap.orient() -------------------------------
#' Mapa de calor alineado a un árbol filogenético con orientación personalizada
#'
#' Esta función genera un mapa de calor alineado con un árbol filogenético, permitiendo ajustes en la orientación,
#' tamaño de fuente, estandarización de valores y otras opciones de personalización para la visualización.
#'
#' Es una versión modificada de `phylo.heatmap` que permite una orientación alternativa para ajustar la
#' visualización de datos junto a árboles filogenéticos.
#'
#' @param tree Un objeto de clase `phylo` que representa el árbol filogenético.
#' @param X Una matriz de valores numéricos con nombres de fila correspondientes a las etiquetas del árbol.
#' @param fsize Vector numérico de longitud 1 o 3 que indica el tamaño de las fuentes (árbol, etiquetas, leyenda).
#' @param colors Vector de colores para la escala del mapa de calor. Por defecto es `heat.colors(n = 20)[20:1]`.
#' @param standardize Lógico. Si es `TRUE`, los datos se estandarizan por columnas antes de graficar.
#' @param ... Argumentos adicionales para personalizar la visualización (e.g., `legend`, `labels`, `split`, `xlim`, `ylim`, `mar`, `grid`).
#'
#' @return Esta función no retorna un valor, sino que produce una visualización en el dispositivo gráfico.
#'
#' @examples
#' # Suponiendo que tienes un árbol y una matriz de datos con los mismos nombres de especies:
#' library(ape)
#' tree <- rtree(10)
#' X <- matrix(rnorm(100), nrow = 10)
#' rownames(X) <- tree$tip.label
#' phylo.heatmap.orient(tree, X)
#'
#' @seealso [phytools::phylo.heatmap] para la versión original.
#'
#' @importFrom phytools phylogram add.color.bar
#' @importFrom ape untangle
#' @export

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
