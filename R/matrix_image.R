#' Visualize a Matrix as an Image
#'
#' Plots a matrix as a color image with an optional legend, contour lines,
#' and customizable color palette.
#'
#' @param A Matrix. The matrix to visualize.
#' @param x Numeric vector. x-axis values (column positions). Default is
#'   \code{1:ncol(A)}.
#' @param y Numeric vector. y-axis values (row positions). Default is
#'   \code{1:nrow(A)}.
#' @param col Character vector. Colors for the image. Default is a rainbow
#'   palette from blue to red.
#' @param bw Logical. If \code{TRUE}, use a greyscale palette. Default is
#'   \code{FALSE}.
#' @param do.contour Logical. If \code{TRUE}, add contour lines. Default is
#'   \code{FALSE}.
#' @param do.legend Logical. If \code{TRUE}, add a color legend. Default is
#'   \code{TRUE}.
#' @param ... Additional arguments passed to \code{\link[graphics]{image}}.
#'
#' @return Invisibly returns \code{NULL}. Called for its side effect of
#'   producing a plot.
#' @export
matrix.image <- function(A,
                         x = NULL,
                         y = NULL,
                         col = rainbow(100, start = 0.67, end = 0),
                         bw = FALSE,
                         do.contour = FALSE,
                         do.legend = TRUE,
                         ...)
{
 if(do.legend) layout(mat = cbind(matrix(1, 5, 5), rep(2, 5)))
 par(mar = c(6, 5, 3, 2))
 if(is.null(x)) x = 1:ncol(A)
 if(is.null(y)) y = 1:nrow(A)
 nx = length(x)
 ny = length(y)
 x1 = c(1.5 * x[1] - 0.5 * x[2], 1.5 * x[nx] - 0.5 * x[nx - 1])
 y1 = rev(c(1.5 * y[1] - 0.5 * y[2], 1.5 * y[ny] - 0.5 * y[ny - 1]))
 if(bw) col = grey((200:50)/200)
 image(list(x = x, y = y, z = t(A)),
       xlim = x1, ylim = rev(y1), col = col, cex.axis = 1.5, cex.lab = 1.5, bty = "u", ...)
 abline(v = range(x1))
 abline(h = range(y1))
 if(do.contour) contour(x, y, t(A), nlevels = 5, labcex = 1.2, add = TRUE)

 if(do.legend) {
    l.y = seq(min(A), max(A), length.out = 100)
    par(mar = c(6, 2, 3, 1))
    image(list(x = 1:2, y = l.y, z = rbind(l.y, l.y)),
          col = col, bty = "o", xaxt = "n", yaxt = "n")
    axis(side = 2, cex.axis = 1.5, at = pretty(seq(min(A), max(A), length.out = 10)))
 }
}
