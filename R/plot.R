# Panel of pairs() function
#' \code{panel.cor} to plot absolute value of correlations
#' 
#' @export
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, color.threshold = 0.5, ...) {
    usr <- par("usr")
    on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    z = na.omit(data.frame(x = x, y = y))
    r <- cor(z$x, z$y)
    txt <- format(c(r, 0.123456789), digits = digits)[1]
    txt <- paste0(prefix, txt)
    if (missing(cex.cor)) 
        cex.cor <- 0.9/strwidth(txt)
    color.cor = ifelse(abs(r) > color.threshold, "red", "black")
    text(0.5, 0.5, txt, cex = 2, col = color.cor)
    # text(0.5, 0.5, txt, cex = 7*cex.cor * r)
}

# Panel of pairs() function
#' \code{panel.hist} to plot absolute value of correlations
#' 
#' @export
panel.hist <- function(x, ...) {
    usr <- par("usr")
    on.exit(par(usr))
    par(usr = c(usr[1:2], 0, 1.5))
    h <- hist(x, plot = FALSE)
    breaks <- h$breaks
    nB <- length(breaks)
    y <- h$counts
    y <- y/max(y)
    rect(breaks[-nB], 0, breaks[-1], y, col = "lightblue", ...)
}
