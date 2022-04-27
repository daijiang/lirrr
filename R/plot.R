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

#' Highlight all branches linking a subset of species in a phylogeny
#' 
#' @param phy The large phylogeny to plot.
#' @param subset_sp A vector of species names to be highlighted on the plot of phylogeny.
#' @param highlight_color The color to highlight the branches, the default is red.
#' @return A `ggplot2` object.
#' 
highlight_subset_sp_in_phylogeny = function(phy, subset_sp, highlight_color = "red"){
  phy_df = tidytree::as_tibble(phy)
  sp_node = filter(phy_df, label %in% subset_sp)$node
  # all possible combinations
  sp_comb = as.data.frame(t(combn(as.character(sp_node), 2))) %>% 
    set_names(c("s1", "s2")) %>% 
    mutate(s1 = as.integer(s1), s2 = as.integer(s2))
  nodes_to_connect = vector("list", length = nrow(sp_comb))
  for(i in 1:nrow(sp_comb)){
    nodes_to_connect[[i]] = ggtree::get.path(phy, sp_comb[i, 1], sp_comb[i, 2])
  }
  nodes_to_highlight = unique(unlist(nodes_to_connect))
  
  phy_df2 = mutate(phy_df, member = ifelse(node %in% nodes_to_highlight, "y", "n")) %>% 
    tidytree::as.treedata()
  
  p = ggtree(phy_df2, aes(color = member)) +
    scale_color_manual(values = c("n" = "black","y" = highlight_color)) +
    theme(legend.position = "none")
  
  return(p)
}
