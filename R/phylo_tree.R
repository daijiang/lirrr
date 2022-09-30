#' Calculate number of polytomies
#' 
#' Calculate the number of polytomies a phylogeny has.
#' 
#' @param tree A phylogeny with 'phylo' as class.
#' @return A data frame with number of polytomies as columns.
#' @export
#' @examples 
#' library(lirrr) 
#' get_n_polytomy(tree)
#' 
get_n_polytomy = function(tree){
  d = tidytree::as_tibble(tree)
  root_node = which(d$parent == d$node)
  n_terminal_polytomy = d[1:(root_node - 1), ] %>% 
    dplyr::group_by(parent) %>%
    dplyr::count() %>% 
    dplyr::filter(n > 2) %>% 
    nrow()
  n_total_polytomy = d[-root_node, ] %>% 
    dplyr::group_by(parent) %>%
    dplyr::count() %>% 
    dplyr::filter(n > 2) %>% 
    nrow()
  data.frame(n_terminal = n_terminal_polytomy,
             n_internal = n_total_polytomy - n_terminal_polytomy,
             n_total = n_total_polytomy)
}
