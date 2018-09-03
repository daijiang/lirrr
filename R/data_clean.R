#' Randomization tests
#' 
#' Perform randomization test between two vectors.
#' 
#' @param x A numeric vector.
#' @param y A numeric vector.
#' @param n The number of randomization, default is 1000.
#' @export
#' @return  a data frame with mean, rank, p.value, etc.
rand_test = function(x, y, n = 1000) {
    delta.obs = mean(x, na.rm = T) - mean(y, na.rm = T)
    # cat(delta.obs, '\n')
    pool = c(x, y)
    delta.null = plyr::rdply(.n = n, function() {
        x1.index = sample(length(pool), length(x), replace = F)
        x1 = pool[x1.index]
        y1 = pool[-x1.index]
        data.frame(delta.null = mean(x1, na.rm = T) - mean(y1, na.rm = T))
    })$delta.null
    obs.rank = rank(c(delta.obs, delta.null))[1]
    data.frame(mean_x = mean(x, na.rm = T), 
               mean_y = mean(y, na.rm = T), 
               delta = delta.obs, 
               rank = obs.rank, 
               n = n, 
               p.value = ifelse(obs.rank/(n + 1) < 0.5, obs.rank/(n + 1), 1 - obs.rank/(n + 1)))
}

#' Change a distance matrix or a matrix to a data frame
#' 
#' Convert a distance class object or a matrix to a data frame.
#' 
#' @param x An object of class 'dist' or a matrix.
#' @export
#' @return A data frame with three columns: Var1, Var2, distance.
dist_to_df = function(x) {
    mat = as.matrix(x)
    df = reshape2::melt(mat)
    df$Var1 = as.character(df$Var1)
    df$Var2 = as.character(df$Var2)
    df = as.data.frame(t(combn(colnames(mat), 2))) %>% 
      rename(Var1 = V1, Var2 = V2) %>%
      dplyr::left_join(df, by = c("Var1", "Var2"))
    names(df) = c("Var1", "Var2", "distance")
    rownames(df) = 1:nrow(df)
    df
}

#' Function to add a column as row names, and remove it from columns
#' 
#' Convert a column of a data frame to be its row name.
#' 
#' @param df A data frame.
#' @param var The column name (as character) you want to put as row name.
#' @export
#' @return  A data frame.
var_to_rownames = function(df, var = "site") {
    df = as.data.frame(df)
    row.names(df) = df[, var]
    df[, var] = NULL
    df
}

#' Remove sp that not observed at any site (site by sp matrix)
#' 
#' Remove species that not observed in any site.
#' 
#' @param df A data frame in wide form, i.e. site by species data frame, with site names as row name.
#' @export
#' @return  A data frame.
rm_sp_noobs = function(df) {
    if (any(colSums(df) == 0)) {
        df = df[, -which(colSums(df) == 0)]
    }
    df
}

#' Remove site that has no observations (site by sp matrix)
#' 
#' Remove site that has no obsrevations of any species.
#' 
#' @param df A data frame in wide form, i.e. site by species data frame, with site names as row name.
#' @export
#' @return  A data frame.
rm_site_noobs = function(df) {
    if (any(rowSums(df) == 0)) {
        df = df[-which(rowSums(df) == 0), ]
    }
    df
}

#' Logit transformation 
#' 
#' Conduct logit-transformation for proportion data.
#' 
#' @param x A vector of proportions (either in form of 0.20 or 20).
#' @param add_num The number to add for 0s, the dafult is 0.01.
#' @export
#' @return  A vector has been logit transformed.
logit_tran = function(x, add_num = 0.01) {
    if (any(x > 1)) 
        x = x/100  # convert to proportion
    log((x + add_num) / (1 - x + add_num))
}

#' Make the first letter upper case
#' 
#' Make the first letter to be upper case.
#' 
#' @param x A vector of species names.
#' @export
#' @return  A vector.
cap_first_letter <- function(x) {
  sub('^([a-z])', '\\U\\1', x, perl = TRUE)
}

#' Match taxa names with the Open Tree of Life
#' 
#' This function will try to do batch match based on rotl::tnrs_match_names(), 
#'   which only allows <= 250 species each time. 
#' 
#' @param taxa A vector of species names to match.
#' @param n_per_req Number of species to match per requery, must be <= 250.
#' @return A data frame.
#' @export
#' 
tnrs_match_names_2 = function(taxa, n_per_req = 20, ...) {
    n = length(taxa)
    stopifnot(n_per_req <= 250)
    if (n < n_per_req) 
        return(rotl::tnrs_match_names(taxa, ...))
    x = data.frame(sp = taxa, 
                   nitem = c(rep(1:floor(n/n_per_req), each = n_per_req), 
                             rep(ceiling(n/n_per_req), n - (n_per_req * floor(n/n_per_req)))), 
        stringsAsFactors = FALSE)
    out = vector("list", max(x$nitem))
    for (i in 1:(max(x$nitem))) {
        cat(i, " of ", max(x$nitem), "\n")
        out[[i]] = try(rotl::tnrs_match_names(x$sp[x$nitem == i], ...))
    }
    plyr::ldply(out)
}
