#' Cut-off points for densities and fever
#' 
#' Generate the cutoffs at every change of density in the fever, but
#' first category is for density 0, and last category if possible have
#' no subjects with no fever.
#' 
#' @param v.fever numeric vector of 0/1 indicating fever or equivalent
#' @param v.density numeric vector of values >= 0 indicating the density
#' 
#' @importFrom magrittr `%>%`
#' @importFrom dplyr filter
#' @importFrom dplyr mutate
#' @importFrom dplyr arrange
#' @importFrom dplyr select
#' @importFrom dplyr group_by
#' @importFrom dplyr ungroup
#' @importFrom dplyr left_join
#' @importFrom dplyr tally
#' @importFrom tidyr pivot_wider
#' 
#' @export
make_cutoffs <- function(v.fever, v.density, add1 = TRUE) {
  stopifnot(!any(is.na(v.fever)))
  stopifnot(!any(is.na(v.density)))
  stopifnot(length(v.fever) == length(v.density))
  stopifnot(all(v.fever == 1 | v.fever == 0))
  stopifnot(all(v.density >= 0))

  df <- data.frame(
    fever = v.fever, 
    density = v.density) %>% 
    group_by(density, fever) %>%
    tally() %>%
    mutate(category = ifelse(fever ==1,"fever","no_fever")) %>%
    select(-fever) %>%
    pivot_wider(names_from = "category", values_from = "n", values_fill = list(n = 0)) %>%
    ungroup() %>%
    mutate(cumnofev = cumsum(no_fever)) %>%
    mutate(cumnofev_lag = lag(cumnofev, default = 0)) %>%
    mutate(chgnofev = cumnofev-cumnofev_lag) %>%
    mutate(cumfev = cumsum(fever)) %>%
    mutate(cumfev_lag = lag(cumfev, default = 0)) %>%
    mutate(chgfev = cumfev-cumfev_lag) %>%
    mutate(chgnofev_lag = lag(chgnofev, default = 0)) %>%
    mutate(chgfev_lag = lag(chgfev,default = 0)) %>%
    mutate(chg2fev = chgfev - chgfev_lag) %>%
    mutate(chg2nofev = chgnofev - chgnofev_lag) %>%
    filter(chg2fev != 0 & chg2nofev != 0)
    
  if (add1) {
    sort(unique(c(1,df$density)))
  } else
  df$density  
}



#' Make a defined number of categories having similar number of positives in each category
#' 
#' Generate the a define number of categories with similar number of positive
#' in each category
#' 
#' @param v.density numeric vector of values >= 0 indicating the density
#' @param ncat number of category
#' @return a vector with the proposed cutoff points
#' @export
make_n_cutoffs <- function(v.fever, v.density, mintot, add1 = TRUE) {
  
  cutpoints <- make_cutoffs(v.fever,v.density, add1)
  df <- 
      data.frame(
        fever = v.fever, 
        density = v.density) %>%
    mutate(k = cut(density,c(cutpoints,Inf), include.lowest =T, labels = cutpoints)) %>%
    group_by(k,fever) %>%
    tally() %>%
    mutate(category = ifelse(fever ==1,"fever","no_fever")) %>%
    select(-fever) %>%
    pivot_wider(names_from = "category", values_from = "n", values_fill = list(n = 0)) %>%
    mutate(total = fever + no_fever) %>%
    ungroup() %>%
    mutate(cumtotal = cumsum(total)) %>%
    mutate(intpart = floor(cumtotal/mintot)) %>%
    mutate(lagpart = lag(intpart, default = 0)) %>%
    mutate(difpart = intpart - lagpart) %>%
    filter(difpart > 0)
    res <- as.numeric(levels(droplevels(df$k)))
    
  if (add1)
    sort(unique(c(1,res)))
  else
    res
}