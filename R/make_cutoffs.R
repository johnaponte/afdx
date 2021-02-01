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
make_cutoffs <- function(v.fever, v.density) {
  stopifnot(!any(is.na(v.fever)))
  stopifnot(!any(is.na(v.density)))
  stopifnot(length(v.fever) == length(v.density))
  stopifnot(all(v.fever == 1 | v.fever == 0))
  stopifnot(all(v.density >= 0))

  df <- data.frame(fever = v.fever, density = v.density) %>% arrange(density)
 
  cut1 <- 
    df %>%
    filter(fever == 1) %>%
    select(density) %>%
    unique() %>%
    unlist() %>%
    as.numeric() %>%
    c(., 0, max(df$density)) %>%
    unique() %>%
    sort()
  
 
  df2 <-
    df %>%
    mutate(denscat = cut(density,cut1, include.lowest =T)) %>%
    group_by(denscat) %>%
    mutate(mincat = min(density)) %>%
    ungroup()
  
 
  df3 <-
    df2 %>%
    group_by(mincat, fever) %>% 
    tally() %>%
    pivot_wider(names_from = "fever", values_from = "n", values_fill=list(n = 0)) %>%
    filter(`0` == 0) %>%
    mutate(todel = 1) %>%
    select(mincat, todel) %>%
    ungroup()
  
  # Second attempt to have at least one suject in each category
  cutoffs <- 
    df2 %>%
    filter(fever == 1 ) %>%
    left_join(df3, by = "mincat") %>%
    filter(is.na(todel)) %>%
    select(density) %>%
    unlist() %>%
    as.numeric() %>%
    c(., 0,1, max(df$density[df$fever == 0]), max(df$density)) %>%
    unique() %>%
    sort() 
  
  cutoffs <- 
    cutoffs[
      !(cutoffs > max(df$density[df$fever == 0]) & 
          cutoffs < max(df$density))]
  
  cutoffs
}