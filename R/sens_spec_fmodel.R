# Methods to estimate specificy and sensitivity of an afmodel object
## 20200502 by JJAV
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

#' S3 methods to estimate diagnosis performance of an afmodel
#' 
#' Estimate sensitivity, specificity, positive predicted value and
#' negative predicted value negative predictive value from an afmodel.
#' The estimated "true" negative and "true" positive are estimated using
#' the estimated overall attributable fraction and the predictive positive value
#' associated with each cut-off point as described by
#' Smith, T., Schellenberg, J.A., Hayes, R., 1994. 
#' Attributable fraction estimates and case definitions for malaria 
#' in endemic areas. Stat Med 13, 2345–2358.
#' 
#' @export
#' @param object with the data to calculate the sensitivity and specificity
#' @param ... other parameters for the implementing functions
senspec <- function(object,...) UseMethod("senspec", object)

#' Default implementation method
#' 
#' @rdname senspec
#' @export
senspec.default <- function(object, ...) {
  stop("I don`t know what to do with an object of class", class(object), "\n")
}

#' Estimate sensitivity and specificity of an afmodel
#' 
#' @param cutoff vector of cut-off points to make the estimations
#' @rdname senspec
#' @return a matrix with the columns sensitivity and specificity,
#'   ppv (positive predicted value) and npv (negative predicted value)
#' @importFrom magrittr `%>%`
#' @importFrom dplyr filter
#' @importFrom dplyr mutate
#' @importFrom dplyr summarize
#' @importFrom dplyr arrange
#' @importFrom dplyr select
#' @export
senspec.afmodel <- function(object, cutoff, ...) {
  if (missing(cutoff)) 
    cutoff <- 
      object$data %>%
      dplyr::select(v.density) %>%
      unique() %>% 
      dplyr::filter(v.density > 0) %>%
      dplyr::arrange(v.density) %>% 
      .$v.density %>% 
      as.numeric()
  else 
    cutoff <- 
      sort(cutoff) %>%
      unique()
  
  # Sanity checks
  stopifnot(inherits(object, "afmodel"))
  stopifnot(is.numeric(cutoff))
  stopifnot(length(cutoff) >= 1)
  stopifnot( all(cutoff > 0))
  stopifnot(!any(is.na(cutoff)))
  
  # Order and unique cutoff points
  
  
  
  # Get values from the afmodel object
  l <- object$af
  N <-  
    object$data %>% 
    dplyr::filter(v.fever == 1) %>% 
    nrow()
  
  # Loop to estimate the sensitivities and specificites
  res <- lapply(cutoff, function(x) {
    
    lc <-
      object$data %>%
      dplyr::filter(v.fever == 1) %>%
      dplyr::mutate(tnc =  v.density >= x) %>%
      dplyr::mutate(rrc = ifelse(tnc, rr, NA)) %>%
      dplyr::summarize(lc = mean(rrc, na.rm = T)) %>%
      dplyr::mutate(lc = ifelse(is.nan(lc),0,lc)) %>%
      .$lc %>%
      as.numeric() 
    
    
    nc <- 
      object$data %>% 
      dplyr::filter(v.fever == 1 & v.density >= x) %>% 
      nrow()
    
    #cat("N:", N, " AF:", l, " nc:", nc, " AFc:", lc, "\n")
    
    sens = (nc * lc) / (N * l)
    names(sens) <- NULL
    spec = 1 - (nc * (1 - lc) / (N * (1 - l)))
    names(spec) <- NULL
    #ppv <- (nc*lc) / nc
    ppv <- lc
    npv <- ((N-N*l)-(nc-nc*lc))/(N-nc)
    names(npv) <-NULL
    c(cutoff = x, sensitivity = sens, specificity = spec, ppv = ppv, npv = npv)
  })
  
  res <- do.call(rbind, res)
  res
}

