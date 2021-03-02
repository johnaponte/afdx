# Run the shinyApp
# code from https://stackoverflow.com/questions/37830819/developing-shiny-app-as-a-package-and-deploying-it-to-shiny-server
# 20210301 by JJAV
###############################################################################


#' Run the shiny app to make a report from a malaria cross-sectional
#' 
#' @export
#' @importfrom shiny runApp
runShiny <- function() {
  appDir <- system.file("shinyapp", package = "afdx")
  if (appDir == "") {
    stop("Could not find shinyapp. Try re-installing `mypackage`.", call. = FALSE)
  }
  
  shiny::runApp(appDir, display.mode = "normal")
}