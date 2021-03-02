# UI section of the shiny app
# 20200301 by JJAV
################################

library(shiny)

ui <- fluidPage(
  titlePanel("Cross sectional for malaria"),
  sidebarLayout(
    sidebarPanel(
        h3("Parameters"),
        textInput("description","Description",placeholder = "Short description"),
        textInput("country","Country", placeholder = "Country of crossectional"),
        textInput("location","Location", placeholder = "Specific location"),
        textInput("population","Population",placeholder = "Age group of the population included"),
        textInput("timeset","Time", placeholder = "Year(s) when the crossectional was made"),
        h3("Data")

    ),
    mainPanel(
      "Here the report"
    )
  )
)
