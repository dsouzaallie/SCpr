library(shiny)
library("shinythemes")

shinyUI(fluidPage(theme=shinytheme("darkly"),
                  tabsetPanel(type="tabs",
                              #Output Plots onto actual Tabs.
                              tabPanel("SCREEN.",
                                       fileInput('ImputedFile', 'Single Cell Filename'),
                                       actionButton('runSCREEN', 'SCREEN')),
                         
                  )
)

)