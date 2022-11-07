library(shiny)
library("shinythemes")

shinyUI(fluidPage(theme=shinytheme("darkly"),
                  tabsetPanel(type="tabs",
                              #Output Plots onto actual Tabs.
                              tabPanel("SCREEN.",
                                       fileInput('input$filename$datapath', 'Single Cell Filename'),
                                       actionButton('runPCA', 'PCA.'),
                                       actionButton('normMatrix', 'Normalization.'),
                                       actionButton('createHeat', 'Heatmap.'),
                                       actionButton('createFrequency', 'Frequency.'),
                                       actionButton('createVenn', 'Venn.'),
                                       actionButton('createCorrelation', 'Correlation.')),
                         
                  )
)

)