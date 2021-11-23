fluidPage(
  img(src="FishSampling.png", height = 200/2, width = 715/2),
  #img(src="QR code.png", height = 150, width = 150),
  navbarPage("",
             tabPanel(strong("HOME"),                      
                      p(""),
                      strong("This tool informs fish sampling strategies, we present an open web-based application which estimates sampling accuracy related to sea louse abundance, given inputs (such as presumed abundance or clustering effect) provided by the user."),
                      p(""),
                      img(src="salmon farm.png", height = 500/2, width = 715/2),
                      p(""),
                      strong("Scenario 1. How probable does the estimated abundance through sampling fall within a certain range?"),
                      tags$i("   Example usage: to figure out the overall accuracy of esatimated abundance"),
                      p(""),
                      strong("Scenario 2. How probable is abundance estimated to be higher than a lice limit?"),
                      tags$i("   Example usage: to determine if abundance is over the lice limit or treatment threshold"),
                      p(""),
                      strong("Scenario 3. How probable is the relative difference between the two abundances correctly determined?"),
                      tags$i("   Example usage: to confirm the difference in abundance between before treatment and after treatment"),
                      p(""),
                      img(src="UPEI Logo_no background.png", height = 130/3, width = 270/3),
                      img(src="CERC.png", height = 100/3, width = 270/3),
                      img(src="ofi.jpg", height = 130/3, width = 270/3),
                      p("Please send any comments, questions or suggestions to Jaewoon Jeong (jjeong@upei.ca)"),
                      p(""),
                      p("R shiny code is available from here, https://github.com/jaewoonjeong/FishSampling")),
             
             tabPanel(strong("SCENARIO 1"), 
                      h3("Accuracy of estimated abundance"),
                      fluidPage(
                        column(2,radioButtons("clustering1", "Clustering Effect:",c(No = "no", Yes = "yes"))),
                        column(2,numericInput("iteration1","Iteration",value=100,min=1,max=100000))
                      ),
                      conditionalPanel(condition = "input.clustering1 == 'no'",div(HTML("<u><b>Total number of fish (A)</b></u> are sampled from a farm with <u><b>a true abundance of (B)</b></u>. 
                                                                                         The abundance estimated through this sampling will fall in <u><b>target limits of (F)</b></u> with <b><u>a probability of (OUTPUT)</b></u>."))),
                      conditionalPanel(condition = "input.clustering1 == 'yes'",div(HTML("<u><b>Total number of fish (A)</b></u> are sampled from a farm consisting of <u><b>(D) pens</b></u> with <u><b>true abundances of (B) in each pen</b></u>. 
                      The fish are sampled only from <u><b>selected numbers of pens of (C)</b></u>. 
                                                                                          The abundance estimated through this sampling will fall in <u><b>target limits of (F)</b></u> with <b><u>a probability of (OUTPUT)</b></u>."))),
                      sidebarLayout(
                        sidebarPanel(fluidRow(
                          column(6,
                                 selectInput("TSS1", "(A) Total number of sampled fish from a farm:",choices = 1:300,selected=c(60),multiple = TRUE),
                                 numericInput("ab1", "(B) Assumed True Abundance:", value = 3.0, min=0, max=100, step=0.01),
                                 numericInput("u.m","(F-1) Higher Target Limit:",3.3,min=0,max=10,step=0.01),
                                 numericInput("l.m","(F-2) Lower Target Limit:",2.7,min=0,max=10,step=0.01)),
                          column(6,
                                 conditionalPanel(condition = "input.clustering1 == 'yes'",selectInput("choice1", "(C) Numbers of selected pens for sampling:", choices =1:20,selected=c(3,5,10), multiple = TRUE)),
                                 conditionalPanel(condition = "input.clustering1 == 'yes'",sliderInput("pen1", "(D) Total Number of pens in a farm:", value=10, min=2, max=20)),
                                 conditionalPanel(condition = "input.clustering1 == 'yes'",sliderInput("icc1", "(E) ICC:", value=0.25, min=0, max=0.35,step = 0.01)))),
                          conditionalPanel(condition = "input.clustering1 == 'yes'",plotOutput("Plot.1B.2")),
                          conditionalPanel(condition = "input.clustering1 == 'yes'",h4(""))),
                        mainPanel(
                          actionButton("goButton1", "Go!", class = "btn-primary"),
                          h4(""),
                          h4("Table 1.1. Probability of estimating the abundance between two target limits"),
                          useShinyalert(),  # Set up shinyalert
                          actionButton("btn1.1", "?"),
                          tableOutput("Table.1.1"),
                          h4(""),
                          conditionalPanel(condition = "input.clustering1 == 'yes'",h4("Table 1.2. Number of sampled fish per pen")),
                          conditionalPanel(condition = "input.clustering1 == 'yes'",
                                           useShinyalert(),  # Set up shinyalert
                                           actionButton("btn1.2", "?")),
                          conditionalPanel(condition = "input.clustering1 == 'yes'",tableOutput("Table.1.2")),
                          h4(""),
                          plotOutput("Plot1"),
                          h4("Figure 1. Distribution of abundance"),
                          useShinyalert(),  # Set up shinyalert
                          actionButton("btn1.3", "?")))),
             
             tabPanel(strong("SCENARIO 2"), 
                      h3("Comparison of an estimated abundance with a certain level of abundance"),
                      fluidPage(
                        column(2,radioButtons("clustering3", "Clustering Effect:",c(No = "no", Yes = "yes"))),
                        column(2,numericInput("iteration3","Iteration",value=100,min=1,max=100000))),
                      conditionalPanel(condition = "input.clustering3 == 'no'",div(HTML("<u><b>Total number of fish (A)</b></u> are sampled from a farm with <u><b>a true abundance of (B)</b></u>. The estimated abundance is higher than <u><b>a certain level of abundance (F)</b></u> with <u><b>a probability of (OUTPUT)</b></u>."))),
                      conditionalPanel(condition = "input.clustering3 == 'yes'",div(HTML("<u><b>Total number of fish (A)</b></u> are sampled from a farm consisting of <u><b>(D) pens</b></u> with <u><b>true abundance of (B) in each pen</b></u>. The fish are sampled only from selected numbers of <u><b>pens of (C)</b></u>. 
                                                                                          The estimated abundance is higher than <u><b>a certain level of abundance (F)</b></u> with <u><b>a probability of (OUTPUT)</b></u>."))),
                      sidebarLayout(
                        sidebarPanel(fluidRow(
                          column(6,
                                 selectInput("TSS3", "(A) Total number of  fish to sample from a farm:",choices = 1:300,selected=c(60),multiple = TRUE),
                                 numericInput("ab3", "(B) Assumed True Abundance:", value =2.6, min=0, max=10, step=0.01),
                                 numericInput("th3","(F) Certain level of abundance:",value=3,min=0,max=10,step=0.01)),
                          column(6,
                                 conditionalPanel(condition = "input.clustering3 == 'yes'",selectInput("choice3", "(C) Numbers of selected pens for sampling:",choices =1:20,selected=c(3,5,10),multiple = TRUE)),
                                 conditionalPanel(condition = "input.clustering3 == 'yes'",sliderInput("pen3", "(D) Number of pens in a farm:", value=10, min=2, max=20)),
                                 conditionalPanel(condition = "input.clustering3 == 'yes'",sliderInput("icc3", "(E) ICC:", value=0.25, min=0, max=0.35,step = 0.01)))),
                          conditionalPanel(condition = "input.clustering3 == 'yes'",plotOutput("Plot.3B.2")),
                          conditionalPanel(condition = "input.clustering3 == 'yes'",h4(""))),
                        mainPanel(
                          actionButton("goButton3", "Go!", class = "btn-primary"),
                            h4("Table 2.1. Probability of estimating abundance higher than a certain level of abundance"),
                            useShinyalert(), 
                            actionButton("btn2.1", "?"),
                            tableOutput("Table.3.1"),
                          conditionalPanel(condition = "input.clustering3 == 'yes'",h4("Table 2.2. Number of sampled fish per pen")),
                          conditionalPanel(condition = "input.clustering3 == 'yes'",
                                           useShinyalert(),  # Set up shinyalert
                                           actionButton("btn2.2", "?")),
                          conditionalPanel(condition = "input.clustering3 == 'yes'",tableOutput("Table.3.2")),
                    plotOutput("Plot3"),
                          h4("Figure 2. Distribution of abundance"),
                          useShinyalert(),  # Set up shinyalert
                          actionButton("btn2.3", "?")))),
             
             tabPanel(strong("SCENARIO 3"), 
                      h3("Comparison of two abundances"),
                      fluidPage(
                        column(2, radioButtons("clustering2", "Clustering Effect:",c(No = "no", Yes = "yes"))),
                        column(2, numericInput("iteration2","Iteration",value=100,min=1,max=100000))
                      ),
                      conditionalPanel(condition = "input.clustering2 == 'no'",div(HTML("<u><b>Total number of fish (A)</b></u> are sampled from each of the two farms with <u><b>true abundances of (B) and (F) </b></u>. 
                                                                                         The relative level of the two abundances estimated through this sampling will be correctly determined with <u><b>a probability of (OUTPUT)</b></u>."))),
                      conditionalPanel(condition = "input.clustering2 == 'yes'",div(HTML("<u><b>Total number of fish (A)</b></u> are sampled from a farm consisting of <u><b>(D) pens</b></u> with <u><b>true abundances of (B) and (F) in each pen</b></u>. 
                                                                                          The fish are sampled only from selected numbers of <u><b>pens of (C)</b></u>. 
                                                                                          The relative level of the two abundances estimated through this sampling will be correctly determined with <u><b>a probability of (OUTPUT)</b></u>."))),
                      sidebarLayout(
                        sidebarPanel(fluidRow(
                          column(6,
                                 selectInput("TSS2","(A) Total number of sampled fish from each farm:",choices=1:300,selected=c(60),multiple=TRUE),
                                 numericInput("ab1.2B", "(B) Assumed True Abundance1:", value =2.6, min=0, max=10, step=0.01),
                                 numericInput("ab2.2B", "(F) Assumed True Abundance2:", value =3.0, min=0, max=10, step=0.01)),
                          column(6,
                                 conditionalPanel(condition = "input.clustering2 == 'yes'",selectInput("choice2", "(C) Numbers of selected pens for sampling:", choices =1:20,selected=c(3,5,10), multiple = TRUE)),
                                 conditionalPanel(condition = "input.clustering2 == 'yes'",sliderInput("pen2", "(D) Total Number of pens in a farm:", value=10, min=2, max=20)),
                                 conditionalPanel(condition = "input.clustering2 == 'yes'",sliderInput("icc1.2B", "(E-1) ICC1:", value=0.25, min=0, max=0.35,step = 0.01)),
                                 conditionalPanel(condition = "input.clustering2 == 'yes'",sliderInput("icc2.2B", "(E-2) ICC2:", value=0.25, min=0, max=0.35,step = 0.01))),
                          column(11,
                                 conditionalPanel(condition = "input.clustering2 == 'yes'",plotOutput("Plot.2B.1")),
                                 conditionalPanel(condition = "input.clustering2 == 'yes'",plotOutput("Plot.2B.2"))))),
                        mainPanel(
                          actionButton("goButton2", "Go!", class = "btn-primary"),
                          h4(""),
                          h4("Table 3.1. Probability of detecting a difference"),
                          useShinyalert(),  # Set up shinyalert
                          actionButton("btn3.1", "?"),
                          tableOutput("Table.2.1"),
                          conditionalPanel(condition = "input.clustering2 == 'yes'",h4("Table 3.2. Number of sampled fish per pen")),
                          conditionalPanel(condition = "input.clustering2 == 'yes'",
                                           useShinyalert(),  # Set up shinyalert
                                           actionButton("btn3.2", "?")),
                          conditionalPanel(condition = "input.clustering2 == 'yes'",tableOutput("Table.2.2")),
                          plotOutput("Plot2"),
                          h4("Figure 3. Distribution of abundance"),
                          useShinyalert(),  # Set up shinyalert
                          actionButton("btn3.2", "?")))),
             
             tabPanel(strong("Glossary"), 
                      fluidPage(
                        strong("Abundance:"),
                        p("Average number of sea lice per salmon. Users should assume a probable abundance to use this application"),
                        strong("Clustering"),
                        p("If abundances among pens in a farm are assumed to be different so that you want to include the difference in the outputs, you should choose the clustering effect."),
                        strong("Farm: "),
                        p("A farm consists of a group of pens"),
                        strong("ICC (Intraclass correlation coefficient): "), 
                        p("The degree of clustering of sea lice infections in fish within cages. High ICC means that abundance is highly different from pen to pen, while no ICC means that all pens in a farm have exactly same abundance."),
                        strong("Iteration: "),
                        p("How many times the model is simulated to generate the output. Higher number of iterations can generate more accurate estimation, but takes longer time of simulation."),
                        strong("Pen:"),
                        p("An enclosure for fish. A farm consists of multiple pens"),
                        strong("Numbers of selected pens for sampling:"),
                        p("Users can determine from how many pens they sample salmon. When salmon are sampled from partial number of pens, ICC may affect the modelling results"),
                        strong("Lice Limit"),
                        p("Lice limit is a certain level of abundance, which is not supposed to be exceeded by true abundance. The excess will facilitate an instantaneous treatment of sea lice or harvesting of salmon."),
                        strong("Treatment Threshold"),
                        p("If abundance exceeds treatment threshold, then treatment should be applied."),
                        p(""))))
)
