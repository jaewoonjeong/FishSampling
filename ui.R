ui <-fluidPage(
  img(src="FishSampling.png", height = 200/2, width = 715/2),
  #img(src="QR code.png", height = 150, width = 150),
  navbarPage("",
             tabPanel(strong("HOME"),                      
                      p(""),
                      strong("This tool informs fish sampling strategies, we present an open web-based application which estimates sampling accuracy related to sea louse abundance, given inputs (such as presumed abundance or clustering effect) provided by the user"),
                      p(""),
                      img(src="salmon farm.png", height = 500/2, width = 715/2),
                      p(""),
                      p("Scenario 1. How probable does the estimated abundance through sampling fall within a certain range?"),
                      p("Scenario 2. How probable is the relative difference between the two abundances correctly determined?"),
                      p("Scenario 3. How probable is abundance estimated to be correctly relative to lice limit?"),
                      p(""),
                      img(src="UPEI Logo_no background.png", height = 130/3, width = 270/3),
                      img(src="CERC.png", height = 100/3, width = 270/3),
                      img(src="ofi.jpg", height = 130/3, width = 270/3),
                      p("Please send any comments, questions or suggestions to Jaewoon Jeong (jjeong@upei.ca)"),
                      p(""),
                      p("R shiny code is available from here, https://github.com/jaewoonjeong/FishSampling")),
             
             tabPanel(strong("SCENARIO: 1"), 
                      h3("Accuracy of estimated abundance"),
                      fluidPage(
                        column(5,radioButtons("clustering1B", "Clustering Effect:",c(No = "no", Yes = "yes"))),
                        column(5,numericInput("iteration.1B","Iteration",value=NULL,min=1,max=100000))
                      ),
                      conditionalPanel(condition = "input.clustering1B == 'no'",div(HTML("<u><b>(1) fish</b></u> are sampled from a farm with <u><b>a true abundance of (2)</b></u>. The abundance estimated through this sampling will fall between <u><b>the limits of (3)</b></u> with <b><u>a probability of (OUTPUT)</b></u>."))),
                      conditionalPanel(condition = "input.clustering1B == 'yes'",div(HTML("<u><b>(1) fish</b></u> are sampled from a farm consisting of <u><b>(4) pens</b></u> with <u><b>true abundances of (2) in each pen</b></u>. The fish are sampled only from <u><b>selected numbers of pens of (5)</b></u>. The abundance estimated through this sampling will fall between <u><b>the limits of (3)</b></u> with <b><u>a probability of (OUTPUT)</b></u>."))),
                      sidebarLayout(
                        sidebarPanel(fluidRow(
                          column(6,
                                 selectInput("TSS.1B", "(1) Total number of sampled fish from a farm:",choices = 1:999,multiple = TRUE),
                                 conditionalPanel(condition = "input.clustering1B == 'no'",numericInput("ab.1B", "(2) Abundance:", value = NULL, min=0, max=100, step=0.01)),
                                 conditionalPanel(condition = "input.clustering1B == 'yes'",uiOutput("ui.1B.1"))),
                          column(6,
                                 numericInput("u.m","(3-1) Higher Limit of Abundance:",NULL,min=0,max=10,step=0.01),
                                 numericInput("l.m","(3-2) Lower Limit of Abundance:",NULL,min=0,max=10,step=0.01),
                                 conditionalPanel(condition = "input.clustering1B == 'yes'",sliderInput("pen.1B", "(4) Total Number of pens in a farm:", value=10, min=2, max=20)),
                                 conditionalPanel(condition = "input.clustering1B == 'yes'",selectInput("choice.1B", "(5) Numbers of selected pens for sampling:", choices =1:20, multiple = TRUE)))
                        ),
                        conditionalPanel(condition = "input.clustering1B == 'yes'",tags$i("The number of abundance must be same to the number of pens in a farm")),
                        conditionalPanel(condition = "input.clustering1B == 'yes'",strong("Average Abundance of Pens:")),
                        conditionalPanel(condition = "input.clustering1B == 'yes'",textOutput('Text.1B.1')),
                        conditionalPanel(condition = "input.clustering1B == 'yes'",strong("ICC:")),
                        conditionalPanel(condition = "input.clustering1B == 'yes'",textOutput("Text.1B.2")),
                        conditionalPanel(condition = "input.clustering1B == 'yes'",plotOutput("Plot.1B.2"))
                        ),
                        mainPanel(
                          actionButton("goButton.1B", "Go!"),
                          h4(""),
                          h4("Table 1.1. Probability of estimating the abundance between two supposed limits"),
                          tableOutput("Table.1B.1"),
                          tags$i("(Probability of corresponding to the region of density plot between the two dashed green vertical lines)"),
                          h4(""),
                          conditionalPanel(condition = "input.clustering1B == 'yes'",h4("Table 1.2. Number of sampled fish per pen")),
                          conditionalPanel(condition = "input.clustering1B == 'yes'",tableOutput("TT.1")),
                          conditionalPanel(condition = "input.clustering1B == 'yes'",
                                           tags$i("(The products of 'Number of sampled pen' and 'Number of sampled fish per pen' may be different to 'Number of total sampled fish' due to a rounding issue. 
                                                  The number of total sampled fish that actually used in simulations was the product, not the number of total fish to sample you supposed.)")),
                          h4(""),
                          plotOutput("Plot.1B.3"),
                          h4("Figure 1.1. Density plot"),
                          tags$i("(The black solid and green dashed lines represent abundance and higher and lower limits, respectively)")))),
             
             tabPanel(strong("SCENARIO: 2"), 
                      h3("Comparison of two abundances"),
                      fluidPage(
                        column(5, radioButtons("clustering2B", "Clustering Effect:",c(No = "no", Yes = "yes"))),
                        column(5, numericInput("iteration.2B","Iteration",value=NULL,min=1,max=100000))
                      ),
                      conditionalPanel(condition = "input.clustering2B == 'no'",div(HTML("<u><b>(1) fish</b></u> are sampled from each of the two farms with <u><b>true abundances of (2)</b></u>. The relative level of the two abundances estimated through this sampling will be correctly determined with <u><b>a probability of (OUTPUT)</b></u>."))),
                      conditionalPanel(condition = "input.clustering2B == 'yes'",div(HTML("<u><b>(1) fish</b></u> are sampled from a farm consisting of <u><b>(3) pens</b></u> with <u><b>true abundances of (2) in each pen</b></u>. The fish are sampled only from selected numbers of <u><b>pens of (4)</b></u>. The relative level of the two abundances estimated through this sampling will be correctly determined with <u><b>a probability of (OUTPUT)</b></u>."))),
                      sidebarLayout(
                        sidebarPanel(fluidRow(
                          column(6,
                                 selectInput("TSS.2B","(1) Total number of sampled fish from a farm:",choices=1:999,multiple=TRUE),
                                 conditionalPanel(condition = "input.clustering2B == 'no'",numericInput("ab1.2B", "(2-1) Abundance1:", value =NULL, min=0, max=10, step=0.01)),
                                 conditionalPanel(condition = "input.clustering2B == 'no'",numericInput("ab2.2B", "(2-2) Abundance2 (Abundance of Interest):", value =NULL, min=0, max=10, step=0.01)),
                                 conditionalPanel(condition = "input.clustering2B == 'yes'",uiOutput("ui.2B.1")),
                                 conditionalPanel(condition = "input.clustering2B == 'yes'",uiOutput("ui.2B.2"))),
                          column(6,
                                 conditionalPanel(condition = "input.clustering2B == 'yes'",sliderInput("pen.2B", "(3) Number of pens in a farm:", value=10, min=2, max=20)),
                                 conditionalPanel(condition = "input.clustering2B == 'yes'",selectInput("choice.2B", "(4) Numbers of selected pens for sampling:", choices =1:20, multiple = TRUE)))
                        ),
                        conditionalPanel(condition = "input.clustering2B == 'yes'",tags$i("The number of abundance must be same to the number of pens in a farm")),
                        conditionalPanel(condition = "input.clustering2B == 'yes'",strong("Average Value of Abundance1 of Pens:")),
                        conditionalPanel(condition = "input.clustering2B == 'yes'",textOutput('Text.2B.1')),
                        conditionalPanel(condition = "input.clustering2B == 'yes'",strong("Average Value of Abundance2 of Pens (Abundance of Interest):")),
                        conditionalPanel(condition = "input.clustering2B == 'yes'",textOutput('Text.2B.2')),
                        conditionalPanel(condition = "input.clustering2B == 'yes'",strong("ICC of Abundance1:")),
                        conditionalPanel(condition = "input.clustering2B == 'yes'",textOutput("Text.2B.3")),
                        conditionalPanel(condition = "input.clustering2B == 'yes'",strong("ICC of Abundance2:")),
                        conditionalPanel(condition = "input.clustering2B == 'yes'",textOutput("Text.2B.4")),
                        conditionalPanel(condition = "input.clustering2B == 'yes'",plotOutput("Plot.2B.1")),
                        conditionalPanel(condition = "input.clustering2B == 'yes'",plotOutput("Plot.2B.2"))
                        ),
                        mainPanel(
                          actionButton("goButton.2B", "Go!"),
                          h4("Table 2.1. Probability of detecting a difference"),
                          tableOutput("Table.2B.1"),
                          #tags$i("(The probability corresponds to the region of density plot to LEFT of the RED dashed vertical line, when Abundance2 is LOWER than Abundance1. The probability corresponds to the region of density plot to RIGHT of the RED dashed vertical line, when Abundance2 is HIGHER than Abundance1.)"),
                          conditionalPanel(condition = "input.clustering2B == 'yes'",h4("Table 2.2. Number of sampled fish per pen")),
                          conditionalPanel(condition = "input.clustering2B == 'yes'",tableOutput("TT.2")),
                          #conditionalPanel(condition = "input.clustering2B == 'yes'",tags$i("(The products of 'Number of sampled pen' and 'Number of sampled fish per pen' may be different to 'Number of total sampled fish' due to a rounding issue. The number of total sampled fish that actually used in simulations was the product, not the number of total fish to sample you supposed.)"))
                          plotOutput("Plot.2B.3"),
                          h4("Figure 2.1. Density plot"),
                          tags$i("(Vertical lines represent the two mean abundance values)")))),
             
             tabPanel(strong("SCENARIO: 3"), 
                      h3("Estimated abundance relative to lice limit"),
                      fluidPage(
                        column(5,radioButtons("clustering3B", "Clustering Effect:",c(No = "no", Yes = "yes"))),
                        column(5,numericInput("iteration.3B","Iteration",value=NULL,min=1,max=100000))
                      ),
                      conditionalPanel(condition = "input.clustering1B == 'no'",div(HTML("<u><b>(1) fish</b></u> are sampled from a farm with <u><b>a true abundance of (2)</b></u>. The relative level of the estimated abundance through this sampling and <u><b>the lice limit of (3)</b></u> will be correctly determined with <u><b>a probability of (OUTPUT)</b></u>."))),
                      conditionalPanel(condition = "input.clustering1B == 'yes'",div(HTML("<u><b>(1) fish</b></u> are sampled from a farm consisting of <u><b>(4) pens</b></u> with <u><b>true abundance of (2) in each pen</b></u>. The fish are sampled only from selected numbers of <u><b>pens of (5)</b></u>. The relative level of the estimated abundance through this sampling and <u><b>the lice limit of (3)</b></u> will be correctly determined with <u><b>a probability of (OUTPUT)</b></u>."))),
                      sidebarLayout(
                        sidebarPanel(fluidRow(
                          column(6,
                                 selectInput("TSS.3B", "(1) Total number of  fish to sample from a farm:",choices = 1:999,multiple = TRUE),
                                 conditionalPanel(condition = "input.clustering3B == 'no'",numericInput("ab.3B", "(2) Abundance:", value =NULL, min=0, max=10, step=0.01)),
                                 conditionalPanel(condition = "input.clustering3B == 'yes'",uiOutput("ui.3B.1"))),
                          column(6,
                                 numericInput("th.3B","(3) Lice Limit:",value=NULL,min=0,max=10,step=0.01),
                                 conditionalPanel(condition = "input.clustering3B == 'yes'",sliderInput("pen.3B", "(4) Number of pens in a farm:", value=10, min=2, max=20)),
                                 conditionalPanel(condition = "input.clustering3B == 'yes'",selectInput("choice.3B", "(5) Numbers of selected pens for sampling:",choices =1:20,multiple = TRUE)))
                        ),
                        conditionalPanel(condition = "input.clustering3B == 'yes'",tags$i("The number of abundance must be same to the number of pens in a farm")),
                        conditionalPanel(condition = "input.clustering3B == 'yes'",strong("Average Abundance of Pens:")),
                        conditionalPanel(condition = "input.clustering3B == 'yes'",textOutput('Text.3B.1')),
                        conditionalPanel(condition = "input.clustering3B == 'yes'",strong("ICC:")),
                        conditionalPanel(condition = "input.clustering3B == 'yes'",textOutput("Text.3B.2")),
                        conditionalPanel(condition = "input.clustering3B == 'yes'",plotOutput("Plot.3B.2"))
                        ),
                        mainPanel(
                          actionButton("goButton.3B", "Go!"),
                          h4("Table 3.1. Probability of estimating abundance correctly relative to threshold"),
                          tableOutput("Table.3B.1"),
                          tags$i("(when Abundance is HIGHER than Threshold, the region of density plot to the RIGHT of the red dashed vertical line. When Abundance is LOWER than Threshold, the region of density plot to the LEFT of the red dashed vertical line)"),
                          p(""),
                          conditionalPanel(condition = "input.clustering3B == 'yes'",h4("Table 3.2. Number of sampled fish per pen")),
                          conditionalPanel(condition = "input.clustering3B == 'yes'",tableOutput("TT.3")),
                          conditionalPanel(condition = "input.clustering3B == 'yes'",tags$i("(The products of 'Number of sampled pen' and 'Number of sampled fish per pen' may be different to 'Number of total sampled fish' due to a rounding issue. The number of total sampled fish that actually used in simulations was the product, not the number of total fish to sample you supposed.)")),
                          plotOutput("Plot.3B.1"),
                          h4("Figure 3.1. Density plot"),
                          tags$i("(The black solid line and red dashed line represent abundance and lice limit, respectively)")))),
             
             tabPanel(strong("Glossary"), 
                      fluidPage(
                        strong("Abundance:"),
                        p("Number of sea lice per salmon"),
                        strong("Clustering"),
                        p("If abundances among pens in a farm are assumed to be different so that you want to include the difference in the outputs, you should choose the clustering effect."),
                        strong("Farm: "),
                        p("A farm consists of a group of pens"),
                        strong("ICC (Intraclass correlation coefficient): "), 
                        p("The degree of clustering of sea lice infections in fish within cages. High ICC means that abundance is highly different from pen to pen, while no ICC means that all pens in a farm have exactly same abundance."),
                        strong("Iteration: "),
                        p("How many times the model is simulated to generate the output"),
                        strong("Pen:"),
                        p("An enclosure for fish. A farm consists of multiple pens"),
                        strong("Numbers of selected pens for sampling:"),
                        p("Users can determine from how many pens they sample salmon. When salmon are sampled from partial number of pens, ICC may affect the modelling results"),
                        strong("Lice Limit"),
                        p("Lice limit is a certain level of abundance, which is not supposed to be exceeded by true abundance. The excess will facilitate an instantaneous treatment of sea lice or harvesting of salmon."),
                        p(""))))
  )