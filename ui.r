library(shiny); 

ui <-fluidPage(
  img(src="FishSampling.png"),
  navbarPage("",
             tabPanel(strong("HOME"),                      
                      p(""),
                      strong("This tool informs fish sampling strategies, we present an open web-based application which estimates sampling accuracy related to sea louse abundance, given inputs (such as presumed abundance or clustering effect) provided by the user"),
                      p("Scenario 1-A. With a given sample size, what range of abundance is expected?"),
                      p("Scenario 1-B. How many fish should be sampled to estimate abundance with a given margin of error?"),
                      p("Scenario 2-A. With a given sample size, how probable is a difference between two levels of abundance detected?"),
                      p("Scenario 2-B. How many fish should be sampled to detect a difference between two levels of abundance?"),
                      p("Scenario 3-B. With a given sample size, how probable is a correct decision made to conduct an intervention?"),
                      p("Scenario 3-B. How many fish should be sampled to make a correct decision for an intervention with a given probability?"),
                      p(""),
                      p(""),
                      strong("Table. Summary of Scenarios"),
                      p(""),
                      img(src="Scenario table.png", height = 300, width = 600),
                      p(""),
                      downloadButton(outputId = "downloadData", label = "Code Download"),
                      p(""),
                      p(""),
                      p(""),
                      p(""),
                      img(src="UPEI Logo_no background.png", height = 130, width = 270),
                      img(src="CERC.png", height = 130, width = 270),
                      img(src="ofi.jpg", height = 130, width = 270),
                      p("Please send any comments, questions or suggestions to Jaewoon Jeong (jjeong@upei.ca)")
             ),
             ####################################################################################################################################
             # 1A
             tabPanel(strong("SCENARIO: 1-A"), 
                      h4("With a given sample size, what range of abundance is expected from the sampling?"),
                      fluidRow(
                        column(3,
                               numericInput("total.sample.size.1A","Total Sample Size",30,min=0,max=10000,step=1),
                               numericInput("pen.1A", "Number of pens in a farm:", value=10, min=1, max=100),
                               selectInput("choice.1A", "Numbers of pens to be sampled:", choices =1:100, multiple = TRUE),
                               numericInput("ab.1A", "Abundance:", value = 1.50, min=0, max=100, step=0.01),
                               numericInput("icc.1A", "ICC (Intraclass Clustering Coefficient:", value=0.1, min=0,max=1,step=0.01),
                               numericInput("iteration.1A","Iteration",1000,min=1,max=100000)
                        ),
                        mainPanel(
                          actionButton("goButton.1A", "Go!"),
                          plotOutput("Plot.1A.1"),
                          tableOutput("Table.1A.1")
                        ))),
             ####################################################################################################################################
             # 1B
             tabPanel(strong("SCENARIO: 1-B"), 
                      h4("How many fish should be sampled to estimate abundance with a given margin of error?"),
                      fluidRow(
                        column(3,
                               sliderInput("range.1B", "Range of number of sampled fish per Pen:",min = 1, max = 20,value = c(1,20)),
                               numericInput("pen.1B", "Number of pens in a farm:", value=10, min=1, max=100),
                               selectInput("choice.1B", "Numbers of pens to be sampled:", choices =1:100, multiple = TRUE),
                               numericInput("ab.1B", "Abundance:", value = 1.50, min=0, max=100, step=0.01),
                               numericInput("icc.1B", "ICC (Intraclass Clustering Coefficient:", value=0.1, min=0,max=1,step=0.01),
                               numericInput("u.m","Upper Margin",0.25,min=0,max=10,step=0.01),
                               numericInput("l.m","Lower Margin",0.25,min=0,max=10,step=0.01),
                               sliderInput("cl.1B","Confidence Level:", min=0, max=1,  value = 0.8, step=0.01),
                               numericInput("iteration.1B","Iteration",1000,min=1,max=100000)
                        ),
                        mainPanel(
                          actionButton("goButton.1B", "Go!"),
                          plotOutput("Plot.1B.1"),
                          tableOutput("Table.1B.1")
                        ))),
             ####################################################################################################################################
             # 2A
             tabPanel(strong("SCENARIO: 2-A"),
                      h4("With a given sample size, how probable is a difference between two levels of abundance detected?"), 
                      fluidRow(
                        column(3,
                               numericInput("total.sample.size.2A","Total Sample Size",30,min=0,max=10000,step=1),
                               numericInput("pen.2A", "Number of pens in a farm:", value=10, min=1, max=100),
                               selectInput("choice.2A", "Numbers of pens to be sampled:", choices =1:100, multiple = TRUE),
                               numericInput("ab1.2A", "Lower Abundance:", value = 1.50, min=0, max=10, step=0.01),
                               numericInput("ab2.2A", "Higher Abundance:", value = 2.00, min=0, max=10, step=0.01),
                               selectInput("ab.standard.2A", "Abundance of Interest:", choices = c("Lower Abundance","Higher Abundance")),
                               numericInput("icc.2A", "ICC (Intraclass Clustering Coefficient:", value=0.1, min=0,max=1,step=0.01),
                               numericInput("iteration.2A","Iteration",1000,min=1,max=100000)
                        ),
                        mainPanel(
                          actionButton("goButton.2A", "Go!"),
                          h4(""),
                          plotOutput("Plot.2A.1"),
                          h4("Statistical Power"),
                          tableOutput("Table.2A.1"),
                          h4("Critical Value"),
                          tableOutput("Table.2A.2")
                        ))),
             ####################################################################################################################################
             # 2B
             tabPanel(strong("SCENARIO: 2-B"), 
                      h4("How many fish should be sampled to detect a difference between two levels of abundance?"),
                      fluidRow(
                        column(3,
                               sliderInput("range.2B", "Range of number of sampled fish per Pen:",min = 1, max = 20,value = c(1,20)),
                               numericInput("pen.2B", "Number of pens in a farm:", value=10, min=1, max=100),
                               selectInput("choice.2B", "Numbers of pens to be sampled:", choices =1:100, multiple = TRUE),
                               numericInput("ab1.2B", "Lower Abundance:", value = 1.5, min=0, max=10, step=0.01),
                               numericInput("ab2.2B", "Higher Abundance:", value = 2.0, min=0, max=10, step=0.01),
                               selectInput("ab.standard.2B", "Abundance of Interest:", choices = c("Lower Abundance","Higher Abundance")),
                               numericInput("icc.2B", "ICC (Intraclass Clustering Coefficient:", value=0.1, min=0,max=1,step=0.01),
                               numericInput("iteration.2B","Iteration",1000,min=1,max=100000),
                               sliderInput("cl.2B","Confidence Level:", min=0, max=1,  value = 0.9, step=0.01)
                        ),
                        mainPanel(
                          actionButton("goButton.2B", "Go!"),
                          plotOutput("Plot.2B.1"),
                          tableOutput("Table.2B.1")
                        ))),
             ####################################################################################################################################
             # 3A
             tabPanel(strong("SCENARIO: 3-A"), 
                      h4("With a given sample size, how probable is a correct decision made to conduct an intervention?"),
                      fluidRow(
                        column(3,
                               numericInput("total.sample.size.3A","Total Sample Size",100,min=0,max=1000,step=1),
                               numericInput("pen.3A", "Number of pens in a farm:", value=10, min=1, max=100),
                               selectInput("choice.3A", "Numbers of pens to be sampled:", choices =1:100, multiple = TRUE),
                               numericInput("ab.3A", "Abundance:", value = 0.4, min=0, max=10, step=0.01),
                               numericInput("th.3A","Threshold:",0.5,min=0,max=10,step=0.01),
                               numericInput("icc.3A", "ICC (Intraclass Clustering Coefficient:",value = 0.1, min=0, max=1, step=0.01),
                               numericInput("iteration.3A","Iteration",1000,min=1,max=100000)
                        ),
                        mainPanel(
                          actionButton("goButton.3A", "Go!"),
                          p("Probability of correct decision"),
                          plotOutput("Plot.3A.1"),
                          tableOutput("Table.3A.1")
                        ))),
             ####################################################################################################################################
             # 3B
             tabPanel(strong("SCENARIO: 3-B"), 
                      h4("How many fish should be sampled to make a correct decision for an intervention with a probability?"),
                      fluidRow(
                        column(3,
                               sliderInput("range.3B", "Range of number of sampled fish per Pen:",min = 1, max = 20,value = c(1,20)),
                               numericInput("pen.3B", "Number of pens in a farm:", value=10, min=1, max=100),
                               selectInput("choice.3B", "Numbers of pens to be sampled:", choices =1:100, multiple = TRUE),
                               numericInput("ab.3B", "Abundance:", value = 0.4, min=0, max=10, step=0.01),
                               numericInput("th.3B","Threshold:",0.5,min=0,max=100,step=0.01),
                               numericInput("icc.3B", "ICC (Intraclass Clustering Coefficient:",value = 0.1, min=0, max=1, step=0.01),
                               numericInput("iteration.3B","Iteration",1000,min=1,max=100000),
                               sliderInput("cl.3B","Confidence Level:", min=0, max=1,  value = 0.8, step=0.01)
                        ),
                        mainPanel(
                          actionButton("goButton.3B", "Go!"),
                          p("Probability of correct decision"),
                          plotOutput("Plot.3B.1"),
                          p("Minimum number of sampled fish per pen"),
                          tableOutput("Table.3B.1")
                        ))),
             ###################################################################
             #ICC
             tabPanel(strong("ICC"), 
                      strong("Intraclass correlation coefficient (ICC)"),
                      p("This section provides an simulation of abundance of each pen in a farm to show how abundance of each pen is determined depending on your assumption of ICC value."),
                      p("ICC you input determine the ICC of prevalence of each pen. ICC of abundance of each pen is determined by the stochastic simulation."),
                      fluidRow(
                        column(3,
                               numericInput("fish.0", "Number of sampled fish per pen",value=5, min=2, max=20, step=1),
                               numericInput("pen.0", "Number of pens", value=10, min=2, max=20),
                               numericInput("ab.0", "Abundance:", value = 1.0, min=0, max=100, step=0.01),
                               numericInput("icc.0", "ICC (Intraclass Clustering Coefficient)", value = 0.15, min=0, max=0.5, step=0.01)
                        ),
                        
                        mainPanel(
                          actionButton("goButton.0", "Go!"),
                          p("Abundance of each pen"),
                          plotOutput("Plot0"),
                          tableOutput("Table0")
                        ))),
             ####################################################################################################################################
             # Glossary
             tabPanel(strong("Glossary"), 
                      fluidPage(
                        strong("Abundance:"),
                        p("Number of sea lice per salmon"),
                        strong("Confidence Level: "),
                        p("A percentage that reveals how confident you can be that the population would select an answer within a certain range. For example, a 80% confidence level means that you can be 80% certain the results lie between x and y numbers"),
                        strong("Critical Value: "),
                        p("A level of abundance that is compared to the assumed abundance to determine whether to determine whether there is a difference in two levels of abundance. If your abundance A is greater than the critical value, you can declare that abundance A are significantly higher than abundance B."),
                        strong("Farm: "),
                        p("A farm consists of a group of pens"),
                        strong("ICC (Intraclass correlation coefficient): "), 
                        p("The degree of clustering of sea lice infections in fish within cages. High ICC means that abundance is highly different from pen to pen, while no ICC means that all pens in a farm have exactly same abundance."),
                        strong("Iteration: "),
                        p("How many times the model is simulated to generate the output"),
                        strong("Margin of Error: "),
                        p("The acceptable range of error plus and minus to abundance. The smaller the margin of error, the closer you are to having the exact answer at a given confidence level but you need larger sample size. For information, see below"),
                        strong("Pen:"),
                        p("An enclosure for fish. A farm consists of multiple pens"),
                        strong("Power: "),
                        p("The ability to detect a difference between groups when a difference actually exists. By convention, most literature uses a value of 80~90% power."),
                        strong("Range of Sample Size: "),
                        p("Users can determine the range of number of sampled fish per pen between 1 and 20"),
                        strong("Numbers of pens to be sampled:"),
                        p("Users can determine from how many pens they sample salmon. When salmon are sampled from partial number of pens, ICC may affect the modelling results"),
                        p("")
                      ))
  ) # navbarPage
) # fluidpage

