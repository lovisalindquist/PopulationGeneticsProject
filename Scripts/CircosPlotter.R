#install libraries
library(circlize)
library(shiny)
library(tibble)
library(reactable)

#read input data
dat <- read.csv("Parsed_Ancient_Data.tsv", sep = "\t", header=F)
colnames(dat) <- dat[1,] #remove header
dat <- dat[c(2:length(dat[,1])),c(1:8)] #remove additional column
dat$cM <- as.numeric(dat$SegmentLengthM)*100 #create new column for centiMorgan based on column with Morgan

# define a function that takes an input subject (individual/group) and one or more selected comparisons, creates a chord diagram
circos_plotter <- function(input_data, subject, comparison, level) { 
  
  finaldata <- data.frame() #holds data following extraction of relevant rows based on subject, comparison
  # the data to be extracted differ depending on the requirements of the final output
  if (level == "Individual-Individual") { # for individual level comparisons
    subdata <- input_data[input_data$ID1==subject,] # create a new data set with only the selected reference individual
    for (i in 1:length(comparison)) { 
      finaldata <- rbind(finaldata, subdata[subdata$ID2==comparison[i],])} #only keep the rows for which the reference is compared to individuals specified in comparison
    #create a data frame to be converted into a matrix 
    df <- data.frame(
      ID = finaldata$ID2,
      chr = paste0("Chr", finaldata$Chr),
      score = finaldata$cM)} 
  
  else if (level == "Population-Population") { # for population level comparisons
    subdata <- input_data[input_data$ID1_Group==subject,] # create a new data set with only the reference population
    for (i in 1:length(comparison)) {
      finaldata <- rbind(finaldata, subdata[subdata$ID2_Group==comparison[i],])} #only keep the rows for which the reference is compared to populations specified in comparison
    #create a data frame to be converted into a matrix 
    df <- data.frame(
      ID = finaldata$ID2_Group,
      chr = paste0("Chr", finaldata$Chr),
      score = finaldata$cM)} 
  
  else if (level == "Individual-Population") { # for individual to population level
    subdata <- input_data[input_data$ID1==subject,] # create a new data set with only the reference individual
    for (i in 1:length(comparison)) {
      finaldata <- rbind(finaldata, subdata[subdata$ID2_Group==comparison[i],])}
     #create a data frame to be converted into a matrix 
    score <- c()
    ID2_Group <- c()
    Chr <- c()
    uniqueID2_Group <- unique(finaldata$ID2_Group)
    uniqueChr <- unique(finaldata$Chr)
    for (ID in uniqueID2_Group) {
      for (chrom in uniqueChr) {
        ID2_Group <- rbind(ID2_Group, ID)
        Chr <- rbind(Chr, chrom)
        score <- rbind(score, mean(finaldata$cM[c(finaldata$ID2_Group == ID & finaldata$Chr == chrom)]))
      }}
     df <- data.frame(
      ID = ID2_Group,
      chr = paste0("Chr", Chr),
      score = score)}

  #Convert the data frame into a matrix, set row and column names
  mat <- matrix(nrow=length(unique(df$ID)), ncol=length(unique(df$chr)))
  rownames(mat) <- unique(df$ID)
  colnames(mat) <- sort(unique(df$chr))
  #create a color matrix of the same size, using a color instead of score, used for color of links
  col_mat <- matrix(rainbow(length(rownames(mat))*length(colnames(mat))), nrow=length(rownames(mat)), ncol=length(colnames(mat)))
  rownames(col_mat) <- unique(df$ID)
  colnames(col_mat) <- sort(unique(df$chr))
  
  #For each pairwise comparison, add the score
  for (i in 1:length(unique(df$ID))) { 
    col_mat[i,] <- rainbow(length(rownames(mat)))[i]
    for (j in 1:length(unique(df$chr))) {
      mat[i,j] <- 0
      id = rownames(mat)[i]
      chr = colnames(mat)[j]
      for (row in 1:length(df$ID)) {
        if (df$ID[row] == id & df$chr[row] == chr) {
          mat[i,j] <- as.numeric(df$score[row])
        }
      } 
    }
  }
  
  circos.clear() #make sure old graphs are cleared
  # call chord diagram function, all grids specified as grey, colors based on previosuly generated matrix
  chordDiagram(mat, grid.col="grey", preAllocateTracks = 0.5, col = col_mat)  
  title(paste("\nSelected:", subject)) # add a title with the reference subject name
}

# define function for user interface
ui <- fluidPage(
  titlePanel("IBD Circos Plotter of ancient individuals and populations"), # create a title of application
  
  sidebarLayout(position = "left",
                # create widget panel with three different types of options (level, reference and comparison used as input in circos_plotter)
                sidebarPanel("",
                             width = 3,
                             selectInput("level", label = "Level", choices = c("Individual-Individual", "Population-Population", "Individual-Population"), selectize = TRUE),
                             selectizeInput("reference", label = "Reference", choices = "", selected = "", options = list(maxOptions = 4104)),
                             checkboxGroupInput("comparison", label = "Comparison", choices = "", selected = "")),
                # create a main panel with a plot and a table
                mainPanel("",
                          width = 8,
                          plotOutput("circos", height = "600px"),
                          reactableOutput("summary_data"))
  )
)

# create function for server function
server <- function(input, output, session) {
  # define new options in widget based on previous selections
  observeEvent(input$level, updateSelectizeInput(session, "reference", "Reference", if (input$level == "Individual-Individual") {choices = sort(unique(dat$ID1))} else if (input$level == "Population-Population") {choices = sort(unique(dat$ID1_Group))} else if (input$level == "Individual-Population") {choices = sort(unique(dat$ID1))}, if (input$level == "Population-Population") {selected = sort(unique(dat$ID1_Group))[1]} else {selected = sort(unique(dat$ID1))[1]}))
  observeEvent(input$reference, updateCheckboxGroupInput(session, "comparison", "Comparison", if (input$level == "Individual-Individual") {choices = sort(unique(dat$ID2[dat$ID1 == input$reference]))} else if (input$level == "Population-Population") {choices = sort(unique(dat$ID2_Group[dat$ID1_Group == input$reference]))} else if (input$level == "Individual-Population") {choices = sort(unique(dat$ID2_Group[dat$ID1 == input$reference]))}, if (input$level == "Individual-Individual") {selected = sort(unique(dat$ID2[dat$ID1 == input$reference]))[1]} else if (input$level == "Population-Population") {selected = sort(unique(dat$ID2_Group[dat$ID1_Group == input$reference]))[1]} else if (input$level == "Individual-Population") {sort(unique(dat$ID2_Group[dat$ID1 == input$reference]))[1]}))
  
  # plot chord diagram if all three options have been used
  output$circos <- renderPlot({
    validate(
      need(input$comparison, "Select at least one subject for comparison.")
    )
    circos_plotter(dat, input$reference, c(input$comparison), input$level)
  })
  
  # generate a summary table from data
  output$summary_data <- renderReactable({
    #define variables to store relevant data in before being displayed
    ID1 = c() #holds name of reference subject
    ID2 = c() #holds name of comparison subject
    total = c() # holds total IBD length (cM)
    avg = c() # holds average IBD length (cM)
    IR = c() # holds inferred relationship based on total IBD length
    
    # The data table depends on the level used for comparison
    if (input$level == "Individual-Individual") { # at an individual level
      for (i in 1:length(input$comparison)) {
        tempdata <- dat[c(dat$ID1 == input$reference & dat$ID2 == input$comparison[i]),] # save the relevant data
        ID1 <- rbind(ID1, tempdata$ID1[1]) 
        ID2 <- rbind(ID2, tempdata$ID2[1])
        total <- rbind(total, round(sum(as.numeric(tempdata$cM)),3))
        IR <- rbind(IR, c(ifelse(sum(as.numeric(tempdata$cM)) > 2007, yes = "1st degree",no = ifelse(sum(as.numeric(tempdata$cM)) > 900, yes = "2nd degree", no = ifelse(sum(as.numeric(tempdata$cM)) > 400, yes = "3rd degree", no = "4-6th degree")))))
      }
      # append all saved data to summary table
      summary_data <- tibble(
        ID1 = ID1,
        ID2 = ID2,
        Total = total,
        InferredRelation = IR
      )
      # create drop down menu with detailed data
      reactable(summary_data, details = function(index) {
        data <- dat[c(dat$ID1 == input$reference & dat$ID2 == summary_data$ID2[index]),c(1:7,9)]
        data$cM <- round(data$cM, 3)
        reactable(data, outlined = T)
        })
    } else if (input$level == "Population-Population") { # ata population level
      for (i in 1:length(input$comparison)) {
        tempdata <- dat[c(dat$ID1_Group == input$reference & dat$ID2_Group == input$comparison[i]),] # save the relevant data
        ID1 <- rbind(ID1, tempdata$ID1_Group[1])
        ID2 <- rbind(ID2, tempdata$ID2_Group[1])
        total <- rbind(total, round(sum(as.numeric(tempdata$cM)),3))
        avg <- rbind(avg, c(sum(as.numeric(tempdata$cM))/length(tempdata$cM)))
      }
      # append all saved data to summary table
      summary_data <- tibble(
        ID1 = ID1,
        ID2 = ID2,
        Total = total,
        Average = avg
      ) 
      # create drop down menu with detailed data
      reactable(summary_data, details = function(index) {
        data <- dat[c(dat$ID1_Group == input$reference & dat$ID2_Group == summary_data$ID2[index]), c(1:7,9)] 
        data$cM <- round(data$cM, 3)
        reactable(data, outlined = T)
      })
    } else if (input$level == "Individual-Population") { # at an individual-poulation level
      for (i in 1:length(input$comparison)) {
        tempdata <- dat[c(dat$ID1 == input$reference & dat$ID2_Group == input$comparison[i]),] # save the relevant data
        ID1 <- rbind(ID1, tempdata$ID1[1])
        ID2 <- rbind(ID2, tempdata$ID2_Group[1])
        total <- rbind(total, round(sum(as.numeric(tempdata$cM)),3))
        avg <- rbind(avg, c(sum(as.numeric(tempdata$cM))/length(tempdata$cM)))
      }
      # append all saved data to summary table
      summary_data <- tibble(
        ID1 = ID1,
        ID2 = ID2,
        Total = total,
        Average = avg
      ) 
      # create drop down menu with detailed data
      reactable(summary_data, details = function(index) {
        data <- dat[c(dat$ID1 == input$reference & dat$ID2_Group == summary_data$ID2[index]), c(1:7,9)]
        data$cM <- round(data$cM, 3)
        reactable(data, outlined = T)
      })
    }
  })
}

# call application
shinyApp(ui, server)



