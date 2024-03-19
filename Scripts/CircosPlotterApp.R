
" This script loads the parsed data file and visualises different types of plots using an R Shiny interface.
Three main functions are used: 
  circos_plotter_all - plots all connections between all individuals or populations, highlights the connections associated with selected individuals or populations
  circos_plotter - takes a main individual or population of interest and plots the connections to chromosomes for a number of selected individuals or populations. 
  circos_table - generates and formats data for output in tabular form
  
The application is developed using R shiny. The function defining the user interface contains several tabs, graphs and tables that are updated in the server. 
The server observes input from the user and updates the seleections and graphs accordingly.


Author: Lovisa Lindquist
Date: 2024-03-19
"

########################################################################################################################3

#install libraries
if (!require("circlize")) install.packages("circlize")
if (!require("shiny")) install.packages("shiny")
if (!require("tibble")) install.packages("tibble")
if (!require("reactable")) install.packages("reactable")
if (!require("shinybusy")) install.packages("shinybusy")
library(circlize)
library(shiny)
library(tibble)
library(reactable)
library(shinybusy)

#read input data
data <- read.csv("Data/Output/Parsed_Ancient_Data.tsv", sep = "\t", header=F)
colnames(data) <- data[1,] #remove header
dat <- data[c(2:length(data[,1])),c(1:10)] #remove additional column
dat$cM <- as.numeric(dat$SegmentLengthM)*100 #create new column for centiMorgan based on column with Morgan
#rename misspelled populations
for (pop in 1:length(dat$ID1_Group)) {
  if (dat$ID1_Group[pop] == "Gernamy") {dat$ID1_Group[pop] <- "Germany"} 
  if (dat$ID2_Group[pop] == "Gernamy") {dat$ID2_Group[pop] <- "Germany"} 
  if (dat$ID1_Group[pop] == "Czechia") {dat$ID1_Group[pop] <- "Czech Republic"} 
  if (dat$ID2_Group[pop] == "Czechia") {dat$ID2_Group[pop] <- "Czech Republic"} 
}
#define a function that prepares data for plotting all connections at once

circos_plotter_all <- function(input_data, selection, level) {
  #define variables to hold information of interest
  ID1 <-c()
  ID2 <-c()
  cM <- c()
  color <- c()
  
  #if at individual level, shorten the data set to reduce computational burden, ideally the entire data set should be used
  if (level == "Individual-Individual") {
    short_data <- input_data[1:100, c(1,4,11)]
    #identify all unique individuals in the data
    unique_id1 <- unique(short_data[,1]) 
    unique_id2 <- unique(short_data[,2])}
  else if (level == "Population-Population") {
    short_data <- input_data[1:1000000, c(2,5,11)]
    unique_id1 <- unique(short_data[,1])
    unique_id2 <- unique(short_data[,2])}
    
    # go through all individuals and generate a temporary dataset for each unique comparison
    for (id1 in unique_id1) {
      subdata <- short_data[short_data[,1] == id1,]
      for (id2 in unique_id2) {
        subsubdata <- subdata[subdata[,2] == id2,]
        if (length(subsubdata$ID1) > 0) {#if there is data, extract the necessary details and assign a colour.
          ID1 <- rbind(ID1, id1)
          ID2 <- rbind(ID2, id2)
          cM <- rbind(cM, sum(subsubdata$cM))
          color <- rbind(color, ifelse(id1 %in% selection, yes = rainbow(5000)[which(selection == id1)*10], no = ifelse(id2 %in% selection, yes = rainbow(5000)[which(selection == id2)*10], no = "grey")))}}}

  # summarise the extracted data in a new data table
  finaldata <- tibble(
    ID1 = ID1,
    ID2 = ID2,
    cM = cM,
    color = color
  )
  return(finaldata)
}

#define a function that prepares data for plotting connections between a subject of interest and subjects of comparison
circos_plotter <- function(input_data, subject, comparison, level) { 
  
  finaldata <- data.frame() #holds data following extraction of relevant rows based on subject, comparison
  score <- c()
  ID2_Group <- c()
  Chr <- c()
  
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
  
  #convert data frame into a matrix and create a colour matrix simultaneously
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
  # return both score matrix and colour matrix
  return(list(mat = mat, col_mat = col_mat))
}

# create a function that prepares data for output in tabular form
circos_table <- function(subject, comparison, level) {
  
  #define variables to hold data
  ID1 = c() #holds name of reference subject
  ID2 = c() #holds name of comparison subject
  total = c() # holds total IBD length (cM)
  avg = c() # holds average IBD length (cM)
  IR = c() # holds inferred relationship based on total IBD length
  
  # The data table depends on the level used for comparison
  if (level == "Individual-Individual") { # at an individual level
    for (i in 1:length(comparison)) {
      tempdata <- dat[c(dat$ID1 == subject & dat$ID2 == comparison[i]),] # save the relevant data
      ID1 <- rbind(ID1, tempdata$ID1[1]) 
      ID2 <- rbind(ID2, tempdata$ID2[1])
      total <- rbind(total, round(sum(as.numeric(tempdata$cM)),3))
      IR <- rbind(IR, c(ifelse(sum(as.numeric(tempdata$cM)) > 2007, yes = "1st degree",no = ifelse(sum(as.numeric(tempdata$cM)) > 900, yes = "2nd degree", no = ifelse(sum(as.numeric(tempdata$cM)) > 400, yes = "3rd degree", no = "4-6th degree")))))
    }
    # append all saved data to summary table
    summary_data <- tibble(
      ID1 = ID1,
      ID2 = ID2,
      color = c(rep("O", length(ID1))),
      Total = total,
      InferredRelation = IR
    )
    
  } else if (level == "Population-Population") { # at a population level
    for (i in 1:length(comparison)) {
      tempdata <- dat[c(dat$ID1_Group == subject & dat$ID2_Group == comparison[i]),] # save the relevant data
      ID1 <- rbind(ID1, tempdata$ID1_Group[1])
      ID2 <- rbind(ID2, tempdata$ID2_Group[1])
      total <- rbind(total, round(sum(as.numeric(tempdata$cM)),3))
      avg <- rbind(avg, c(round(sum(as.numeric(tempdata$cM))/length(tempdata$cM),3)))
    }
    # append all saved data to summary table
    summary_data <- data.frame(
      ID1 = ID1,
      ID2 = ID2,
      color = c(rep("O", length(ID1))),
      Total = total,
      Average = avg
    )

  } else if (level == "Individual-Population") { # at an individual-population level
    for (i in 1:length(comparison)) {
      tempdata <- dat[c(dat$ID1 == subject & dat$ID2_Group == comparison[i]),] # save the relevant data
      ID1 <- rbind(ID1, tempdata$ID1[1])
      ID2 <- rbind(ID2, tempdata$ID2_Group[1])
      total <- rbind(total, round(sum(as.numeric(tempdata$cM)),3))
      avg <- rbind(avg, c(round(sum(as.numeric(tempdata$cM))/length(tempdata$cM),3)))
    }
    # append all saved data to summary table
    summary_data <- tibble(
      ID1 = ID1,
      ID2 = ID2,
      color = c(rep("O", length(ID1))),
      Total = total,
      Average = avg) 

  }
  #return the data table
  return(summary_data)
}


# define function for user interface
ui <- fluidPage(
  titlePanel("IBD Circos Plotter of ancient individuals and populations"), # create a title of application
  
  sidebarLayout(position = "left",
                # create widget panel with three different types of options (plot type, level, reference and comparison used as input in circos_plotter or circos_plotter_all)
                sidebarPanel("",
                             width = 3,
                             
                             selectInput("plot", label = "Plot Type", choices = c("All at once", "One to many"), selected=""),
                             selectInput("level", label = "Level", choices = c("Individual-Individual", "Population-Population", "Individual-Population"), selectize = TRUE),
                             selectizeInput("reference", label = "Selection 1", choices = "", selected = "", options = list(maxOptions = 4104)),
                             checkboxGroupInput("comparison", label = "Selection 2", choices = "", selected = "")),
                
                # create a main panel with a plot and a table
                mainPanel("", width = 8,
                          h4(paste("Welcome to CIRCibd!")),
                          add_loading_state(selector = ".shiny-plot-output", spinner = "circle"),
                          tabsetPanel(id = "tabs",
                                      tabPanel("All at once", plotOutput("circos_1", height = "400px")),
                                      tabPanel("One to many", plotOutput("circos_2", height = "400px"), reactableOutput("summary_data")))
                          
                          
                )
                
  ))


# create function for server function
server <- function(input, output, session) {
  
  # define new options in widget based on previous selections
  observeEvent(input$plot, updateSelectInput(session, "level", "Level", if (input$plot == "All at once") {choices = c("Individual-Individual", "Population-Population")} else if (input$plot == "One to many") {choices = c("Individual-Individual", "Population-Population", "Individual-Population")}, selected = "Individual-Individual"))
  observeEvent(c(input$level, input$plot), updateSelectizeInput(session, "reference", "Selection 1", if (input$plot == "One to many" & input$level == "Individual-Individual") {choices = sort(unique(dat$ID1))} else if (input$plot == "One to many" & input$level == "Population-Population") {choices = sort(unique(dat$ID1_Group))} else if (input$level == "Individual-Population") {choices = sort(unique(dat$ID1))} else {choices = c("NA")}, if (input$plot == "One to many" & input$level == "Population-Population") {selected = sort(unique(dat$ID1_Group))[1]} else if (input$plot == "One to many" & input$level == "Individual-Individual") {selected = sort(unique(dat$ID1))[1]} else if (input$plot == "One to many" & input$level == "Individual-Population") {selected = sort(unique(dat$ID1))[1]} else {selected = "NA"}))
  observeEvent(c(input$plot, input$level, input$reference), updateCheckboxGroupInput(session, "comparison", "Selection 2", if (input$plot == "All at once" & input$level == "Individual-Individual") {choices = sort(unique(dat$ID2))} else if (input$plot == "All at once" & input$level == "Population-Population") {choices = sort(unique(dat$ID2_Group))} else if (input$plot == "One to many" & input$level == "Individual-Individual") {choices = sort(unique(dat$ID2[dat$ID1 == input$reference]))} else if (input$plot == "One to many" & input$level == "Population-Population") {choices = sort(unique(dat$ID2_Group[dat$ID1_Group == input$reference]))} else if (input$plot == "One to many" & input$level == "Individual-Population") {choices = sort(unique(dat$ID2_Group[dat$ID1 == input$reference]))}))
  observeEvent(input$plot, updateTabsetPanel(session, "tabs", if (input$plot == "All at once") {selected = "All at once"} else if (input$plot == "One to many") {selected = "One to many"} else {selected = "All at once"}))

  # identify which type of plot should be used based on plot type and level
  observeEvent(c(input$plot, input$level, input$comparison, input$reference),
               
               if (input$plot == "All at once" & input$level == "Population-Population" & length(input$comparison) > 0) {
               output$circos_1 <- renderPlot({  
               finaldata <- circos_plotter_all(dat, input$comparison, input$level)
               
               circos.clear()
               circos.par(gap.degree = 0.1)
               chordDiagram(finaldata, grid.col="grey", preAllocateTracks = 0.5, symmetric = F, col = finaldata$color, annotationTrack = "grid") 
               circos.track(track.index = 1, panel.fun = function(x, y, ...) {
                 xlim = get.cell.meta.data("xlim")
                 xplot = get.cell.meta.data("xplot")
                 ylim = get.cell.meta.data("ylim")
                 sector.name = get.cell.meta.data("sector.index")
                 if (sector.name %in% finaldata$ID1[finaldata$color != "grey"] | sector.name %in% finaldata$ID2[finaldata$color != "grey"]) {
                   circos.text(mean(xlim), ylim[1], sector.name, facing = "clockwise",
                               niceFacing = TRUE, adj = c(0, 1), cex = 0.8)}}, bg.border = NA)})}
               
               
               else if (input$plot == "All at once" & input$level == "Individual-Individual" & length(input$comparison) > 0) {
                 output$circos_1 <- renderPlot({  
                  finaldata <- circos_plotter_all(dat, input$comparison, input$level) 
                   circos.clear()
                   circos.par(gap.degree = 0.1)
                   chordDiagram(finaldata, grid.col="grey", preAllocateTracks = 0.5, symmetric = F, col = finaldata$color, annotationTrack = "grid") 
                   circos.track(track.index = 1, panel.fun = function(x, y, ...) {
                     xlim = get.cell.meta.data("xlim")
                     xplot = get.cell.meta.data("xplot")
                     ylim = get.cell.meta.data("ylim")
                     sector.name = get.cell.meta.data("sector.index")
                     if (sector.name %in% finaldata$ID1[finaldata$color != "grey"] | sector.name %in% finaldata$ID2[finaldata$color != "grey"]) {
                       circos.text(mean(xlim), ylim[1], sector.name, facing = "clockwise",
                                   niceFacing = TRUE, adj = c(0.5, 1), cex = 0.8)}}, bg.border = NA)})}
               
               
               else if (input$plot == "One to many" & input$level == "Individual-Individual" & length(input$comparison) > 0) {
                output$circos_2 <- renderPlot({
                  validate(
                    need(input$comparison, "Oops! Look like you have not selected an individual of comparison yet!")
                  )
                res <- circos_plotter(dat, subject = input$reference, comparison = input$comparison, level = input$level)
                circos.clear() #make sure old graphs are cleared
                # call chord diagram function, all grids specified as grey, colors based on previously generated matrix
                chordDiagram(res$mat, grid.col="grey", preAllocateTracks = 0.5, col = res$col_mat)  
                title(paste("\nSelected:", input$reference)) # add a title with the reference subject name
                summary_data <- circos_table(input$reference, input$comparison, input$level)
                output$summary_data <- renderReactable({
                  # create drop down menu with detailed data
                  reactable(summary_data, columns = list(color = colDef(name = "", style = function(value, index) {
                    
                    if (summary_data$ID2[index] %in% rownames(as.matrix(res$col_mat))) {
                      color <- as.list(res$col_mat)[index][1]
                    }
                    list(color = color)
                  })), 
                            details = function(index) {
                    data <- dat[c(dat$ID1 == input$reference & dat$ID2 == summary_data$ID2[index]), c(1:2, 10, 3:4, 11, 5:7,9)]
                    data$cM <- round(data$cM, 3)
                    reactable(data, outlined = T)
                  })})})}
               
               else if (input$plot == "One to many" & input$level == "Population-Population" & length(input$comparison) > 0) {
                 output$circos_2 <- renderPlot({
                   validate(
                     need(input$comparison, "Oops! Look like you have not selected a population of comparison yet!")
                   )
                   res <- circos_plotter(dat, subject = input$reference, comparison = input$comparison, level = input$level)
                   circos.clear() #make sure old graphs are cleared
                   # call chord diagram function, all grids specified as grey, colors based on previously generated matrix
                   chordDiagram(res$mat, grid.col="grey", preAllocateTracks = 0.5, col = res$col_mat)  
                   title(paste("\nSelected:", input$reference)) # add a title with the reference subject name
                   summary_data <- circos_table(input$reference, input$comparison, input$level)
                   output$summary_data <- renderReactable({
                     # create drop down menu with detailed data
                     reactable(summary_data, columns = list(color = colDef(name = "", style = function(value, index) {
                       
                       if (summary_data$ID2[index] %in% rownames(as.matrix(res$col_mat))) {
                         color <- as.list(res$col_mat)[index][1]
                       }
                       list(color = color)
                     })), 
                               details = function(index) {
                       data <- dat[c(dat$ID1_Group == input$reference & dat$ID2_Group == summary_data$ID2[index]), c(1:2, 10, 3:4, 11, 5:7,9)] 
                       data$cM <- round(data$cM, 3)
                       reactable(data, outlined = T)
                     })})})}
               
               else if (input$plot == "One to many" & input$level == "Individual-Population" & length(input$comparison) > 0) {
                 output$circos_2 <- renderPlot({
                   validate(
                     need(input$comparison, "Oops! Look like you have not selected a population of comparison yet!")
                   )
                   res <- circos_plotter(dat, subject = input$reference, comparison = input$comparison, level = input$level)
                   
                   circos.clear() #make sure old graphs are cleared
                   # call chord diagram function, all grids specified as grey, colors based on previously generated matrix
                   chordDiagram(res$mat, grid.col="grey", preAllocateTracks = 0.5, col = res$col_mat)  
                   title(paste("\nSelected:", input$reference)) # add a title with the reference subject name
                   summary_data <- circos_table(input$reference, input$comparison, input$level)
                   output$summary_data <- renderReactable({
                     reactable(summary_data, columns = list(color = colDef(name = "", style = function(value, index) {
                       
                       if (summary_data$ID2[index] %in% rownames(as.matrix(res$col_mat))) {
                         color <- as.list(res$col_mat)[index][1]
                       }
                       list(color = color)
                     })), 
                               details = function(index) {
                       data <- dat[c(dat$ID1 == input$reference & dat$ID2_Group == summary_data$ID2[index]), c(1:2, 10, 3:4, 11, 5:7,9)]
                       data$cM <- round(data$cM, 3)
                       reactable(data, outlined = T)})})})}
               
)}


# call application
shinyApp(ui, server)


