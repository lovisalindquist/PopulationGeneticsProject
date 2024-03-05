#install libraries
library(circlize)
library(shiny)
library(tibble)
library(DT)
library(reactable)

#read input data
dat <- read.csv("/Users/lovisalindquist/Documents/BINP29/Assignments/Pop_Genetics_Project/Data/IBD/ancient_short_data_transformed.tsv", sep = "\t", header=F)
colnames(dat) <- dat[1,]
dat <- dat[c(2:length(dat[,1])),c(1:6)]
dat$SegmentLengthcM <- as.numeric(dat$SegmentLengthM)*100


# define a function that takes an input individual/group (subject) and the selected comparisons, creates a chord diagram

circos_plotter <- function(input_data, subject, comparison, level) { # requires the input data set and the individual of choice
  
  finaldata <- data.frame()
  if (level == "Individual") {
    subdata <- input_data[input_data$ID1==subject,] # create a new data set with only this individual
    for (i in 1:length(comparison)) {
      finaldata <- rbind(finaldata, subdata[subdata$ID2==comparison[i],])
    }
  } else if (level == "Population") {
    subdata <- input_data[input_data$ID1_Group==subject,] # create a new data set with only this individual
    for (i in 1:length(comparison)) {
      finaldata <- rbind(finaldata, subdata[subdata$ID2_Group==comparison[i],])}}

  # results in a final data set consisting of all rows containing both the selected reference subject and each comparison
    
  # Create a data frame with the relevant information
  if (level == "Individual") {
    df <- data.frame(
      ID = finaldata$ID2,
      chr = paste0("Chr", finaldata$Chr),
      score = finaldata$SegmentLengthcM)
  } else {
    df <- data.frame(
      ID = finaldata$ID2_Group,
      chr = paste0("Chr", finaldata$Chr),
      score = finaldata$SegmentLengthcM)
  }
  # results in a new data frame with only ID of comparison, chromosome and score
  
  #Convert the data frame into a matrix, set row and column names
  mat <- matrix(nrow=length(unique(df$ID)), ncol=length(unique(df$chr)))
  rownames(mat) <- unique(df$ID)
  colnames(mat) <- sort(unique(df$chr))
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
  
  circos.clear()
  chordDiagram(mat, grid.col="grey", preAllocateTracks = 0.5, col = col_mat)
  
  title(paste("\nSelected:", subject))
  
  }


# define function for user interface

ui <- fluidPage(
  titlePanel("IBD Circos Plotter of ancient individuals and populations"),
  
  sidebarLayout(position = "left",
                sidebarPanel("",
                             width = 3,
                             selectInput("level", label = "Level", choices = c("Individual", "Population"), selectize = TRUE),
                             selectInput("reference", label = "Reference", choices = "", selected = ""),
                             checkboxGroupInput("comparison", label = "Comparison", choices = "", selected = "")),
                
                mainPanel("",
                          width = 8,
                          plotOutput("circos"),
                          reactableOutput("summary_data"))
  )
)

# create function for server function

server <- function(input, output, session) {
   
  observeEvent(input$level, updateSelectInput(session, "reference", "Reference", if (input$level == "Individual") {choices = unique(dat$ID1)} else {choices = unique(dat$ID1_Group)}, selected = "I13504"))
  observeEvent(input$reference, updateCheckboxGroupInput(session, "comparison", "Comparison", if (input$level == "Individual") {choices = unique(dat$ID2[dat$ID1 == input$reference])} else {choices = unique(dat$ID2_Group[dat$ID1_Group == input$reference])}, selected = "RISE497_noUDG.SG"))

  output$circos <- renderPlot({
      circos_plotter(dat, input$reference, c(input$comparison), input$level)
  })
  
  output$summary_data <- renderReactable({
    
    
    
    if (input$level == "Individual") {
      ID1 = c()
      ID2 = c()
      total = c()
      for (i in 1:length(input$comparison)) {
        tempdata <- dat[c(dat$ID1 == input$reference & dat$ID2 == input$comparison[i]),]
        ID1 <- rbind(ID1, tempdata$ID1[1])
        ID2 <- rbind(ID2, tempdata$ID2[1])
        total <- rbind(total, sum(as.numeric(tempdata$SegmentLengthcM)))
      }
      summary_data <- tibble(
        ID1 = ID1,
        ID2 = ID2,
        Total = total
      ) 
    } else {
      ID1 = c()
      ID2 = c()
      total = c()
      avg = c()
      for (i in 1:length(input$comparison)) {
        tempdata <- dat[c(dat$ID1_Group == input$reference & dat$ID2_Group == input$comparison[i]),]
        ID1 <- rbind(ID1, tempdata$ID1_Group[1])
        ID2 <- rbind(ID2, tempdata$ID2_Group[1])
        total <- rbind(total, sum(as.numeric(tempdata$SegmentLengthcM)))
        avg <- rbind(avg, c(total/length(tempdata$SegmentLengthcM)))
      }
      print(avg)
      summary_data <- tibble(
        ID1 = ID1,
        ID2 = ID2,
        Total = total,
        Average = avg
      ) 
    }
    

  
    if (input$level == "Individual") {
      reactable(summary_data, details = function(index) {
        data <- dat[c(dat$ID1 == input$reference & dat$ID2 == summary_data$ID2[index]),c(1:5,7)]
        #htmltools::div(style = "padding: 1rem",
        reactable(data, outlined = T)
        
      })
    } else {
      reactable(summary_data, details = function(index) {
        data <- dat[c(dat$ID1_Group == input$reference & dat$ID2_Group == summary_data$ID2[index]), c(1:5,7)]
        reactable(data, outlined = T)
    })
    }
    
  })
}
shinyApp(ui, server)
  
  
  