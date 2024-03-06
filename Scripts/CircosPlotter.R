#install libraries
library(circlize)
library(shiny)
library(tibble)
library(DT)
library(reactable)

#read input data
dat <- read.csv("/Users/lovisalindquist/Documents/BINP29/Assignments/Pop_Genetics_Project/Data/IBD/Ancient_Parsed_Data.tsv", sep = "\t", header=F)
colnames(dat) <- dat[1,]
dat <- dat[c(2:length(dat[,1])),c(1:8)]
dat$cM <- as.numeric(dat$SegmentLengthM)*100

# define a function that takes an input individual/group (subject) and the selected comparisons, creates a chord diagram
subdata <- dat[dat$ID1 == "I18159",]
temp_data <- subdata[subdata$ID2_Group=="Armenia",]
temp_data %>% group_by(ID2_Group, Chr) %>% summarise(cM = sum(cM))
sum(temp_data$cM)


circos_plotter <- function(input_data, subject, comparison, level) { # requires the input data set and the individual of choice
  
  finaldata <- data.frame()
  if (level == "Individual-Individual") {
    subdata <- input_data[input_data$ID1==subject,] # create a new data set with only this individual
    for (i in 1:length(comparison)) {
      finaldata <- rbind(finaldata, subdata[subdata$ID2==comparison[i],])}
    df <- data.frame(
      ID = finaldata$ID2,
      chr = paste0("Chr", finaldata$Chr),
      score = finaldata$cM)} 
  
  else if (level == "Population-Population") {
    subdata <- input_data[input_data$ID1_Group==subject,] # create a new data set with only this individual
    for (i in 1:length(comparison)) {
      finaldata <- rbind(finaldata, subdata[subdata$ID2_Group==comparison[i],])}
    df <- data.frame(
      ID = finaldata$ID2_Group,
      chr = paste0("Chr", finaldata$Chr),
      score = finaldata$cM)} 
  
  else if (level == "Individual-Population") {
    subdata <- input_data[input_data$ID1==subject,] # create a new data set with only this individual
    for (i in 1:length(comparison)) {
      temp_data <- subdata[subdata$ID2_Group==comparison[i],]
      temp_data %>% group_by(ID2_Group, Chr) %>% summarise(cM = sum(cM))
      finaldata <- rbind(finaldata, temp_data)}
    df <- data.frame(
      ID = finaldata$ID2_Group,
      chr = paste0("Chr", finaldata$Chr),
      score = finaldata$cM)}

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
                             selectInput("level", label = "Level", choices = c("Individual-Individual", "Population-Population", "Individual-Population"), selectize = TRUE),
                             selectInput("reference", label = "Reference", choices = "", selected = ""),
                             checkboxGroupInput("comparison", label = "Comparison", choices = "", selected = "")),
                
                mainPanel("",
                          width = 8,
                          plotOutput("circos", height = "600px"),
                          reactableOutput("summary_data"))
  )
)

# create function for server function

server <- function(input, output, session) {
  
  observeEvent(input$level, updateSelectInput(session, "reference", "Reference", if (input$level == "Individual-Individual") {choices = sort(unique(dat$ID1))} else if (input$level == "Population-Population") {choices = sort(unique(dat$ID1_Group))} else if (input$level == "Individual-Population") {choices = sort(unique(dat$ID1))}, if (input$level == "Population-Population") {selected = sort(unique(dat$ID1_Group))[1]} else {selected = sort(unique(dat$ID1))[1]}))
  observeEvent(input$reference, updateCheckboxGroupInput(session, "comparison", "Comparison", if (input$level == "Individual-Individual") {choices = sort(unique(dat$ID2[dat$ID1 == input$reference]))} else if (input$level == "Population-Population") {choices = sort(unique(dat$ID2_Group[dat$ID1_Group == input$reference]))} else if (input$level == "Individual-Comparison") {choices = sort(unique(dat$ID2_Group[dat$ID1 == input$reference]))}, if (input$level == "Individual-Individual") {selected = sort(unique(dat$ID2[dat$ID1 == input$reference]))[1]} else if (input$level == "Population-Population") {selected = sort(unique(dat$ID2_Group[dat$ID1_Group == input$reference]))[1]} else if (input$level == "Individual-Population") {sort(unique(dat$ID2_Group[dat$ID1 == input$reference]))[1]}))
  
  output$circos <- renderPlot({
    validate(
      need(input$comparison, "Select at least one subject for comparison.")
    )
    circos_plotter(dat, input$reference, c(input$comparison), input$level)
  })
  
  output$summary_data <- renderReactable({
    ID1 = c()
    ID2 = c()
    total = c()
    avg = c()
    IR = c()
    if (input$level == "Individual-Individual") {
      
      for (i in 1:length(input$comparison)) {
        tempdata <- dat[c(dat$ID1 == input$reference & dat$ID2 == input$comparison[i]),]
        ID1 <- rbind(ID1, tempdata$ID1[1])
        ID2 <- rbind(ID2, tempdata$ID2[1])
        total <- rbind(total, round(sum(as.numeric(tempdata$cM)),3))
        IR <- rbind(IR, c(ifelse(sum(as.numeric(tempdata$cM)) > 2007, yes = "1st degree",no = ifelse(sum(as.numeric(tempdata$cM)) > 900, yes = "2nd degree", no = ifelse(sum(as.numeric(tempdata$cM)) > 400, yes = "3rd degree", no = "4-6th degree")))))
      }
      summary_data <- tibble(
        ID1 = ID1,
        ID2 = ID2,
        Total = total,
        InferredRelation = IR
      )
      reactable(summary_data, details = function(index) {
        data <- dat[c(dat$ID1 == input$reference & dat$ID2 == summary_data$ID2[index]),c(1:7,9)]
        data$cM <- round(data$cM, 3)
        reactable(data, outlined = T)
        })
    } else if (input$level == "Population-Population") {
      for (i in 1:length(input$comparison)) {
        tempdata <- dat[c(dat$ID1_Group == input$reference & dat$ID2_Group == input$comparison[i]),]
        ID1 <- rbind(ID1, tempdata$ID1_Group[1])
        ID2 <- rbind(ID2, tempdata$ID2_Group[1])
        total <- rbind(total, round(sum(as.numeric(tempdata$cM)),3))
        avg <- rbind(avg, c(sum(as.numeric(tempdata$cM))/length(tempdata$cM)))
      }
      summary_data <- tibble(
        ID1 = ID1,
        ID2 = ID2,
        Total = total,
        Average = avg
      ) 
      reactable(summary_data, details = function(index) {
        data <- dat[c(dat$ID1_Group == input$reference & dat$ID2_Group == summary_data$ID2[index]), c(1:7,9)]
        data$cM <- round(data$cM, 3)
        reactable(data, outlined = T)
      })
    } else if (input$level == "Individual-Population") {
      for (i in 1:length(input$comparison)) {
        tempdata <- dat[c(dat$ID1 == input$reference & dat$ID2_Group == input$comparison[i]),]
        ID1 <- rbind(ID1, tempdata$ID1_Group[1])
        ID2 <- rbind(ID2, tempdata$ID2_Group[1])
        total <- rbind(total, round(sum(as.numeric(tempdata$cM)),3))
        avg <- rbind(avg, c(sum(as.numeric(tempdata$cM))/length(tempdata$cM)))
      }
      summary_data <- tibble(
        ID1 = ID1,
        ID2 = ID2,
        Total = total,
        Average = avg
      ) 
      reactable(summary_data, details = function(index) {
        data <- dat[c(dat$ID1 == input$reference & dat$ID2_Group == summary_data$ID2[index]), c(1:7,9)]
        data$cM <- round(data$cM, 3)
        reactable(data, outlined = T)
      })
    }
  })
}
shinyApp(ui, server)



