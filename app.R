library(shiny)
library(ggplot2)
library(colourpicker) 
library("stringr")
libs <- c("tidyverse", "ggVennDiagram", "BiocManager",
          "DESeq2", "edgeR", "limma","gplots","heatmap3")
library(RColorBrewer)
for (package in libs) {
    suppressPackageStartupMessages(require(package, 
                                           quietly = T, 
                                           character.only = T))
    require(package, character.only = T)
}
options(shiny.maxRequestSize=30*1024^2)
library(DT)
library(shinycssloaders)

# Define UI for application that draws a histogram
ui <- fluidPage(
    titlePanel("BF591-R Final project"),
    br(),
    p("By: Meghana Peshwa"),
    br(),
tabsetPanel(
                tabPanel("Samples",
                     sidebarPanel(
                             fileInput("metadata", "Load metadata", accept=".csv"),
                             actionButton("submit1","Submit", width = "100%")
                         ),
                         mainPanel(
                             tabsetPanel(
                                 tabPanel("Summary",
                                          br(),
                                          dataTableOutput("sample_summary")
                                 ),
                                 tabPanel("Table",
                                          br(),
                                          dataTableOutput("sample_info",width=850)
                                 ),
                                 tabPanel("Plots",
                                          br(),
                                          fluidRow(
                                              column(6,
                                                     plotOutput("plot7")),
                                              column(6,
                                                     plotOutput("plot8")),
                                          ),
                                          fluidRow(
                                              column(6,
                                          plotOutput("plot1")),
                                          column(6,
                                          plotOutput("plot2")),
                                              ),
                                          fluidRow(
                                              column(6,
                                          plotOutput("plot3")),
                                          column(6,
                                          plotOutput("plot4")
                                              )),
                                          
                            )
                        )
                    )
                ),
                tabPanel("Counts",
                     sidebarPanel(
                         fileInput("counts_data", "Load normalized counts matrix", accept=".csv"),
                         sliderInput(inputId = "slider1", min = 0, max = 100,
                        label = "Select the percentile of variance to filter genes:", value = 50, step = 10),
                        sliderInput(inputId = "slider2", min = 0, max = 80,
                         label = "Slider to include genes with at least X samples that are non-zero:", value = 50, step = 1),
                        actionButton("submit2","Plot", icon = icon("chart-line"), width = "100%"),
                     ),
                     mainPanel(
                         tabsetPanel(
                             tabPanel("Table",
                                      br(),
                                      dataTableOutput("counts_summary")
                             ),
                             tabPanel("Scatter plots",
                                      br(),
                                      plotOutput("plot5"),
                                      plotOutput("plot6"),
                             ),
                             tabPanel("Heatmap",
                                      br(),
                                      br(),
                                      plotOutput("heatmap"),
                             ),
                             tabPanel("PCA",
                                      br(),
                                      sidebarPanel(
                                          numericInput(inputId="pc1",label="First principal component",min=1,max=80,step = 1,value=1),
                                          numericInput(inputId="pc2",label="Second principal component",min=1,max=80,step = 1,value=2),
                                      ),
                                      mainPanel(
                                          br(),
                                        plotOutput("pca_plot")
                                      )
                             ),
                        )
                    )
                ),
                tabPanel("DE",
                         sidebarPanel(
                             fileInput("de_results", "Load DE results", accept=".csv"),
                             radioButtons("xvar", "Choose the column for the x-axis",
                                          choices = c("baseMean","log2FoldChange",
                                                      "lfcSE", "stat", "pvalue", "padj"),
                                          selected = "log2FoldChange"),
                             radioButtons("yvar", "Choose the column for the y-axis",
                                          choices = c("baseMean","log2FoldChange",
                                                      "lfcSE", "stat", "pvalue", "padj"),
                                          selected = "padj"),
                             colourInput("col1", "Base point color", "#22577A"),
                             colourInput("col2", "Highlight point color", "#FFCF56"),
                             sliderInput(inputId = "slider3", min = -300, max = 0,
                                         label = "Select the magnitude of the p adjusted coloring:", value = -150, step = -30),
                             actionButton("submit3","Plot", icon = icon("chart-line"), width = "100%")
                         ),
                         mainPanel(
                             tabsetPanel(
                                 tabPanel("Table",
                                          br(),
                                          dataTableOutput("de_table")
                                 ),
                                 tabPanel("Plot",
                                          br(),
                                          br(),
                                          plotOutput("volcano")
                                 )
                             )
                         )
                ),
                tabPanel("Individual gene expression",
                         sidebarPanel(
                             fileInput("counts_data1", "Load normalized counts matrix", accept=".csv"),
                             fileInput("metadata1", "Load metadata", accept=".csv"),
                             actionButton("submit4","Submit", width = "100%"),
                         ),
                         mainPanel(
                             tabsetPanel(
                                 tabPanel("Plot",
                                          br(),
                                          fluidRow(
                                              column(3,
                                          selectizeInput(
                                              'gene', 'Select gene', choices = NULL
                                          )),
                                          column(3,
                                          selectInput("categorical_input", "Select categorical field", choices = c("timepoint", "treatment", "sex","lifestage"))),
                                          column(6,
                                          radioButtons("plot_type", "Choose the type of plot",
                                                       choices = c("Boxplot","Bar plot",
                                                                   "Violin plot", "Beeswarm plot"),
                                                       selected = "Violin plot",inline=T))),
                                          fluidRow(column(12, align="center",
                                          actionButton("submit5","Plot", icon = icon("chart-line"), width = "30%"))),
                                          br(),
                                          plotOutput("de")
                                 ),
                             )
                         )
                )
        
    )
)

# Define server logic required to draw a histogram
server <- function(input, output, session) {
    load_metadata <- reactive({
        if(is.null(input$metadata$datapath)) {
            return(NULL)
        }
        metadata <- read.csv(input$metadata$datapath,header=TRUE)
        metadata_updated <- transform(metadata, AvgSpotLen = as.numeric(AvgSpotLen), 
                                      Bases = as.numeric(Bases),
                                      Bytes = as.numeric(Bytes))
        return(metadata_updated)
    })
    
    load_counts_data <- reactive({
        if(is.null(input$counts_data$datapath)) {
            return(NULL)
        }
        counts_data <- read.csv(input$counts_data$datapath,header=TRUE)
        return(counts_data[,-c(1)])
    })
    load_de_results <- reactive({
        if(is.null(input$de_results$datapath)) {
            return(NULL)
        }
        de_results <- read.csv(input$de_results$datapath,header=TRUE)
        return(de_results[,-c(1)])
    })
    load_metadata1 <- reactive({
        if(is.null(input$metadata1$datapath)) {
            return(NULL)
        }
        metadata <- read.csv(input$metadata1$datapath,header=TRUE)
        metadata_updated <- transform(metadata, AvgSpotLen = as.numeric(AvgSpotLen), 
                                      Bases = as.numeric(Bases),
                                      Bytes = as.numeric(Bytes))
        return(metadata_updated)
    })
    
    load_counts_data1 <- reactive({
        if(is.null(input$counts_data1$datapath)) {
            return(NULL)
        }
        counts_data <- read.csv(input$counts_data1$datapath,header=TRUE)
        return(counts_data[,-c(1)])
    })
    
    sample_summary <-
        function(metadata) {
            if(is.null(metadata))
                return(NULL)
            metadata_updated = metadata[,-1]
            metadata_updated = metadata_updated[,-5]
            metadata_updated = metadata_updated[,-11]
            metadata_updated = metadata_updated[,-11]
            metadata_updated = metadata_updated[,-18]
            metadata_updated = metadata_updated[,-18]
            column_names = names(metadata_updated)
            column_type = vector()
            unique_values = vector()
            metadata_updated <- transform(metadata_updated, AvgSpotLen = as.numeric(AvgSpotLen), 
                                          Bases = as.numeric(Bases),
                                          Bytes = as.numeric(Bytes))
            for(i in names(metadata_updated)){
                column_type = append(column_type,class(metadata_updated[[i]]))
                if(class(metadata_updated[[i]])=="character") {
                    unique_values <- append(unique_values,paste(unique(metadata_updated[[i]]), collapse=";"))
                }
                else{
                    unique_values <- append(unique_values,mean(metadata_updated[[i]]))
                }
            }
            sample_summary <- data.frame(
                "Column names"=column_names,
                "Column type"=column_type,
                "Unique values"=unique_values,
                check.names = FALSE
            )
            return(sample_summary)
        }
    filter_counts <- function(counts_data, zero_filter, variance_filter) {
        if(is.null(counts_data))
            return(NULL)
        counts_data_filtered <- counts_data[rowSums(counts_data[-1]!=0)>=zero_filter,]
        data <- counts_data_filtered[,-c(1)]
        variance <- apply(data, 1, var)
        counts_data_filtered_final <- counts_data_filtered[variance>=quantile(variance,variance_filter/100,na.rm=TRUE),]
        summary <- data.frame(
            "Summary type"=c("No. of samples","Total number of genes","% genes passing filter","% genes not passing filter"),
            "Value"=c(length(counts_data_filtered)-1,nrow(counts_data),round(nrow(counts_data_filtered_final)/nrow(counts_data)*100,2),
                         round((nrow(counts_data)-nrow(counts_data_filtered_final))/nrow(counts_data)*100,2)),
            check.names = FALSE
        )
        return(list(counts_data_filtered_final,summary))
    }
    plot_variance_vs_median <- function(data, scale_y_axis=FALSE, title="") {
        median <- apply(data[, -c(ncol(data))], 1, median)
        variance <- apply(data[, -c(ncol(data))], 1, var)
        plot_data <- tibble(median=median, variance=variance, filtered = data$filtered)
        plot_data$rank <- rank(plot_data$median)
        plot <- ggplot(plot_data, aes(x=rank, y=variance)) +
            geom_point(mapping=aes(x=rank, y=variance, color=filtered)) +
            xlab("Rank(median)") +
            ylab("Variance") +
            ggtitle(title)
        if (scale_y_axis) {
            plot <- plot + scale_y_log10()
        }
        return(plot)
    }
    plot_zeros_vs_median <- function(data, title="") {
        zeros <- rowSums(data[, -c(ncol(data))] == 0)
        median <- apply(data[, -c(ncol(data))], 1, median)
        plot_data <- tibble(zeros=zeros, median=median, filtered = data$filtered)
        plot_data$rank <- rank(plot_data$median)
        plot <- ggplot(plot_data, aes(x=rank, y=zeros)) +
            geom_point(mapping=aes(x=rank, y=zeros, color=filtered)) +
            xlab("Rank(median)") +
            ylab("Number of zeros") + 
            ggtitle(title)
        return(plot)
    }
    timepoint_from_sample <- function(str) {
        return(substr(str,3,4))
    }
    
    sample_treatment <- function(str) {
        if(substr(str,1,1)=='F')
            return("Fluctuating")
        else
            return("Control")
    }
    
    meta_info_from_labels <- function(sample_names) {
        timepoint <- sapply(sample_names, timepoint_from_sample)
        treatment <- sapply(sample_names, sample_treatment)
        lifestage <- sapply(sample_names, lifestage_from_sample)
        result <- tibble(
            sample=sample_names,
            timepoint=timepoint,
            treatment=treatment,
            lifestage=lifestage
        )
        return(result)
    }
    plot_pca <- function(data, meta, pc1, pc2, title="") {
        pca <- prcomp(t(data))
        meta$PC1 <- pca$x[ , pc1]
        meta$PC2 <- pca$x[ , pc2]
        variance <- pca$sdev^2 / sum( pca$sdev^2 )
        pca_plot <- ggplot(meta, aes(x=PC1, y=PC2, col=lifestage)) +
            geom_point() +
            xlab(paste0("PC",toString(pc1),": ",round(variance[pc1] * 100),"% variance")) +
            ylab(paste0("PC",toString(pc2),": ",round(variance[pc2] * 100),"% variance")) +
            ggtitle(title)
        return(pca_plot)
    }
    volcano_plot <-
        function(dataf, x_name, y_name, slider, color1, color2) {
            if(is.null(dataf))
                return(NULL)
            plot <- ggplot(dataf, aes(x = !!sym(x_name),
                                      y = -log10(!!sym(y_name)))) +
                geom_point(aes(color = !!sym(y_name) < 1 * 10 ^ (as.numeric(slider)))) +
                theme_minimal() +
                theme_bw() +
                scale_color_manual(values=c(color1, color2)) 
            plot$labels$colour <- paste0(y_name," < 1 x 10^",toString(slider))
            return(plot)
        }
    draw_table <- function(dataf, slider) {
        if(is.null(dataf)) {
            return(NULL)
        }
        #dataf<-na.omit(dataf)
        dataf <- dataf[which(dataf$padj < 1 * 10 ^ (as.numeric(slider))),]
        #dataf<-na.omit(dataf)
        if(nrow(dataf)==0) {
            return(dataf)
        }
        dataf$pvalue <- formatC(dataf$pvalue,digits = -2)
        dataf$padj <- formatC(dataf$padj,digits = -2)
        return(dataf)
    }
    sex_from_sample <- function(str) {
        if(nchar(str)==4){
            return('Empty (larvae)')
        }
        if(substr(str,nchar(str),nchar(str))=='M')
            return("Male")
        else
            return("Female")
    }
    
    lifestage_from_sample <- function(str) {
        if(substr(str,2,2)=='A')
            return("Adult")
        else
            return("Larval")
    }
    individual_gene_expression <- function(counts_data,sample_data, gene) {
        gene_data <- t(counts_data[counts_data$gene_id == gene,][,-c(1)])
        cols <- colnames(counts_data)[colnames(counts_data) != "gene_id"]
        timepoint <- sapply(cols, timepoint_from_sample)
        treatment <- sapply(cols, sample_treatment)
        sex <- sapply(cols, sex_from_sample)
        lifestage <- sapply(cols, lifestage_from_sample)
        result <- tibble(
            sample=cols,
            counts=gene_data[,1],
            timepoint=c(timepoint),
            treatment=c(treatment),
            sex=c(sex),
            lifestage=c(lifestage)
        )
        return(result)
    }
    observeEvent(input$submit1,{
        dataf<-load_metadata()
        sample_summary <- sample_summary(dataf)
        coul <- brewer.pal(5, "Set2") 
        if(!is.null(dataf)) {
        output$plot1 <- renderPlot(barplot(table(dataf$lifestage),
                                           main="Bar chart for lifestage",
                                           ylab="Count",col=coul,xlab="Lifestage")
        )
        output$plot2 <- renderPlot(barplot(table(dataf$timepoint),
                                           main="Bar chart for timepoint",
                                           ylab="Count",col=coul,xlab="Timepoint")
        )
        output$plot3 <- renderPlot(barplot(table(dataf$treatment),
                                           main="Bar chart for treatment",
                                           ylab="Count",col=coul,xlab="Treatment")
        )
        output$plot4 <- renderPlot(barplot(table(dataf$sex),
                                           main="Bar chart for sex",
                                           ylab="Count",col=coul,xlab="Sex",names.arg=c("Empty","female","male"))
        )
        output$plot7 <- renderPlot(ggplot(dataf)+geom_violin(aes(x=treatment,y=Bases,fill=treatment))+
                                       ggtitle("Violin plot of bases by treatment type"))
        output$plot8 <- renderPlot(ggplot(dataf)+geom_violin(aes(x=treatment,y=Bytes,fill=treatment))+
                                       ggtitle("Violin plot of bytes by treatment type"))
        }
        output$sample_info <- DT::renderDataTable(DT::datatable(dataf, 
                                    extensions = 'Buttons', class = "display",
                                    options = list(paging = TRUE, searching = TRUE,
                                                   fixedColumns = TRUE, autoWidth = TRUE,
                                                   ordering = TRUE, dom = 'Bfrtip',
                                                   buttons = c('copy', 'csv'),scrollX = TRUE)))
        output$sample_summary <-  DT::renderDataTable(DT::datatable(sample_summary, 
                                    extensions = 'Buttons', class = "display",
                                    options = list(paging = TRUE, searching = TRUE,
                                                   fixedColumns = TRUE, autoWidth = TRUE,
                                                   ordering = TRUE, dom = 'Bfrtip',
                                                   buttons = c('copy', 'csv'))))
    })
    observeEvent(input$submit2,{
        counts_data <- load_counts_data()
        results <- filter_counts(counts_data,input$slider1,input$slider2)
        
        output$counts_summary <- DT::renderDataTable(DT::datatable(results[[2]], 
                                   extensions = 'Buttons', class = "display",
                                   options = list(paging = TRUE, searching = TRUE,
                                                  fixedColumns = TRUE, autoWidth = TRUE,
                                                  ordering = TRUE, dom = 'Bfrtip',
                                                  buttons = c('copy', 'csv'))))
        
        df1 <- counts_data %>%
            left_join(results[[1]] %>% transmute(gene_id, filtered = 'no')) %>%
            replace_na(list(filtered = 'yes'))
        
        df2 <- df1[,-c(1)]
        output$plot5 <- renderPlot(plot_variance_vs_median(df2,TRUE,"Median count vs variance"))
        
        output$plot6 <- renderPlot(plot_zeros_vs_median(df2,"Median count vs number of zeros"))
        
        output$heatmap <- renderPlot({
            counts_filtered <- results[[1]]
            df3 <- log10(counts_filtered[,-c(1)])
            df3[df3=="-Inf"] <- 0
            heatmap3(as.matrix(df3),main="\r\nHeatmap of counts remaining after filtering\r\n")
            }
        )
        
        output$pca_plot <- renderPlot({
            metadata1 <- meta_info_from_labels(colnames(counts_data)[colnames(counts_data) != "gene_id"])
            plot_pca(counts_data[c(-1)],metadata1,input$pc1,input$pc2,paste0("PCA plot of PC",input$pc1," vs PC",input$pc2))
        })
    })
    observeEvent(input$submit3,{
        de_results <- load_de_results()
        volcano <- volcano_plot(de_results,input$xvar,input$yvar,input$slider3,input$col1,input$col2)
        table <- draw_table(de_results,input$slider3)
        output$volcano <- renderPlot(volcano,height = 600)
        output$de_table <-  DT::renderDataTable(DT::datatable(table, 
                           extensions = 'Buttons', class = "display",
                           options = list(paging = TRUE, searching = TRUE,
                                          fixedColumns = TRUE, autoWidth = TRUE,
                                          ordering = TRUE, dom = 'Bfrtip',
                                          buttons = c('copy', 'csv'))))
    })
    observeEvent(input$submit4,{
        counts_data <- load_counts_data1()
        dataf<-load_metadata1()
        first_column <- as.list(counts_data[,1])
        updateSelectizeInput(session, 'gene', choices = first_column, server = TRUE)
        observeEvent(input$submit5,{
            result <- individual_gene_expression(counts_data,dataf,input$gene)
            categorical_input <- input$categorical_input
            coul <- brewer.pal(5, "Set2") 
            if(input$plot_type=="Bar plot") {
                output$de <- renderPlot(barplot(table(result[,categorical_input]),
                                                main=paste0("Bar chart for ",categorical_input),
                                                ylab="Count",col=coul))
            }
            if(input$plot_type=="Beeswarm plot") {
                plot <- ggplot(result, aes(x=eval(as.name(categorical_input)),y=counts,color=eval(as.name(categorical_input)))) +
                    geom_beeswarm()+
                    ggtitle(paste0("Beeswarm plot of ",categorical_input))+
                    labs(color=categorical_input,x=categorical_input,y='Counts') 
                output$de <- renderPlot(plot)
            }
            if(input$plot_type=="Boxplot") {
                output$de <- renderPlot(ggplot(result, aes(x=eval(as.name(categorical_input)),y=counts,fill=eval(as.name(categorical_input)))) +
                                            geom_boxplot()+
                                            ggtitle(paste0("Boxplot of ",categorical_input))+ labs(fill=categorical_input,x=categorical_input,y='Counts') )
            }
            if(input$plot_type=="Violin plot") {
                output$de <- renderPlot({g<-ggplot(result, aes(x=eval(as.name(categorical_input)),y=counts,fill=eval(as.name(categorical_input)))) +
                                            geom_violin()+
                    ggtitle(paste0("Violin plot of ",categorical_input))+labs(fill=categorical_input,x=categorical_input,y='Counts') 
                g})
            }
        })
    })
}

# Run the application
shinyApp(ui = ui, server = server)
