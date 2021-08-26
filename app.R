library(shiny)
library(shinyFiles)
library(shinydashboard)
#library(xlsx)
library(survival)
library(survminer)
library(plotly)
library(shinyjs)

jsCode <- "shinyjs.resetSel = function() { Plotly.restyle(plot, {selectedpoints: [null]});}"


#setwd("/srv/shiny-server")
anno     <- read.csv("LGG_umap_2-PCA5-10k.txt", sep="\t", row.names=2, header=T);
#anno     <- read.csv("umap_2shiny.txt", sep="\t", row.names=2, header=T);
mycolors <- readRDS("mycols.rds")

## Classification uniformly match colors
anno$CNS.MCF      <- gsub(" |,|/|__", "_", anno$CNS.MCF) 
anno$CNS.MCF      <- gsub("__", "_", anno$CNS.MCF) 
anno$CNS.MCF      <- gsub("__", "_", anno$CNS.MCF) 
anno$CNS.Subclass <- gsub(" |,|/|__", "_", anno$CNS.Subclass ) 
anno$CNS.Subclass <- gsub("__", "_", anno$CNS.Subclass)
anno$CNS.Subclass <- gsub("__", "_", anno$CNS.Subclass)
anno$CNS.Subclass <- gsub("_$", "", anno$CNS.Subclass, perl=TRUE)


options(max.print=1000000)
ui <- dashboardPage(skin = "red",
                    dashboardHeader(
                      # Set height of dashboardHeader
                      tags$li(class = "dropdown",
                              tags$style(".main-header {max-height: 70px}"),
                              tags$style(".main-header .logo {height: 70px;}"),
                              tags$style(".sidebar-toggle {height: 70px; padding-top: 0px !important;}"),
                              tags$style(".navbar {min-height:70px !important}"),
                      ),
                       title = "National Cancer Institute"),
                    ## Sidebar content
                    dashboardSidebar(
                      sidebarMenu(
                        menuItem("UMAP",  tabName = "umap"),
                        menuItem("densMAP",     tabName = "umapd3"),
                        #menuItem("densMAP", tabName = "densMAP"),
                        #menuItem("t-SNE", tabName = "tsne"),
                        #menuItem("purity correction test", tabName = "correction"),
                        # menuItem("Copy number", tabName = "copy"),
                        # menuItem("MGMT promoter", tabName = "MGMT"),
                        #menuItem("DKFZ MNP Classifier", 
                        #         href = "https://www.molecularneuropathology.org/mnp"),
                        radioButtons("radio", h3("Choose groups", style = "font-size:17px;"),
                                     choices = list("Group 1" = 1, "Group 2" = 2,
                                                    "Group 3" = 3), selected = 1),
                        actionButton("action", "Reset all Groups"),
                        br(),
                        #uiOutput("currentSelections"),
                        #shinySaveButton("save", "Save file", "Save file as ...", filetype=list(xlsx="xlsx")),
                        tags$style(".topimg {
                            margin-left:10px;
                            position: absolute;
                            bottom:0;   
                          }"),
                        tags$style(".left-side, .main-sidebar {padding-top: 90px}"),
                        div(class="topimg",img(src = "/srv/shiny-server/www/nih-logo-footer.png", height = "30%", width = "30%"), style="text-align: center;")
                      )
                    ),
                    ## Body content
                    dashboardBody(
                      tags$head(tags$style(HTML('
                        .main-header .logo {
                          font-family: "Arial", Helvetica, sans-serif;
                          #font-weight: bold;
                          #font-size: 30px;
                          #float: left!important;
                          #line-height: 200px !important;
                          padding: 10px 0px;}
                        .content-wrapper, .right-side {
                                background-color: white;
                        }
                        .form-group{
                        margin-bottom: 0px;
                        }
                          .main-sidebar { font-size: 20px; }
                          .radio label {
                            font-size: 17px;
                          }
                      '))),
                      useShinyjs(),
                      extendShinyjs(text = jsCode, functions = c("resetSel")),
                      tabItems(
                        # First tab content
                        tabItem(tabName = "umap",
                                # fluidRow(
                                #   scatterD3Output("scatterPlot", width = "1000px", height = "500px"),
                                # )
                                fluidRow(
                                  column(12,
                                         fluidRow(
                                           column(width = 6, offset = 0, div(plotlyOutput("plot1"), align = "left"),
                                                  # fluidRow(id="fluidrow1",
                                                  #   column(12,
                                                  #          plotOutput("plot2",brush = brushOpts(id = "plot2_brush",resetOnNew = TRUE)))
                                                  #          #DT::dataTableOutput("brush"),
                                                  #          #tags$head(tags$style("#brush{font-size:11px;}"))) #end column
                                                  # )
                                           ),
                                           column(width = 5, offset = 1,
                                                  div(DT::dataTableOutput("brush"), align = "right"),
                                                  tags$head(tags$style("#brush{font-size:12px;}")),
                                                  fluidRow(
                                                    column(12, div(plotOutput("survivalCurve", width = "75%"), align = "bottom"))
                                                  )
                                                  )
                                                  #plotlyOutput("plot1")) #end column
                                         )
                                  )
                                )
                        ), #end tabItem
                       
                        
                          # Second tab content
                        tabItem(tabName = "umapd3",
                                fluidRow(
                                    plotlyOutput("scatterPlot_umap", width = "1600px", height = "1000px")
                                )
                        ) #end tabItem
                      )
                    )
)


server <- shinyServer(function(input, output, session) {
  
  ## data.frame with an extra group column (initially set to NA)  
  rv <- reactiveValues(data_df = anno %>% mutate(group = NA))
  
  ## when a selection is made, assign group values to data_df based on selected radio button
  observeEvent(
    event_data("plotly_selected"), {
      d <- event_data("plotly_selected")
      ## reset values for this group
      rv$data_df$group <- ifelse(rv$data_df$group == input$radio, NA, rv$data_df$group)
      ## then re-assign values:
      rv$data_df[d$key,"group"] <- input$radio
    }
  )
  
  ## when reset button is pressed, reset the selection rectangle 
  ## and also reset the group column of data_df to NA
  observeEvent(input$action, {
    js$resetSel()
    rv$data_df$group <- NA
  })
  
  ## when radio button changes, reset the selection rectangle and reset plotly_selected
  ## (otherwise selecting the same set of points for two groups consecutively will 
  ## not register the selection the second time)
  observeEvent(input$radio, {
    js$resetSel()
    runjs("Shiny.setInputValue('plotly_selected-A', null);")
  })
  
  # Linked plots (middle and right)
  ranges2 <- reactiveValues(x = NULL, y = NULL)
  
   
  ##  draw the first main plot
  output$plot1 <- renderPlotly({
    # use the key aesthetic/argument to help uniquely identify selected observations
    key1 <- row.names(anno) 
    key2 <- anno$Sample
#    cnames <- aggregate(cbind(x1, y1) ~ CNS.Subclass, data=anno, FUN=function(x)mean((x)))    
    p <- ggplot(data=anno, aes(x=x1,y=y1, key=key1)) +
       geom_point(aes(color = CNS.Subclass),  size=2, alpha=1) +
#      geom_text(data=cnames, aes(x=x1, y=y1, label = CNS.Subclass), position = position_dodge(width=0.5),  size=2.5, inherit.aes = FALSE) +       ## Text for clusters
       theme_classic() +
       theme(axis.text.y = element_text(size=9, color="black"),
            axis.text.x = element_text(size=9, color="black"),
            axis.ticks.y = element_line(color="black", size = 0.5),
            axis.ticks.x = element_line(color="black", size = 0.5),
            axis.ticks.length = unit(2,"mm"),
            panel.border = element_rect(colour = "black", fill=NA, size=1),
            #panel.border = element_blank(),
            panel.grid.major = element_line(colour="grey", size=0.5),
            axis.line = element_blank(),
            legend.text=element_text(size=10),
            legend.title=element_text(size=11),
            legend.key.size = unit(0.5, 'lines'),
            axis.title.x = element_blank(),
            axis.title.y = element_blank()) +
#     labs(x = "umap1", 
#           y = "umap2") +
       theme(legend.position="right") +
       theme(plot.margin=unit(c(0,0,0,0),"cm")) +
       scale_color_manual("Methylation class", values = mycolors[sort(names(mycolors))], guide = guide_legend(override.aes = list(shape = 15))) +
       scale_x_continuous(breaks = seq(-100, 100, by=5)) +
       scale_y_continuous(breaks = seq(-100, 100, by=5)) +
       coord_fixed(xlim = ranges2$x, ylim = ranges2$y, expand = FALSE)

    ggplotly(p) %>% layout(height = 1000, width = 1000, dragmode = "select",
                  xaxis = list(autorange = TRUE), yaxis = list(autorange = TRUE)) %>%
    
    add_annotations(x = subset(p$data, Study == "compass")$x1,
                      y = subset(p$data, Study == "compass")$y1,
                      text = subset(p$data, Study == "compass")$Sample,
                      showarrow = TRUE,
                      arrowcolor='black',
                      arrowhead = 6,
                      arrowsize = 1,
                      xref = "x",
                      yref = "y",
                      font = list(color = 'black',
                                  family = 'arial',
                                  size = 16))
    

  })  
  
  ## for each group, show the number of selected points
  ## (not required by the rest of the app but useful for debugging)
  output$currentSelections <- renderUI({
    number_by_class <- summary(factor(rv$data_df$group, levels = c("1","2","3")))
    tagList(
      h5("Current Selections:"),
      p(paste0("Group 1: ",number_by_class[1], " points selected")),
      p(paste0("Group 2: ",number_by_class[2], " points selected")),
      p(paste0("Group 3: ",number_by_class[3], " points selected"))
    )
  })

  output$brush <- DT::renderDataTable({
    d <- event_data("plotly_selected")
    req(d)
      DT::datatable(anno[unlist(d$key), c(
        "Sample",
        "Age",
        "Sex", 
        "material_prediction",
        "RFpurity.ABSOLUTE",
#      "LUMP",
        "Study",
        "Location_general",  
        "Location_specific",  
        "Histology",
        "Molecular", 
        "CNS.MCF",
        "CNS.MCF.score", 
        "CNS.Subclass",
        "CNS.Subclass.score")],
				options = list(lengthMenu = c(5, 30, 50), pageLength = 10, scrollY = '450px')) %>%
      DT::formatStyle(columns = c(1, 2, 3, 4, 5, 6, 7), fontSize = '100%')
  })
  
  # # When a double-click happens, check if there's a brush on the plot.
  # # If so, zoom to the brush bounds; if not, reset the zoom.
  # observe({
  #   brush <- input$plot2_brush
  #   if (!is.null(brush)) {
  #     ranges2$x <- c(brush$xmin, brush$xmax)
  #     ranges2$y <- c(brush$ymin, brush$ymax)
  #     
  #   } else {
  #     ranges2$x <- NULL
  #     ranges2$y <- NULL
  #   }
  # })
  
  ## draw survival curves if a point has been selected
  ## if none have been selected then draw a blank plot with matching background color
  output$survivalCurve <- renderPlot({
    if (any(c(1,2,3) %in% rv$data_df$group)) {
      fit <- survfit(Surv(OS_months, OS_status) ~ group,
                     data = rv$data_df)
      ggsurvplot(fit, data = rv$data_df, 
                 censor.size=1, size = 1, 
                 conf.int = FALSE,
                 conf.int.style = "ribbon",
                 risk.table = FALSE, 
                 pval = TRUE,
                 pval.size = 10,
                 #palette = "jco",
                 legend = "none",
                 legend.title = "",
                 xlab = "Overall survival (months)",
                 palette = c("red", "blue", "green"),
                 ggtheme = theme_bw(base_size = 20,
                                    base_family = "Arial"))
    } else {
      par(bg = "#ecf0f5")
      plot.new()
    }
  }, height = 600, width = 800)
  
  ## draw the second tab plot 

    output$scatterPlot_umap <- renderPlotly({
    # use the key aesthetic/argument to help uniquely identify selected observations
    key2 <- anno$Sample
    
    p1 <- plot_ly(data   = anno, x = ~x2, y = ~y2,
                  type   = "scatter", mode = "markers",
                  color  = ~anno$CNS.Subclass,
                  colors = ~mycolors,
                  text   = ~paste(rownames(anno), 
                                 "<br>Sample: ", Sample, " Sex: ", Sex, " Age: ", Age, 
                                 "<br>DKFZ Subclass:",  CNS.Subclass, " | ",CNS.MCF, 
                                 "<br>", Location_general, Location_specific, Histology, Molecular)) %>% 
#      add_trace(CNS.Subclass, mode = "markers")      %>%
          layout(yaxis = list(title = ""), xaxis = list(title = ""), showlegend = T )    
    
    ggplotly(p1) %>% layout(height = 1000, width = 1300, dragmode = "select",
                     xaxis = list(autorange = TRUE),yaxis = list(autorange = TRUE))   %>% 
     
    add_annotations(x = anno[anno$Study == "compass", "x2"],
                        y = anno[anno$Study == "compass", "y2"],
                        text = anno[anno$Study == "compass", "Sample"],
                        showarrow = TRUE,
                        arrowcolor='black',
                        arrowhead = 6,
                        arrowsize = 1,
                        xref = "x",
                        yref = "y",
                        font = list(color = 'black',
                                     family = 'arial',
                                     size = 16))
    
  })  
})

shinyApp(ui = ui, server = server)
