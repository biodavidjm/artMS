# ------------------------------------------------------------------------------
#' @title Select Reference protein(s) for normalization
#'
#' NEEDS TO BE FIXED
#'
#' @description This function creates an interactive environment where the
#' user can select which prein(s) to use as a reference when normalizing
#' the different samples prior to analysis. Only proteins that are present
#' in each sample are considered. The user can remove the unnecessary/unwanted
#' proteins by checking the box next to the protein name. One can use the
#' general view of the proteins plotted across the samples, or at a clustering
#'   of the proteins.
#' @param dat_file (char) The filepath to the MaxQuant searched data (evidence)
#' file (txt tab delimited file).
#' @param keys_file (char) The filepath to the MSStats formatted keys file
#' (txt tab delimited file).
#' @return (shiny) app to choose a protein to be use for normalization
#' @keywords internal, msstat, reference protein, normalization
.select_ref <- function(keys_file, dat_file) {
  cat("Reading in files...\n")
  keys <- read.delim(keys_file, stringsAsFactors = FALSE)
  dat <- x <- read.delim(dat_file, stringsAsFactors = FALSE)
  
  cat("Removing MaxQuant described contaminants...\n")
  x <- x[-grep("__|;", x$Proteins),]
  if (length(which(x$Proteins == "")) > 0)
    x <- x[-which(x$Proteins == ""),]
  
  # Keep only the files in this experiment showing up in the keys file
  x <- unique(x[, c("Raw.file", "Proteins")])
  # rename the Rawfile name in the keys columsn
  names(keys)[grep("(r|R)aw.*(f|F)ile", names(keys))] <- "RawFile"
  x <- x[x$Raw.file %in% keys$RawFile,]
  
  # Counting how often a protein is found
  length(unique(x$Raw.file))
  tmp <- data.table::dcast(x, Proteins ~ 1)
  
  tmp <- tmp[order(tmp$`1`, decreasing = TRUE),]
  tmp <- tmp[tmp$`1` == length(unique(x$Raw.file)), "Proteins"]
  
  # keep only the proteins that appear in ALL the samples
  temp <-
    dat[(dat$Proteins %in% tmp) & (dat$Raw.file %in% keys$RawFile),]
  
  # getting into ggplot compatible format
  tmp <-
    data.table::dcast(temp,
                      Proteins ~ Raw.file,
                      value.var = "Intensity",
                      median,
                      na.rm = TRUE)
  tmp <- melt(tmp)
  tmp <- tmp[order(tmp$variable, tmp$value),]
  cat("Finished Pre-processing steps\n")
  
  
  # ~~~~~~~~~~~~~~~~~~~~~~~ BEGIN CLUSTERING ~~~~~~~~~~~~~~~~~~
  
  # want to cluster on intensity profile for each sample
  cat("Parsing data\n")
  x <- dcast(tmp, Proteins ~ variable, value.var = "value")
  row.names(x) <- x[, 1]
  x <- x[,-1]
  x <- data.matrix(x)
  
  cat("Removing NA's\n")
  x2 <- apply(x, 1, sum)
  if (length(which(is.na(x2))) > 0)
    x <- x[-which(is.na(x2)),]
  
  cat("Clustering data...\n")
  kmods <- kmeans(x = x, centers = 10)
  
  # apply clusters as colors and plot the results
  tmp2 <- data.frame((kmods$cluster))
  tmp2 <- cbind(tmp2, row.names(tmp2))
  names(tmp2) <- c("cluster", "Proteins")
  
  x.clust <- merge(tmp, tmp2, by = "Proteins")
  
  ## annotate the data a little better for legibility in the plots add in gene names
  prots <- unique(dat[, c("Proteins", "Gene.names")])
  x.clust <- merge(x.clust, prots, by = "Proteins")
  x.clust$prot_names <-
    paste(x.clust$Gene.names, x.clust$Proteins, sep = "-")
  
  
  # add in bait names
  x.clust <-
    merge(x.clust, keys[, c("RawFile", "Condition")], by.x = "variable", by.y = "RawFile")
  x.clust$sample_name <-
    paste(x.clust$variable, x.clust$Condition, sep = "-")
  
  
  
  cat("Finished CLustering!\n")
  
  ################################### FINISHED PREPROCESSING #####################
  
  
  
  ## Active shiny scripts ~~~~~~~~~~~~~~~~~~~~~~~~
  serv <- shiny::shinyServer(function(input, output) {
    # output dynamic range
    output$text1 <- shiny::renderText({
      idx <-
        intersect(
          grep(
            paste(input$check_proteins, collapse = "|"),
            x.clust$prot_names
          ),
          which(x.clust$cluster %in%
                  input$check_clusters)
        )
      paste("You have chosen a the following proteins:  \n",
            paste(unique(x.clust$Proteins[idx]), collapse = ","))
    })
    
    # output the dynamically filtered data
    output$plot1 <- plotly::renderPlotly({
      idx <-
        intersect(
          grep(
            paste(input$check_proteins, collapse = "|"),
            x.clust$prot_names
          ),
          which(x.clust$cluster %in%
                  input$check_clusters)
        )
      p <-
        ggplot(data = x.clust[idx,], aes(
          x = sample_name,
          y = value,
          colour = prot_names
        )) + geom_line(aes(group = prot_names)) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
      # p <- ggplot(data=tmp[grep(paste(input$check_proteins,collapse='|'), tmp$Proteins),], aes(x=variable, y=(value),
      # colour=Proteins )) + geom_line(aes(group=Proteins)) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
      plotly::ggplotly(p)
    })
    
    output$plot_clusts <- plotly::renderPlotly({
      # p <- ggplot(data=tmp[grep(paste(input$check_proteins,collapse='|'), tmp$Proteins),], aes(x=variable, y=(value),
      # colour=Proteins )) + geom_line(aes(group=Proteins)) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
      idx <-
        intersect(
          grep(
            paste(input$check_proteins, collapse = "|"),
            x.clust$prot_names
          ),
          which(x.clust$cluster %in%
                  input$check_clusters)
        )
      p <-
        ggplot(data = x.clust[idx,], aes(
          x = sample_name,
          y = value,
          colour = cluster
        )) + geom_line(aes(group = prot_names)) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
      plotly::ggplotly(p)
    })
  })
  
  
  
  
  #################################################
  
  
  ui <-
    shiny::shinyUI(shiny::basicPage(
      shiny::titlePanel("Reference Protein Selection"),
      shiny::sidebarLayout(
        shiny::div(style = "width: 1000px",
                   shiny::sidebarPanel(
                     shiny::checkboxGroupInput(
                       "check_proteins",
                       "Proteins displayed:",
                       sort(unique(x.clust$prot_names)),
                       selected = unique(x.clust$prot_names)
                     )
                   )),
        shiny::div(align = "left", shiny::mainPanel(
          shiny::tabsetPanel(
            shiny::tabPanel("Main",
                            plotly::plotlyOutput("plot1", height = "800px")),
            shiny::tabPanel(
              "Clusters",
              shiny::checkboxGroupInput(
                "check_clusters",
                "Proteins displayed:",
                sort(unique(x.clust$cluster)),
                selected = unique(x.clust$cluster),
                inline = TRUE
              ),
              plotly::plotlyOutput("plot_clusts",
                                   height = "800px")
            ),
            shiny::tabPanel("Proteins", shiny::textOutput("text1"))
          )
        ))
      )
    ))
  shiny::shinyApp(ui = ui, server = serv)
}
