
#' @title Look for Overlap between Two Groups
#' @description Will generate an overlap analysis on a pairwise basis of multiple datasets, checking for significance with Fischers Exact Test, as well as outputting plots to help visualize the overlaps.
#' @param results The filepath to OR a dataframe containing the following columns:
#' \itemize{
#'  \item{Protein : }{a unique identifier for the protein, such as Uniprot Accession. The annotation feature only works with Uniprot Accession codes.}
#'  \item{dataset : }{The name of the group/dataset/bait associated with this protein. (interaction)}
#'  \item{thresh : }{Whether the protein is to be considered a confident interactor in this dataset.}
#'  \item{plot.venn : }{Whether or not to plot and save the venn diagrams of the overlaps. May not want to if there are a ton of comparisons.}
#' }
#' @param out_file A file path to where the files will be saved.
#' @param annotate Whether to add in human annotations to the final results or not. (Default = TRUE)
#' @param plot.venn Whether to plot the pairwise Venn Diagrams of each group (default = TRUE). This will result in many plots if there are a lot of comparisons. This is recommended for a smaller number of groups in the dataset.
#' @keywords overlap, overlap analysis
#' interaction_overlap()
#' @export
interaction_overlap <- function(results, out_file, annotate = T, plot.venn = T) {
    
    cat("READING IN DATA...")
    results = as.data.frame(checkIfFile(results, is.evidence = FALSE), stringsAsFactors = F)
    cat("DONE\n")
    
    # create directory if it doesn't exists
    if (!dir.exists(dirname(out_file))) 
        dir.create(dirname(out_file), recursive = T)
    
    x = data.frame(t(combn(unique(results$dataset), 2)), stringsAsFactors = F)  # get all combos of data sets
    # remove negatives from this
    x <- x[!unlist(apply(x, 1, function(y) {
        return(any(grepl("MOCK", y)))
    })), ]
    
    cat("CALCULATING SIGNIFICANCE OF OVERLAP...")
    overlap.table = c()
    fishers.results = c()
    intersections = c()
    for (i in 1:dim(x)[1]) {
        # simplify the data breakdown
        dat1.all = unique(results[results$dataset == x[i, 1], "Protein"])
        dat2.all = unique(results[results$dataset == x[i, 2], "Protein"])
        dat1.thresh = unique(results[(results$dataset == x[i, 1]) & (results$thresh == 1), "Protein"])
        dat2.thresh = unique(results[(results$dataset == x[i, 2]) & (results$thresh == 1), "Protein"])
        
        # get overlapping proteins
        m11 <- length(intersect(dat1.thresh, dat2.thresh))
        # get significant dat1 proteins not in dat2 hits
        m12 <- length(setdiff(dat1.thresh, dat2.thresh))
        # get significant dat2 proteins not in dat1 hits
        m21 <- length(setdiff(dat2.thresh, dat1.thresh))
        # get all the proteins not making threshold
        m22 <- length(setdiff(unique(c(dat1.all, dat2.all)), unique(c(dat1.thresh, dat2.thresh))))
        
        # get names for matrix. Only really for pretty print
        mat_names = list()
        mat_names[[x[i, 1]]] = c("In", "Not In")
        mat_names[[x[i, 2]]] = c("In", "Not In")
        
        dat1_v_dat2 <- matrix(c(m11, m12, m21, m22), nrow = 2, dimnames = mat_names)
        
        # calculate p-value using 2 tails
        pval.dat1_v_dat2 <- fisher.test(dat1_v_dat2, alternative = "greater")
        
        i.results <- data.frame(dataset1 = x[i, 1], dataset2 = x[i, 2], pval = pval.dat1_v_dat2$p.value, overlap = m11, dat1_not_dat2 = m12, 
            dat2_not_dat1 = m21, dataset1.thresh = length(dat1.thresh), dataset2.thresh = length(dat2.thresh), dataset1.size = length(dat1.all), 
            dataset2.size = length(dat2.all))
        # backup results
        fishers.results = rbind(fishers.results, i.results)
        
        # backup intersecting proteins
        if (m11 > 0) {
            i.intersections = data.frame(dataset1 = x[i, 1], dataset2 = x[i, 2], Protein = intersect(dat1.thresh, dat2.thresh), 
                pval = pval.dat1_v_dat2$p.value, stringsAsFactors = F)
            intersections = rbind(intersections, i.intersections)
        }
        
        
        if (plot.venn) {
            # plot venn diagram of the two datasets <<------- turning off since so many plots will be made
            # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            venn_out <- gsub(basename(out_file), paste0(x[i, 1], "_vs_", x[i, 2], ".pdf"), out_file)
            pdf(venn_out)
            venn.plot = plotVenn2(m11, m12, m22, pval = NULL)
            grid.draw(venn.plot)
            dev.off()
            grid.newpage()
            
            ## Deprecated tmp = results[which( (results$dataset %in% x[i,]) & (results$thresh==1) ),] venn_out <- gsub(basename(out_file),
            ## paste0(x[i,1],'_vs_', x[i,2],'.pdf'), out_file) pdf(venn_out) plotVenn(tmp, pval=pval.dat1_v_dat2$p.value) dev.off()
        }
        
        
        # update overlap table from results
        overlap.table <- rbind(overlap.table, data.frame(dataset = paste(x[i, 1], x[i, 2], sep = "-"), overlap = m11, set1only = m12, 
            set2only = m21, p.overlap = m11/sum(m11, m12, m21), p.set1only = m12/sum(m11, m12, m21), p.set2only = m21/sum(m11, 
                m12, m21), stringsAsFactors = F))
    }
    cat("DONE\n")
    
    if (sum(overlap.table$overlap > 0) > 1) {
        # plot combined overlap plot with all the sets (with overlap) together
        proportionPlot(overlap.table[overlap.table$overlap > 0, ], out_file)
        
        # annotate intersections
        if (annotate == 1) {
            keys <- suppressMessages(select(org.Hs.eg.db, keys = unique(as.character(intersections$Protein)), columns = c("ENTREZID", 
                "GENENAME", "SYMBOL", "UNIPROT"), keytype = "UNIPROT"))
            intersections <- merge(intersections, keys, by.x = "Protein", by.y = "UNIPROT", all.x = T)
        }
    }
    
    # write out intersections
    write.table(intersections, gsub(".txt", "_intersections.txt", out_file), quote = F, row.names = F, sep = "\t")
    # write out fishers results
    fishers.results <- fishers.results[order(fishers.results$pval, decreasing = F), ]
    write.table(fishers.results, out_file, quote = F, row.names = F, sep = "\t")
    return(fishers.results)
}


# create a venn diagram of a group of datasets
#' @export
plotVenn <- function(x, pval = NULL) {
    x.wide = dcast(x, Protein ~ dataset)
    idx <- (x.wide[-1] == 0) | is.na(x.wide[-1])
    x.wide[, -1][idx] = 0
    
    # remove proteins that arent ever confident hits
    if (length(which(apply(x.wide[, -1], 1, sum) == 0)) > 0) {
        x.wide <- x.wide[-which(apply(x.wide[, -1], 1, sum) == 0), ]
    }
    
    if (exists("pval") & !is.null(pval)) {
        if (round(pval, 4) == 0) {
            pval = format(pval, scientific = TRUE)
            vennDiagram(x.wide[, -1], main = paste0("Fisher's Exact P-value : ", pval))
        } else {
            vennDiagram(x.wide[, -1], main = paste0("Fisher's Exact P-value : ", round(pval, 4)))
        }
    } else {
        vennDiagram(x.wide[, -1])
    }
}


# This will create and save a venn digram that displays both the number of members ber group, but also the percentage of the
# whole group
#' @export
plotVenn2 <- function(m11, m12, m22, pval = NULL) {
    # Format the title
    if (exists("pval") & !is.null(pval)) {
        if (round(pval, 4) == 0) {
            pval = format(pval, scientific = TRUE)
            main_title = paste0("Fisher's Exact P-value : ", pval)
        } else {
            main_title = paste0("Fisher's Exact P-value : ", round(pval, 4))
        }
    } else {
        main_title = paste0(x[i, 1], "_vs_", x[i, 2])
    }
    
    # Plot the pairwise venn diagram
    venn.plot = draw.pairwise.venn(area1 = m12 + m11, area2 = m21 + m11, cross.area = m11, c(x[i, 1], x[i, 2]), print.mode = c("raw", 
        "percent"), cex = 2, cat.cex = 2, cat.pos = c(200, 160), main = main_title)
    return(venn.plot)
}


# Make a horizontal protein plot of the proportion
#' @export
proportionPlot <- function(x, out_file) {
    # make long version of proportion matrix
    idx <- grep("p\\.", names(x))
    tmp1 <- melt(x[, c(1, idx)], id = "dataset")
    tmp1$variable <- gsub("p\\.", "", tmp1$variable)
    # make long version of the count matrix
    tmp2 <- melt(x[, c(-idx)], id = "dataset")
    # bring the counts and percentages together and make better names
    x.long = merge(tmp1, tmp2, by = c("dataset", "variable"))
    names(x.long)[3:4] = c("percentage", "count")
    
    # overlap.table.long <- melt(overlap.table, id=c('dataset'))
    x.long$variable <- factor(x.long$variable, levels = c("set2only", "overlap", "set1only"))
    # re-order by factors in order to calculate location of where to print text
    x.long = x.long[order(x.long$dataset, x.long$variable, decreasing = T), ]
    
    # adjust the percentages
    x.long$percentage = x.long$percentage * 100
    # add in a coordinate to print out the
    x.long <- ddply(x.long, .(dataset), transform, pos = cumsum(percentage) - (0.5 * percentage))
    x.long[, c("percentage", "pos")] = round(x.long[, c("percentage", "pos")], 3)
    
    # create a plot based on the PROPORTION of overlap.  Numbers on the plot are actual counts of proteins in that section
    p <- ggplot(x.long, aes(x = dataset, y = percentage)) + geom_bar(stat = "identity", aes(fill = variable)) + scale_fill_manual(values = c("red", 
        "purple", "royalblue")) + labs(y = "Percentage", x = "Compared Datasets", fill = "Dataset", title = "Number of Genes Found") + 
        geom_text(data = x.long, aes(x = dataset, y = pos, label = count), size = 4) + theme(plot.title = element_text(hjust = 0.5)) + 
        coord_flip()
    
    # Save the plot
    ggsave(filename = gsub(".txt", "_count_all.pdf", out_file), plot = p, width = 7, height = max(7, (2 * 7 * length(unique(x.long$dataset))/100)), 
        limitsize = FALSE)
    p
    
    
    # create a plot based on the PERCENTAGE of overlap.  Numbers on the plot are actual PERCENTAGE of proteins in that section
    p <- ggplot(x.long, aes(x = dataset, y = percentage)) + geom_bar(stat = "identity", aes(fill = variable)) + scale_fill_manual(values = c("red", 
        "purple", "royalblue")) + labs(y = "Percentage", x = "Compared Datasets", fill = "Dataset", title = "Number of Genes Found") + 
        geom_text(data = x.long, aes(x = dataset, y = pos, label = round(percentage, 1)), size = 4) + theme(plot.title = element_text(hjust = 0.5)) + 
        coord_flip()
    
    # Save the plot
    ggsave(filename = gsub(".txt", "_percentage_all.pdf", out_file), plot = p, width = 7, height = max(7, (2 * 7 * length(unique(x.long$dataset))/100)), 
        limitsize = FALSE)
    p
    
}


# plotting a heatmap from the correlations of each dataset
#' @export
corrPlot <- function(dat, out_file, plot_title = "") {
    # keep only things confident hits
    dat <- dat[which(dat$thresh > 0), ]
    dat$thresh <- NULL
    
    # put df into matrix format for the cor function
    dat.wide <- suppressMessages(dcast(data = dat, Protein ~ dataset))
    # add in rwo names
    row.names(dat.wide) = paste(dat.wide$Gene, dat.wide$Protein, sep = "-")
    dat.wide$Gene <- dat.wide$Protein <- NULL
    dat.wide[is.na(dat.wide)] <- 0
    dat.wide[dat.wide != 0] <- 1  #~!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Should we care if the gene shows up multiple times for a dataset
    
    dat.mat = data.matrix(dat.wide)
    dat.cor = cor(dat.mat, use = "pairwise.complete.obs", method = "pearson")
    # calculate the size of the plot to be outputted so we can read the labels
    sz <- max(floor(dim(dat.cor)[1]/4), 7)
    
    # order by hierarchical clustering
    pdf(gsub(".txt", "_heatmap_hclust.pdf", out_file), width = sz, height = sz)
    # pdf('~/Desktop/heatmap_hclust.pdf', width = sz, height = sz)
    corrplot(dat.cor, method = "color", order = "hclust", bg = "grey", title = plot_title, mar = c(0, 0, 1, 0))
    dev.off()
}


# function to create a directory if it doesn't exist
#' @export
checkDir <- function(output_dir) {
    if (!dir.exists(dirname(output_dir))) 
        dir.create(output_dir, recursive = T)
    return(output_dir)
}







# INPUT: dataframe: Protein = unique identifier for the protein, such as Uniprot Accession dataset = The name of the
# group/dataset/bait associated with this protein. (interaction) thresh = Whether the protein is to be considered a confident
# interactor in this dataset.  plot.venn = Whether or not to plot and save the venn diagrams of the overlaps. May not want to
# if there are a ton of comparisons.  head(results) out_file = '~/Desktop/overlap_test/overlap.txt' fishersBatch(results,
# out_file, annotate=T, plot.venn=T)






