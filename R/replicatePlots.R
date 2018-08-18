# ------------------------------------------------------------------------------
#' @title Create the replicate plots based on the pairings from 
#' the replicate plot file
#' 
#' @description Outputs a replicate plots based on a user provied file 
#' containing the replicates to be compared. Values are based on the log2 
#' value of the maximum intensities per modified sequence. 
#' The "replicate file" should describe which replicates of which conditions 
#' should be compared against each other. Each row represents a replicate 
#' plot to be created. The file should be structured using the following 
#' format and column names:
#' 
#' | condition1 | rep1_1 | rep1_2 | condition2 | rep2_1 | rep2_2 |
#' |---|---|---|---|---|---|
#' |  Infected | Infected_Rep1_name | Infected_Rep2_name | Negative | Negative_Rep1_name | Negative_Rep2_name |
#' | etc... |   |   |   |   |   |
#' | etc... |   |   |   |   |   |
#' 
#' @param input_file MaxQuant evidence file and location
#' @param keys_file Keys file with the experimental details
#' @param replicate_file Replicate file
#' @param out_file Output text file with the intensity values of every feature
#' @return The output file of the summary of features and intensity values
#' @keywords evidence, replica, plots
#' artms_replicatePlots()
#' @export
artms_replicatePlots <- function(input_file, keys_file, replicate_file, out_file){
  cat(">> READING IN FILES...\n")
  # read in data
  dat <- read.delim(input_file, stringsAsFactors=F)
  # keys
  keys <- read.delim(keys_file, stringsAsFactors=F)
  # profile plot list
  repplot <- read.delim(replicate_file, stringsAsFactors=F)
  
  # remove negatives from MaxQuant
  if( length(grep("__", dat$Proteins)) >0 ) dat <- dat[-grep("__", dat$Proteins),]
  
  # remove blank protein names
  if(any(dat$Proteins == "")){ dat <- dat[-which(dat$Proteins == ""),]}
  
  # for UB, jj suggests using unique peptide and charge for distinguishing 
  # numbers
  # NOTE: dimensions betwen x and dat may differ if there is data in dat that 
  # isn't in the keys file
  names(dat)[grep("Raw.file", names(dat))] <- 'RawFile'
  # x <- merge(dat, keys[,c('RawFile','Condition','BioReplicate','IsotopeLabelType')], by=c('RawFile', 'IsotopeLabelType') )  ## !!!!!!! DIfferent for SILAC
  x <- merge(dat, keys[,c('RawFile','Condition','BioReplicate')], by=c('RawFile') )
  
  # Put into a data matrix format
  x <- dcast(data=x, Proteins+Modified.sequence+Charge~Condition+BioReplicate, value.var="Intensity", max, na.rm=T)
  # remove cases where -Inf  is introduced
  x[x==-Inf] = 0   ###### May cause problems? Check.
  write.table(x, out_file, quote=F, row.names=F, sep='\t')
  
  # cycle through the condition pairs in the file and plot each pair
  for(i in 1:dim(repplot)[1]){
    
    cat(">> PLOTTING REPLICATE PLOT ", i, "\n")
    
    # check if the replicate combination exists in the plots
    rep1_1 <- paste(repplot$condition1[i], repplot$rep1_1[i], sep="_")
    rep1_2 <- paste(repplot$condition1[i], repplot$rep1_2[i], sep="_")
    rep2_1 <- paste(repplot$condition2[i], repplot$rep2_1[i], sep="_")
    rep2_2 <- paste(repplot$condition2[i], repplot$rep2_2[i], sep="_")
    reps <- c(rep1_1, rep1_2, rep2_1, rep2_2)
    if( !any( !(reps %in% names(x)) ) ){
      
      # prep 1st replicate comparison for plot
      rep1 <- log2( x[,paste(repplot$condition1[i], repplot$rep1_1[i], sep="_")] / x[,paste(repplot$condition2[i], repplot$rep2_1[i], sep="_")] )
      # prep 2nd replicate comparison for plot
      rep2 <- log2( x[,paste(repplot$condition1[i], repplot$rep1_2[i], sep="_")] / x[,paste(repplot$condition2[i], repplot$rep2_2[i], sep="_")] )
      
      # remove NA pairs
      idx <- which( is.na(rep1) | is.na(rep2) | is.infinite(rep1) | is.infinite(rep2))
      rep1 <- rep1[-idx]
      rep2 <- rep2[-idx]
      reps.cor <- cor(rep1, rep2, use="pairwise.complete.obs", method="pearson")
      
      # set up a square plot centered at 0
      x.lim <- ceiling(max(abs(c(rep1, rep2 )), na.rm=T))
      y.lim <- c( -x.lim, x.lim)
      x.lim <- c( -x.lim, x.lim)
      
      # name axes labels
      y.label <- paste0(repplot$condition1[i]," vs. ",repplot$condition2[i], "  (", repplot$rep1_1[i],"/", repplot$rep2_1[i],")")
      x.label <- paste0(repplot$condition1[i]," vs. ",repplot$condition2[i], "  (", repplot$rep1_2[i],"/", repplot$rep2_2[i],")")
      # make plot name
      plot.name <- paste( repplot$condition1[i], " vs "  , repplot$condition2[i], " R = ",round(reps.cor,3), sep="" )
      plot.name2 <- paste( repplot$condition1[i], " vs "  , repplot$condition2[i], " R ",round(reps.cor,3), sep="" )
      #       pdf( paste( dirname(out_file), "/", gsub(" ","_",plot.name) ,"_", repplot$rep1_1[i], "_", repplot$rep1_2[i], ".pdf", sep="") )
      #       plot(rep1, rep2, main=plot.name, xlab=repplot$rep1_1[i], ylab=repplot$rep1_2[i], xlim=x.lim, ylim=y.lim, pch=".")
      #       dev.off()
      tmp <- data.frame(rep1, rep2, stringsAsFactors=F)
      p <- ggplot(tmp, aes(x=rep1, y=rep2)) + 
        geom_point() +
        xlim(x.lim[1],x.lim[2]) + 
        ylim(x.lim[1],x.lim[2]) + 
        ggtitle(plot.name) + 
        labs(x=x.label, y=y.label)
      ggsave(filename = paste( dirname(out_file), "/", gsub(" ","_",plot.name2) ,"_", repplot$rep1_1[i], "_", repplot$rep1_2[i], ".pdf", sep=""), 
             plot=p, 
             width = 10, 
             height = 10)
    }else{
      warning("REPLICATE PLOT ",i," NOT MADE -- MISSING DATA FROM \n", paste("\t", reps[!(reps %in% names(x))],"\n", collapse=""))
    }
  }
}




