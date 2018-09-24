
# ------------------------------------------------------------------------------
#' @title Quality Control of the MaxQuant summary.txt file
#' 
#' @description Performs quality control based on the information available in
#' the MaxQuant summary.txt file.
#' @param summary_file (char or data.frame) The evidence file path and name, or
#' data.frame
#' @param keys_file (char or data.frame) The keys file path and name or 
#' data.frame
#' @param isFractions (logical) `TRUE` if it is a 2D experiment (fractions). 
#' Default: `FALSE`
#' @param plotMS1SCANS (logical) It plots whatever
#' @param plotMS2 (logical) It plots whatever
#' @param plotMSMS (logical) It plots whatever
#' @param plotISOTOPE (logical) It plots whatever
#' @return A number of plots from the summary file
#' @keywords qc, summary, keys
#' @examples
#' # Testing warning if files are not submitted
#' test <- artms_qualityControlSummaryExtended(summary_file = NULL, 
#' keys_file = NULL)
#' @export
artms_qualityControlSummaryExtended <- function(summary_file, 
                                                keys_file,
                                                isFractions = FALSE,
                                                plotMS1SCANS = TRUE,
                                                plotMS2 = TRUE,
                                                plotMSMS = TRUE,
                                                plotISOTOPE = TRUE){
  
  cat("EXTENDED QUALITY CONTROL ANALYSIS (summary.txt based)---------------\n")
  
  if(is.null(summary_file) & is.null(keys_file)){
    return("You need to provide both evidence and keys")
  }
  
  # Getting data ready
  summarykeys <- artms_mergeEvidenceAndKeys(summary_file, keys_file, 
                                            isSummary = TRUE)  
  colnames(summarykeys) <- tolower(colnames(summarykeys))
  
  if("fraction" %in% colnames(summarykeys)){
    summary2fx <- summarykeys
  }
  
  if("fraction" %in% colnames(summarykeys)){
    summarykeys <- data.table(summarykeys[,c("condition", "bioreplicate", "ms", "ms.ms", "ms.ms.identified....", "isotope.patterns", "isotope.patterns.sequenced..z.1.")])
    summarykeys <- summarykeys[,list(ms=sum(ms), ms.ms=sum(ms.ms), ms.ms.identified....=mean(ms.ms.identified....), isotope.patterns=sum(isotope.patterns), isotope.patterns.sequenced..z.1.=sum(isotope.patterns.sequenced..z.1.)),
                               by=list(condition, bioreplicate)]
  }else{
    summarykeys <- data.table(summarykeys[,c("condition", "bioreplicate", "ms", "ms.ms", "ms.ms.identified....", "isotope.patterns", "isotope.patterns.sequenced..z.1.")])
  }
  
  summary2 <- plyr::ddply(summarykeys, c("condition"), 
        plyr::summarise, 
        num.MS1.mean = mean(ms), num.MS1.max = max(ms), 
        num.MS1.min = min(ms), 
        num.MS1.sem = sd(ms)/sqrt(length(ms)),
        num.MS2.mean = mean(ms.ms), 
        num.MS2.max = max(ms.ms), 
        num.MS2.min = min(ms.ms), 
        num.MS2.sem = sd(ms.ms)/sqrt(length(ms.ms)),
        num.IsotopePatterns.mean = mean(isotope.patterns), 
        num.IsotopePatterns.max = max(isotope.patterns),
        num.IsotopePatterns.min = min(isotope.patterns), 
        num.IsotopePatterns.sem = sd(isotope.patterns)/sqrt(length(isotope.patterns)),
        num.IsotopePatternsSeq.mean = mean(isotope.patterns.sequenced..z.1.), 
        num.IsotopePatternsSeq.max = max(isotope.patterns.sequenced..z.1.),
        num.IsotopePatternsSeq.min = min(isotope.patterns.sequenced..z.1.), 
        num.IsotopePatternsSeq.sem = sd(isotope.patterns.sequenced..z.1.)/sqrt(length(isotope.patterns.sequenced..z.1.)),
        pct.MS2Id.mean = mean(ms.ms.identified....), 
        pct.MS2Id.max = max(ms.ms.identified....),
        pct.MS2Id.min = min(ms.ms.identified....), 
        pct.MS2Id.sem = sd(ms.ms.identified....)/sqrt(length(ms.ms.identified....))
  ) 
  
  # PLOTS
  cat(">> GENERATING QC PLOTS\n")
  
  nsamples <- length(unique(summarykeys$bioreplicate))
  nconditions <- length(unique(summarykeys$condition))
  
  # Check the largest number so that width and height of pdf is not too large
  if(nsamples > 20){
    nsamples <- 20
  }
  
  if(nconditions > 7){
    nconditions <- 7
  }
  
  ## NUMBER OF MS1 SCANS
  if(plotMS1SCANS){
    cat("--- Plot NUMBER OF MS1 SCANS")
    pdf('QC_Plots_summary_MS1SCANS.pdf', width=nsamples*3, height=20, onefile = TRUE)
      aa <- ggplot(summarykeys, aes(x = bioreplicate, y = ms, fill = condition)) +
        geom_bar(stat="identity", alpha=0.7) +
        geom_text(aes(label=round(ms, digits=0)), vjust=1 , size = 15) +
        xlab("Experiment") + ylab("Counts") +
        ggtitle("Number of MS1 scans") +
        theme(legend.text = element_text(size=20)) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, size=20)) + 
        theme(axis.text.y = element_text(size=20)) + 
        theme(axis.title.x = element_text(size=30)) + 
        theme(axis.title.y = element_text(size=30)) + 
        theme(plot.title = element_text(size = 40)) + 
        scale_fill_brewer(palette="Spectral")
      print(aa)
      
      ab <- ggplot(summary2, aes(x = condition, y = num.MS1.mean, fill = factor(condition))) +
        geom_bar(stat="identity", position = position_dodge(width = 1), alpha=0.7) +
        geom_errorbar(aes(ymin=num.MS1.mean-num.MS1.sem, ymax=num.MS1.mean+num.MS1.sem), width=.2, position=position_dodge(.9)) +
        geom_text(aes(label=round(num.MS1.mean, digits=0)), hjust=0.5, vjust=-0.5, size = 10, position = position_dodge(width = 1)) +
        xlab("Condition") + ylab("Counts") +
        ggtitle("Mean number of MS1 scans per condition, error bar= std error of the mean") +
        theme(legend.text = element_text(size=20)) +
        theme(axis.text.x = element_text(angle = 0, hjust = 1, size=20)) + 
        theme(axis.text.y = element_text(size=20)) + 
        theme(axis.title.x = element_text(size=30)) + 
        theme(axis.title.y = element_text(size=30)) + 
        theme(plot.title = element_text(size = 40)) + 
        scale_fill_brewer(palette="Spectral")
      print(ab)
      
      if(isFractions){
        ac <- ggplot(summary2fx, aes(x = bioreplicate, y = ms, fill = factor(fraction))) +
          geom_bar(stat="identity", alpha=0.7) +
          geom_text(aes(label=round(ms, digits=0)), hjust=0.5, vjust=1.5, size = 7, position = position_stack()) +
          xlab("Experiment") + ylab("Counts") +
          ggtitle("Number of MS1 scans per Fraction") +
          theme(legend.text = element_text(size=20)) +
          theme(axis.text.x = element_text(angle = 90, hjust = 0, size=20)) + 
          theme(axis.text.y = element_text(size=20)) + 
          theme(axis.title.x = element_text(size=30)) + 
          theme(axis.title.y = element_text(size=30)) + 
          theme(plot.title = element_text(size = 40)) 
        print(ac)
      }
    garbage <- dev.off()
    cat(" done\n")
  }

  
  ## Number of MS2 scans
  if(plotMS2){
    cat("--- Plot Number of MS2 scans")
    pdf('QC_Plots_summary_MS2.pdf', width=nsamples*3, height=20, onefile = TRUE)
      ba <- ggplot(summarykeys, aes(x = bioreplicate, y = ms.ms, fill = condition)) +
        geom_bar(stat="identity", alpha=0.7) +
        geom_text(aes(label=round(ms.ms, digits=0)), vjust=1 , size = 15) +
        xlab("Experiment") + ylab("Counts") +
        ggtitle("Number of MS2 scans") +
        theme(legend.text = element_text(size=20)) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, size=20)) + 
        theme(axis.text.y = element_text(size=20)) + 
        theme(axis.title.x = element_text(size=30)) + 
        theme(axis.title.y = element_text(size=30)) + 
        theme(plot.title = element_text(size = 40)) + 
        scale_fill_brewer(palette="Spectral")
      print(ba)
      
      bb <- ggplot(summary2, aes(x = condition, y = num.MS2.mean, fill = factor(condition))) +
        geom_bar(stat="identity", position = position_dodge(width = 1), alpha=0.7) +
        geom_errorbar(aes(ymin=num.MS2.mean-num.MS2.sem, ymax=num.MS2.mean+num.MS2.sem), width=.2, position=position_dodge(.9)) +
        geom_text(aes(label=round(num.MS2.mean, digits=0)), hjust=0.5, vjust=-0.5, size = 10, position = position_dodge(width = 1)) +
        xlab("Condition") + ylab("Counts") +
        ggtitle("Mean number of MS2 scans per condition, error bar= std error of the mean") +
        theme(legend.text = element_text(size=20)) +
        theme(axis.text.x = element_text(angle = 0, hjust = 1, size=20)) + 
        theme(axis.text.y = element_text(size=20)) + 
        theme(axis.title.x = element_text(size=30)) + 
        theme(axis.title.y = element_text(size=30)) + 
        theme(plot.title = element_text(size = 40)) + 
        scale_fill_brewer(palette="Spectral")
      print(bb)
      
      if(isFractions){
        bc <- ggplot(summary2fx, aes(x = bioreplicate, y = ms.ms, fill = factor(fraction))) +
          geom_bar(stat="identity", alpha=0.7) +
          geom_text(aes(label=round(ms.ms, digits=0)), hjust=0.5, vjust=1.5, size = 7, position = position_stack()) +
          xlab("Experiment") + ylab("Counts") +
          ggtitle("Number of MS2 per Fraction") +
          theme(legend.text = element_text(size=20)) +
          theme(axis.text.x = element_text(angle = 90, hjust = 0, size=20)) + 
          theme(axis.text.y = element_text(size=20)) + 
          theme(axis.title.x = element_text(size=30)) + 
          theme(axis.title.y = element_text(size=30)) + 
          theme(plot.title = element_text(size = 40)) 
        print(bc)
      }
      
      summarykeys.scans <- melt(summarykeys[,1:4], id.vars=1:2)
      bd <- ggplot(summarykeys.scans, aes(x = interaction(variable, bioreplicate), y = value, fill = condition)) +
        geom_bar(stat="identity", alpha=0.7) +
        geom_text(aes(label=round(value, digits=0)), angle = 90, vjust=0.5 , size = 10) +
        xlab("Experiment") + ylab("Counts") +
        ggtitle("Number of MS1 and MS2 scans") +
        theme(legend.text = element_text(size=20)) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, size=20)) + 
        theme(axis.text.y = element_text(size=20)) + 
        theme(axis.title.x = element_text(size=30)) + 
        theme(axis.title.y = element_text(size=30)) + 
        theme(plot.title = element_text(size = 40)) + 
        scale_fill_brewer(palette="Spectral")
      print(bd)
    garbage <- dev.off()
    cat(" done\n")
  }

  
  # Number of msms.identification rate
  if(plotMSMS){
    cat("--- Plot Number of msms.identification rate")
    pdf('QC_Plots_summary_MSMS.pdf', width=nsamples*3, height=20, onefile = TRUE)
      ca <- ggplot(summarykeys, aes(x = bioreplicate, y = ms.ms.identified...., fill = condition)) +
        geom_bar(stat="identity", alpha=0.7) +
        geom_text(aes(label=round(ms.ms.identified...., digits=2)), vjust=1 , size = 15) +
        xlab("Experiment") + ylab("Rate") +
        ggtitle("MS2 Identification rate") +
        theme(legend.text = element_text(size=20)) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, size=20)) + 
        theme(axis.text.y = element_text(size=20)) + 
        theme(axis.title.x = element_text(size=30)) + 
        theme(axis.title.y = element_text(size=30)) + 
        theme(plot.title = element_text(size = 40)) + 
        scale_fill_brewer(palette="Spectral")
      print(ca)
      
      cb <- ggplot(summary2, aes(x = condition, y = pct.MS2Id.mean, fill = factor(condition))) +
        geom_bar(stat="identity", position = position_dodge(width = 1), alpha=0.7) +
        geom_errorbar(aes(ymin=pct.MS2Id.mean-pct.MS2Id.sem, ymax=pct.MS2Id.mean+pct.MS2Id.sem), width=.2, position=position_dodge(.9)) +
        geom_text(aes(label=round(pct.MS2Id.mean, digits=2)), hjust=0.5, vjust=-0.5, size = 10, position = position_dodge(width = 1)) +
        xlab("Condition") + ylab("Rate") +
        ggtitle("Mean MS2 Identification rate across bioreplicates and fractions") +
        theme(legend.text = element_text(size=20)) +
        theme(axis.text.x = element_text(angle = 0, hjust = 1, size=20)) + 
        theme(axis.text.y = element_text(size=20)) + 
        theme(axis.title.x = element_text(size=30)) + 
        theme(axis.title.y = element_text(size=30)) + 
        theme(plot.title = element_text(size = 40)) + 
        scale_fill_brewer(palette="Spectral")
      print(cb)
      
      if(isFractions){
        cc <- ggplot(summary2fx, aes(x = bioreplicate, y = ms.ms.identified...., fill = factor(fraction))) +
          geom_bar(stat="identity", alpha=0.7) +
          geom_text(aes(label=round(ms.ms.identified...., digits=1)), hjust=0.5, vjust=1.5, size = 7, position = position_stack()) +
          xlab("Experiment") + ylab("Counts") +
          ggtitle("MS2 Identification rate per Fraction") +
          theme(legend.text = element_text(size=20)) +
          theme(axis.text.x = element_text(angle = 90, hjust = 0, size=20)) + 
          theme(axis.text.y = element_text(size=20)) + 
          theme(axis.title.x = element_text(size=30)) + 
          theme(axis.title.y = element_text(size=30)) + 
          theme(plot.title = element_text(size = 40)) 
        print(cc)
      }
    garbage <- dev.off()
    cat(" done\n")
  }

  
  # Number of isotope patterns
  if(plotISOTOPE){
    cat("--- Plot Number of isotope patterns")
    pdf('QC_Plots_summary_ISOTOPE.pdf', width=nsamples*3, height=20, onefile = TRUE)
      da <- ggplot(summarykeys, aes(x = bioreplicate, y = isotope.patterns, fill = condition)) +
        geom_bar(stat="identity", alpha=0.7) +
        geom_text(aes(label=round(isotope.patterns, digits=0)), vjust=1 , size = 15) +
        xlab("Experiment") + ylab("Counts") +
        ggtitle("Number of detected Isotope Patterns") +
        theme(legend.text = element_text(size=20)) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, size=20)) + 
        theme(axis.text.y = element_text(size=20)) + 
        theme(axis.title.x = element_text(size=30)) + 
        theme(axis.title.y = element_text(size=30)) + 
        theme(plot.title = element_text(size = 40)) + 
        scale_fill_brewer(palette="Spectral")
      print(da)
      
      db <- ggplot(summary2, aes(x = condition, y = num.IsotopePatterns.mean, fill = factor(condition))) +
        geom_bar(stat="identity", position = position_dodge(width = 1), alpha=0.7) +
        geom_errorbar(aes(ymin=num.IsotopePatterns.mean-num.IsotopePatterns.sem, ymax=num.IsotopePatterns.mean+num.IsotopePatterns.sem), width=.2, position=position_dodge(.9)) +
        geom_text(aes(label=round(num.IsotopePatterns.mean, digits=0)), hjust=0.5, vjust=-0.5, size = 10, position = position_dodge(width = 1)) +
        xlab("Condition") + ylab("Counts") +
        ggtitle("Mean number of detected Isotope Patterns") +
        theme(legend.text = element_text(size=20)) +
        theme(axis.text.x = element_text(angle = 0, hjust = 1, size=20)) + 
        theme(axis.text.y = element_text(size=20)) + 
        theme(axis.title.x = element_text(size=30)) + 
        theme(axis.title.y = element_text(size=30)) + 
        theme(plot.title = element_text(size = 40)) + 
        scale_fill_brewer(palette="Spectral")
      print(db)
      
      if(isFractions){
        dc <- ggplot(summary2fx, aes(x = bioreplicate, y = isotope.patterns, fill = factor(fraction))) +
          geom_bar(stat="identity", alpha=0.7) +
          geom_text(aes(label=round(isotope.patterns, digits=0)), hjust=0.5, vjust=1.5, size = 7, position = position_stack()) +
          xlab("Experiment") + ylab("Counts") +
          ggtitle("Number of detected Isotope Patterns per Fraction") +
          theme(legend.text = element_text(size=20)) +
          theme(axis.text.x = element_text(angle = 90, hjust = 0, size=20)) + 
          theme(axis.text.y = element_text(size=20)) + 
          theme(axis.title.x = element_text(size=30)) + 
          theme(axis.title.y = element_text(size=30)) + 
          theme(plot.title = element_text(size = 40)) 
        print(dc)
      }
      
      # Number of sequenced isotope patterns with charge = 2 or more
      dd <- ggplot(summarykeys, aes(x = bioreplicate, y = isotope.patterns.sequenced..z.1., fill = condition)) +
        geom_bar(stat="identity", alpha=0.7) +
        geom_text(aes(label=round(isotope.patterns.sequenced..z.1., digits=0)), vjust=1 , size = 15) +
        xlab("Experiment") + ylab("Counts") +
        ggtitle("Number of sequenced Isotope Patterns with charge state greater than 1") +
        theme(legend.text = element_text(size=20)) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, size=20)) + 
        theme(axis.text.y = element_text(size=20)) + 
        theme(axis.title.x = element_text(size=30)) + 
        theme(axis.title.y = element_text(size=30)) + 
        theme(plot.title = element_text(size = 40)) + 
        scale_fill_brewer(palette="Spectral")
      print(dd)
      
      de <- ggplot(summary2, aes(x = condition, y = num.IsotopePatternsSeq.mean, fill = factor(condition))) +
        geom_bar(stat="identity", position = position_dodge(width = 1), alpha=0.7) +
        geom_errorbar(aes(ymin=num.IsotopePatternsSeq.mean-num.IsotopePatternsSeq.sem, ymax=num.IsotopePatternsSeq.mean+num.IsotopePatternsSeq.sem), width=.2, position=position_dodge(.9)) +
        geom_text_repel(aes(label=round(num.IsotopePatternsSeq.mean, digits=0)), hjust=0.5, vjust=-0.5, size = 10, position = position_dodge(width = 1)) +
        xlab("Condition") + ylab("Counts") +
        ggtitle("Mean number of sequenced Isotope Patterns with charge state greater than 1") +
        theme(legend.text = element_text(size=20)) +
        theme(axis.text.x = element_text(angle = 0, hjust = 1, size=20)) + 
        theme(axis.text.y = element_text(size=20)) + 
        theme(axis.title.x = element_text(size=30)) + 
        theme(axis.title.y = element_text(size=30)) + 
        theme(plot.title = element_text(size = 40)) + 
        scale_fill_brewer(palette="Spectral")
      print(de)
      
      if(isFractions){
        df <- ggplot(summary2fx, aes(x = bioreplicate, y = isotope.patterns.sequenced..z.1., fill = factor(fraction))) +
          geom_bar(stat="identity", alpha=0.7) +
          geom_text(aes(label=round(isotope.patterns.sequenced..z.1., digits=0)), hjust=0.5, vjust=1.5, size = 7, position = position_stack()) +
          xlab("Experiment") + ylab("Counts") +
          ggtitle("Number of sequenced Isotope Patterns  with charge state greater than 1 per Fraction") +
          theme(legend.text = element_text(size=20)) +
          theme(axis.text.x = element_text(angle = 90, hjust = 0, size=20)) + 
          theme(axis.text.y = element_text(size=20)) + 
          theme(axis.title.x = element_text(size=30)) + 
          theme(axis.title.y = element_text(size=30)) + 
          theme(plot.title = element_text(size = 40)) 
        print(df)
      }
    garbage <- dev.off()
    cat(" done\n")
  }

} # END OF SUMMARY



