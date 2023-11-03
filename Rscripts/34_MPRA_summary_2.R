
suppressMessages(library("plyr", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("data.table", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("crayon", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("withr", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("ggplot2", lib.loc = "/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("farver", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("labeling", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("optparse", lib.loc = "/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("dplyr", lib.loc = "/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("withr", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("backports", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("broom", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("rstudioapi", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("cli", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("tzdb", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("svglite", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("ggeasy", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("sandwich", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("digest", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("tidyverse", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("BiocGenerics", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("S4Vectors", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("IRanges", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("GenomeInfoDb", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("GenomicRanges", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("Biobase", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("AnnotationDbi", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("GO.db", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("org.Hs.eg.db", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("TxDb.Hsapiens.UCSC.hg19.knownGene", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("rtracklayer", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))

opt = NULL

options(warn = 1)

data_wrangling_MPRA = function(option_list)
{

  opt_in = option_list
  opt <<- option_list
  
  cat("All options:\n")
  printList(opt)
  
  
  #### READ and transform type ----
  
  type = opt$type
  
  cat("TYPE_\n")
  cat(sprintf(as.character(type)))
  cat("\n")

  
  #### READ and transform out ----
  
  out = opt$out
  
  cat("out_\n")
  cat(sprintf(as.character(out)))
  cat("\n")
  
 
  #### READ MPRA_results ----
  
  MPRA_results<-as.data.frame(fread(file=opt$MPRA_results, sep="\t", header=T), stringsAsFactors=F)
  
  
  cat("MPRA_results_0\n")
  cat(str(MPRA_results))
  cat("\n")
  cat(str(unique(MPRA_results$VAR)))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(MPRA_results$ASSAY_CLASS))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(MPRA_results$ASSAY_CLASS)))))
  cat("\n")
  
  MPRA_results$carried_to_VAR<-paste('chr',MPRA_results$carried_variants,sep='')
  
  
  MPRA_results_subset<-unique(MPRA_results[which(MPRA_results$carried_to_VAR == MPRA_results$VAR &
                                            MPRA_results$ASSAY_CLASS == 'Screened variant'),])
  
  cat("MPRA_results_subset_0\n")
  cat(str(MPRA_results_subset))
  cat("\n")
  cat(str(unique(MPRA_results_subset$VAR)))
  cat("\n")
  
  ##### collapse -----
  
  MPRA_results_subset.dt<-data.table(MPRA_results_subset, key="VAR")
  
  MPRA_results_subset_collapsed<-as.data.frame(MPRA_results_subset.dt[,.(string_Per_tile_experimental_class=paste(unique(Per_tile_experimental_class), collapse = '|')), by=key(MPRA_results_subset.dt)], stringsAsFactors=F)
  
  cat("MPRA_results_subset_collapsed_0\n")
  cat(str(MPRA_results_subset_collapsed))
  cat("\n")
  cat(str(unique(MPRA_results_subset_collapsed$VAR)))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(MPRA_results_subset_collapsed$string_Per_tile_experimental_class))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(MPRA_results_subset_collapsed$string_Per_tile_experimental_class)))))
  cat("\n")
  
  MPRA_results_subset_collapsed$consolidated_EA<-NA
  
  MPRA_results_subset_collapsed$consolidated_EA[grep(paste("EA",'EA&ASE', sep='|'), MPRA_results_subset_collapsed$string_Per_tile_experimental_class)]<-'EA'
 
  cat(sprintf(as.character(names(summary(as.factor(MPRA_results_subset_collapsed$consolidated_EA))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(MPRA_results_subset_collapsed$consolidated_EA)))))
  cat("\n")
  
  MPRA_results_subset_collapsed$consolidated_ASE<-NA
  
  MPRA_results_subset_collapsed$consolidated_ASE[grep(paste("ASE",'EA&ASE', sep='|'), MPRA_results_subset_collapsed$string_Per_tile_experimental_class)]<-'ASE'
  
  cat(sprintf(as.character(names(summary(as.factor(MPRA_results_subset_collapsed$consolidated_ASE))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(MPRA_results_subset_collapsed$consolidated_ASE)))))
  cat("\n")
  
  MPRA_results_subset_collapsed$consolidated_EA_and_ASE<-NA
  
  MPRA_results_subset_collapsed$consolidated_EA_and_ASE[grep('EA&ASE', MPRA_results_subset_collapsed$string_Per_tile_experimental_class)]<-'EA&ASE'
  
  cat(sprintf(as.character(names(summary(as.factor(MPRA_results_subset_collapsed$consolidated_EA_and_ASE))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(MPRA_results_subset_collapsed$consolidated_EA_and_ASE)))))
  cat("\n")
  
  
  MPRA_results_subset_collapsed_subset<-unique(MPRA_results_subset_collapsed[,c(which(colnames(MPRA_results_subset_collapsed) == 'VAR'),which(colnames(MPRA_results_subset_collapsed) == 'consolidated_EA'),
                                                                         which(colnames(MPRA_results_subset_collapsed) == 'consolidated_ASE'),which(colnames(MPRA_results_subset_collapsed) == 'consolidated_EA_and_ASE'))])
  
  cat("MPRA_results_subset_collapsed_subset_0\n")
  cat(str(MPRA_results_subset_collapsed_subset))
  cat("\n")
  
  

  #### SAVE ----

  setwd(out)


  write.table(MPRA_results_subset_collapsed_subset, file=paste("MPRA_results_subset_collapsed",".tsv",sep=''),sep="\t", quote=F,row.names = F)
  saveRDS(MPRA_results_subset_collapsed_subset, file=paste('MPRA_results_subset_collapsed',".rds",sep=''))

}




printList = function(l, prefix = "    ") {
  list.df = data.frame(val_name = names(l), value = as.character(l))
  list_strs = apply(list.df, MARGIN = 1, FUN = function(x) { paste(x, collapse = " = ")})
  cat(paste(paste(paste0(prefix, list_strs), collapse = "\n"), "\n"))
}


#### main script ----

main = function() {
  cmd_line = commandArgs()
  cat("Command line:\n")
  cat(paste(gsub("--file=", "", cmd_line[4], fixed=T),
            paste(cmd_line[6:length(cmd_line)], collapse = " "),
            "\n\n"))
  option_list <- list(
    make_option(c("--MPRA_results"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--type"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--out"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required.")
  )
  parser = OptionParser(usage = "140__Rscript_v106.R
                        --subset type
                        --TranscriptEXP FILE.txt
                        --cadd FILE.txt
                        --ncboost FILE.txt
                        --type type
                        --out filename",
                        option_list = option_list)
  opt <<- parse_args(parser)
  
  data_wrangling_MPRA(opt)
  
  
}


###########################################################################

system.time( main() )