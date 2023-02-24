pct.capture <- function(bams, df, threads, outfile = NULL) {
  stat <- utilsFanc::safelapply(seq_along(bams), function(i) {
  	  bam <- bams[i]
  	  sample.name <- names(bams)[i]
  	  total <- Rsamtools::countBam(file = bam)$records
  	  # print(total)
      gr <- GenomicRanges::makeGRangesFromDataFrame(df)
      stat <- Rsamtools::countBam(file = bam, param = Rsamtools::ScanBamParam(which = gr))
      stat$total.reads <- total
      stat$frac.on.target <- stat$records/stat$total
      stat <- utilsFanc::change.name.fanc(stat, "records", "reads.on.target")
      stat$nucleotides <- NULL
      if (is.null(sample.name))
      	sample.name <- bam
      stat <- utilsFanc::add.column.fanc(df1 = stat, df2 = data.frame(sample.name = rep(sample.name, nrow(stat))), pos = 1)
      return(stat)
  	}, threads = threads) %>% Reduce(rbind, .)


  #stat <- data.frame(regional = regional, total = total, pct = regional/total)
  if (!is.null(outfile)) {
  	write.table(stat, outfile, 	sep = "\t", row.names = F, col.names = T, quote = F)
  }
  return(stat)
}

