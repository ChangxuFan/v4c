feed.sureDesign <- function(cutsite.beds, region.gr, flank.size, out.bed) {
  
  cutsite.gr <- lapply(cutsite.beds, function(x) rtracklayer::import.bed(x, which = region.gr))
  cutsite.gr <- Reduce(c, cutsite.gr) %>% GenomicRanges::reduce()
  # take up flank:
  upflank.gr <- flank(cutsite.gr, width = flank.size, start = T, both = F, ignore.strand = T)
  dnflank.gr <- flank(cutsite.gr, width = flank.size, start = F, both = F, ignore.strand = T)
  flank <- c(upflank.gr, dnflank.gr) %>% GenomicRanges::reduce()
  utilsFanc::write.zip.fanc(df = flank, bed.shift = T, out.file = out.bed, zip = T)
  flank.txt <- utilsFanc::gr.get.loci(gr = flank)
  write(flank.txt, sub(".bed$", ".txt", out.bed), sep = "\n")
  utilsFanc::write.zip.fanc(df = flank, bed.shift = T, out.file = sub(".bed$", ".txt", out.bed), zip = F)
  utilsFanc::bash2ftp(paste0(out.bed, ".gz"))
  return(out.bed)
}