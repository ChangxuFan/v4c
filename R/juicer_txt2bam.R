juicer.txt2bam <- function(juicer.txt, mapq, out.bam = NULL, threads, genome, ...) {
  if (is.null(out.bam))
    out.bam <- sub(".txt$", ".bam", juicer.txt)
  out.sam <- sub(".bam$", ".sam", out.bam)
  awk.cmd <- paste0("awk -F \" \" 'BEGIN {OFS = \"\t\"} $9 > ", mapq,
                    "  && $12 > ", mapq, " {print $15, $1, $2, $3, $9, $10, \"*\", \"0\", \"0\", $11, gensub(/./,\"F\",\"g\", $11); print $16, $5, $6, $7, $12, $13, \"*\", \"0\", \"0\", $14, gensub(/./,\"F\",\"g\", $14) } ' ",
                    juicer.txt, " > ", out.sam)
  paste0(awk.cmd)
  system(awk.cmd)
  juicer.sam2bam(sam = out.sam, out.bam = out.bam, threads = threads, genome = genome, ...) %>%
    return()

}

juicer.txt2bam.2 <- function(juicer.txt, mapq, out.bam = NULL, threads, genome, ...) {
  options(scipen = 999)
  if (is.null(out.bam))
    out.bam <- sub(".txt$", ".bam", juicer.txt)

  out.sam <- sub(".bam$", ".sam", out.bam)
  system(paste0("rm -rf ", out.sam))

  f <- function(x, pos) {
    x <- x[x$X9 > mapq & x$X12 > mapq, ]
    if (nrow(x) < 1)
      return(NULL)
    df <- data.frame(rname = c(x$X15, x$X16), flag = c(x$X1, x$X5), chr = c(x$X2, x$X6),
                     pos = c(x$X3, x$X7), mapq = c(x$X9, x$X12), cigar = c(x$X10, x$X13),
                     chr.mate = "*", pos.mate = 0, tlen = 0, seq = c(x$X11, x$X14),
                     qual = c(x$X11, x$X14) %>% gsub(".", "F", .), stringsAsFactors = F)

    df$length <- -1 * df$cigar %>% GenomicAlignments::cigarWidthAlongReferenceSpace(cigar = .) + 1

    df[df$flag != 16,"length"] <- 0

    df$pos <- df$pos %>% as.numeric() + df$length
    #print(df)
    df$length <- NULL
    # print(df)
    write.table(df, out.sam, sep = "\t", quote = F, col.names = F, row.names = F, append = T)
    return()
  }

  reads <- readr::read_delim_chunked(file = juicer.txt, delim = " ", col_names = F, chunk_size = 100000,
                                 callback = f,
                                 col_types = rep("c", 17) %>% paste0(collapse = ""))

  juicer.sam2bam(sam = out.sam, out.bam = out.bam, threads = threads, genome = genome, ...) %>%
    return()

}

juicer.txt2sam <- function(merge.sam=NULL, juicer.sam=NULL, n = 1000000,
                           out.sam = NULL) {
  options(scipen = 999)
  options(stringsAsFactors = F)
  merge.df <- readr::read_tsv(merge.sam,
                              # col_types = rep("c", 11) %>% paste0(collapse = ""),
                              col_names = F)
  rnames <- merge.df$X1

  if (is.null(out.sam)) {
    out.sam <- sub(".sam$",".nodups.sam",juicer.sam)
  }

  system(paste0("rm -rf ", out.sam))

  con <- file(juicer.sam, "r", blocking = T)

  while (1) {
    lines <- readLines(con, n = n)
    r <- gsub("\\t.+", "", lines)
    lines.filtered <- lines[r %in% rnames]
    write(x = lines.filtered, file = out.sam, sep = "\n", append = T)
    if (length(lines) < 1)
      break
  }

  try(close(con = con))

  return(out.sam)


}
