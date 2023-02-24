juicer.txt2bam <- function(juicer.txt, mapq, out.bam = NULL, sam2bam = T, threads, genome, ...) {
  if (is.null(out.bam))
    out.bam <- sub(".txt$", ".bam", juicer.txt)
  out.sam <- sub(".bam$", ".sam", out.bam)
  awk.cmd <- paste0("awk -F \" \" 'BEGIN {OFS = \"\t\"} $9 > ", mapq,
                    "  && $12 > ", mapq, " {print $15, $1, $2, $3, $9, $10, \"*\", \"0\", \"0\", $11, gensub(/./,\"F\",\"g\", $11); print $16, $5, $6, $7, $12, $13, \"*\", \"0\", \"0\", $14, gensub(/./,\"F\",\"g\", $14) } ' ",
                    juicer.txt, " > ", out.sam)
  paste0(awk.cmd)
  system(awk.cmd)
  if (sam2bam == T) {
    trash <- juicer.sam2bam(sam = out.sam, out.bam = out.bam, threads = threads, genome = genome, ...) 
    return(out.bam)
  } else {
    return(out.sam)
  }

}

juicer.txt.getName <- function(juicer.txt, mapq, inverse = F, region = NULL, out.file = NULL, run = T) {
  # region should be written in chr:start-end format.
  cmd.head <- "awk -F \" \" 'BEGIN {OFS = \"\t\"} "
  if (!is.null(region)) {
    if (length(region) != 1)
      stop("length of -region parameter must be 1")
    region <- utilsFanc::loci.2.df(loci.vec = region)

    cmd.filter.pos1 <- paste0("(", "$2 == \"", region$chr, 
                                "\" && ($3 >= ", region$start,
                                " && $3 <= ", region$end, "))")
    cmd.filter.pos2 <- paste0("(", "$6 == \"", region$chr, 
                              "\" && ($7 >= ", region$start,
                              " && $7 <= ", region$end, "))")
    
    cmd.filter.region <- paste0("(", cmd.filter.pos1, " || ", cmd.filter.pos2, ") && ")
    region.string <- paste0("_", region$chr, ":", region$start, "-", region$end)
  } else {
    cmd.filter.region <- NULL
    region.string <- ""
  }
  
  if (inverse == F) {
    cmd.filter.mapq <- paste0("($9 > ", mapq,"  && $12 > ", mapq, ")")
    inverse.string <- ""
  } else {
    cmd.filter.mapq <- paste0("($9 <= ", mapq, "  || $12 <= ", mapq, ")")
    inverse.string <- "_inverse"
  }
  
  if (is.null(out.file))
    out.file <- sub(".txt$", paste0(region.string, "_mapq",mapq,"_rname",inverse.string ,".txt"), juicer.txt)
  
  cmd.tail <- paste0(" {print $15}' ",
                     juicer.txt, " > ", out.file)
  
  cmd <- paste0(cmd.head, cmd.filter.region, cmd.filter.mapq, cmd.tail)
  
  # awk.cmd <- paste0("awk -F \" \" 'BEGIN {OFS = \"\t\"} $9 > ", mapq,
  #                   "  && $12 > ", mapq, " {print $15}' ",
  #                   juicer.txt, " > ", out.file)
  
  utilsFanc::cmd.exec.fanc(cmd, run = run, intern = F)
  return(out.file)
}


juicer.txt.filter.ez <- function(rname.file = NULL, inverse = F,
                                 region = NULL,
                                 juicer.sam, juicer.txt, mapq, out.bam = NULL,
                                 reverse.grep = F,
                                 n.line.chunk = 2000000, overwrite.splits = T,
                                 genome, header.file = NULL,
                                 threads, run = T, samtools = liteRnaSeqFanc::SAMTOOLS) {
  if (is.null(rname.file)) {
    rname.file <- juicer.txt.getName(juicer.txt = juicer.txt, mapq = mapq, run = run,
                                     region = region, inverse = inverse)
  }
  if (is.null(out.bam)) {
    if (inverse == T) {
      inverse.string <- "_inverse"
    } else {
      inverse.string <- ""
    }
    if (!is.null(region)) {
      region.string <- paste0("_", region)
      if (length(region.string) > 1) {
        stop("region must be of length 1")
      }
    } else {
      region.string <- ""
    }
    out.bam <- juicer.sam %>%  
      sub(".sam$", paste0(region.string, "_nodups_mapq", mapq, inverse.string, ".bam"),. )
  }
  
  out.sam <- sub(".bam$", ".sam", out.bam)
  
  # if (!is.null(n.line.chunk)) {
  #   if (overwrite.splits == T) {
  #     system(paste0("rm -rf ", rname.file, ".split*"))
  #     wcl <- system(paste0("wc -l ", rname.file), intern = T) %>% sub(" +.+$", "", .) %>% as.numeric()
  #     cmd <- paste0("split --numeric-suffixes -l ", n.line.chunk, " ", rname.file, " " ,rname.file, ".split")
  #     print(cmd)
  #     system(cmd)
  #   } 
  #   rname.files <- Sys.glob(paste0(rname.file, ".split*"))
  # } else {
  #   rname.files <- rname.file
  # }
  # 
  if (reverse.grep == T) {
    sam.rname.file <- paste0(juicer.sam, ".rnames.txt")
    intersect.file <- paste0(juicer.sam, ".rnames", "_nodups_mapq", mapq,".txt")
    cmd <- paste0(samtools, " view ", juicer.sam, " | awk -F '\\t' '{print $1}' > ", sam.rname.file )
    print(cmd); system(cmd)
    cmd <- paste0("/bin/bash -c ", "'comm -12 <(sort ", rname.file, ")", " ",
                  "<(sort ", sam.rname.file, ")", " > ", intersect.file, "'")
    print(cmd); system(cmd)
    rname.file <- intersect.file
  }
  
  cmd <- paste0("export LC_ALL=C && ",
                "parallel -j ", threads, " --pipepart -a ", juicer.sam,
                " --block 30M fgrep -F -f ", rname.file, " ", ">", " ",out.sam)
  utilsFanc::cmd.exec.fanc(cmd = cmd, intern = F, run = run)
  
  # for (i in 1:length(rname.files)) {
  #   rname.file <- rname.files[i]
  #   if (i == 1) 
  #     pipe <- ">"
  #   else
  #     pipe <- ">>"
  # 
  # }
  
  
  out.bam <- juicer.sam2bam(juicer.dir = NULL, sam = out.sam, add.header = T, header.file = header.file,
                            genome = genome, out.bam = out.bam, threads = threads, 
                            samtools = samtools)
  
  if (!file.exists(out.bam))
    stop(paste0(out.bam, " was not succesfully generated"))
  return(out.bam)
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

juicer.txt.filter <- function(merge.sam=NULL, juicer.sam=NULL, n = 1000000,
                           out.sam = NULL, 
                           sam2bam = T, genome, threads) {
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

  try({
    while (1) {
      lines <- readLines(con, n = n)
      r <- gsub("\\t.+", "", lines)
      lines.filtered <- lines[r %in% rnames]
      write(x = lines.filtered, file = out.sam, sep = "\n", append = T)
      if (length(lines) < 1)
        break
    }
  })

  close(con = con)
  
  if (sam2bam == T) {
    out.bam <- juicer.sam2bam(sam = out.sam, threads = threads, genome = genome, ...) 
    return(out.bam)
  }
  return(out.sam)


}

juicer.txt2sam <- juicer.txt.filter