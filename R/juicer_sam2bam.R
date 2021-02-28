juicer.sam2bam <- function(juicer.dir=NULL, sam = NULL, genome = NULL, space.to.tab = F, 
                           header.file = NULL, out.bam = NULL, threads,
                           samtools = "/bar/cfan/anaconda2/envs/jupyter/bin/samtools") {
  if (is.null(header.file)) {
    header.file <- paste0("/bar/cfan/genomes/", genome, "/bam_header")
  }

  if (is.null(sam)) {
    sam <- paste0(juicer.dir, "/splits/*.fastq.gz.sam") %>% Sys.glob() %>% normalizePath(mustWork = T)
    if (length(sam) > 1)
      stop("more than 1 sam file captured")
  }

  if (is.null(out.bam)) {
    out.bam <- sub(".sam$", ".bam", sam)
  }
  bam.tempt <- paste0(out.bam, ".tempt")
  sam.header <- paste0(sam, ".header")
  if (space.to.tab == T) {
      cmd <- paste0("cat ", header.file, " > ", sam.header)
      print(cmd)
      system(cmd)
      cmd <- paste0("sed 's/ /\\t/g' ", sam, " >> ", sam.header)
      print(cmd)
      system(cmd)

    } else {
      cmd <- paste0("cat ", header.file, " ", sam, " > ", sam.header)
        print(cmd)
        system(cmd)
    }



  cmd <- paste0(samtools, " view -h -b -@ ", threads, " ", sam.header, " > ", bam.tempt)
  print(cmd)
  system(cmd)

  cmd <- paste0(samtools, " sort -O bam -m 4G -@ ", threads, " ", bam.tempt , " > ", out.bam)
  print(cmd)
  system(cmd)

  cmd <- paste0(samtools, " index ", out.bam)
  print(cmd)
  system(cmd)
  return(out.bam)

}
