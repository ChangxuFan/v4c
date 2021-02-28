#
V4C.juicer <- function (pos.list, length, resolution, genome.build, out.root.name=NULL, plot, observed_or_oe, normalization,
                        hic.file, help.only, bdg.dir=NULL, shift=0, add.vline=NULL, liftover.params=NULL,
                        use.1.6 = F, strip.chr = F,...) {

  if (!is.null(bdg.dir))
    plot <- F

  
  genome.bins.directory <- paste0("~/genomes/", genome.build, "/bins/")

  prox.file <- paste0(genome.bins.directory, "/", genome.build, "_", resolution, "_bins.prox")
  if (!file.exists(prox.file)) {
    system(paste0("~/scripts/bin_gen.sh ",
                  paste0("~/genomes/", genome.build, "/",genome.build, ".chrom.sizes "),
                  resolution, " ",
                  paste0("~/genomes/", genome.build, "/", genome.build, "_TSS.bed")))
  }

  npos.list <- as.integer(sub(".*:(.*)", "\\1", pos.list)) %>% resolution.convert(resolution)
  chr.list <- sub("(.*):.*", "\\1", pos.list)
  pos.list <- paste0(chr.list,":", npos.list)
  # these 3 lines seems trivial but they exist to convert resolution.
  if (help.only==T) {
    print("observed_or_oe and normalization should be in the format required by juicer_tools.jar dump")
    print("this is not an error even but R will say it's an error since I used stop()")
    stop()
  }
  if (!observed_or_oe %in% c("observed", "oe")) stop("observed_or_oe must be set as 'observed' or 'oe'")
  if (!normalization %in% c("NONE","VC","VC_SQRT", "KR")) stop('normalization must be one of "NONE","VC","VC_SQRT", "KR"')
  prox.df <- read.table(prox.file, as.is = T)
  prox.df$prox <- paste0(prox.df$V1, ":", prox.df$V2)
  t.dir <- tempdir()

  trash <- lapply(pos.list, function(pos) {
    npos <- as.integer(sub(".*:(.*)", "\\1", pos))
    chr <- sub("(.*):.*", "\\1", pos)
    distalStart=npos-length
    distalEnd=npos+length
    print(paste0(chr, ":", distalStart, ":", distalEnd))
    gene <- filter(prox.df, prox==pos) %>% pull(V4)
    if (length(gene) < 1)
      gene <- "not_TSS"
    juicer.tools <- "~/software/juicer/SLURM/scripts/juicer_tools.jar"
    if (use.1.6 == T)
      juicer.tools <- "~/software/juicer-1.6/SLURM/scripts/juicer_tools.jar"
    chr.juicer <- chr
    if (strip.chr == T)
      chr.juicer <- sub("chr", "", chr)
    pos.juicer <- pos
    if (strip.chr == T)
      pos.juicer <- sub("chr", "", pos)
    cmd <- paste0("/opt/apps/java/jdk1.8.0_92/bin/java -jar ",juicer.tools," dump ", observed_or_oe, " ", normalization,
                  " ", hic.file, " ", paste0(pos.juicer, ":", npos+1), " ", paste0(chr.juicer, ":", distalStart, ":", distalEnd),
                  " ", "BP", " ", resolution, " ", paste0(t.dir, "/", pos, ".txt"))

    print(cmd)
    system(cmd)

    print(paste0("~/scripts/hic/v4/hic2v4.sh ", chr, " ", gene, " ",paste0(t.dir, "/", pos, ".txt") , " ", npos))
    system(paste0("~/scripts/hic/v4/hic2v4.sh ", chr, " ", gene, " ",paste0(t.dir, "/", pos, ".txt") , " ", npos))
    return(NULL)
  })
  v4.list <- lapply(pos.list, function (x) {
    return(V4C(pos = x, length=length, resolution=resolution, plot=plot, normalization=paste0( normalization, ".",observed_or_oe),
               indir=t.dir, ...))
  } )
  if (!is.null(bdg.dir)) {
    v4c2bdg(v4c.list = v4.list, bdg.dir = bdg.dir, hic.file = hic.file,
            length = length, resolution = resolution, shift = shift,
            add.vline = add.vline, liftover.params = liftover.params, 
            root.name = out.root.name, genome = genome.build)
  }

  return(v4.list)
}



V4C.cool <- function (pos.list, length=500000, resolution, prox.file, plot=T,
                      cool.file, ...) {
  stop("this function is not correct at this moment. Something wrong with the normalization")
  # note: this function works only on mcool files, not cool files! the difference between the 2 is that mcool is
  # multiple cool files from different resolutions joined together.
  prox.df <- read.table(prox.file, as.is = T)
  prox.df$prox <- paste0(prox.df$V1, ":", prox.df$V2)
  file.to.py <- lapply(pos.list, function (pos) {
    npos <- as.integer(sub(".*:(.*)", "\\1", pos))
    chr <- sub("(.*):.*", "\\1", pos)
    if (npos-length < 0) {
      length <- floor(npos/resolution)*resolution
    }
    gene <- filter(prox.df, prox==pos) %>% pull(V4)
    df <- data.frame( chr=chr,
                      proxStart=npos,
                      proxEnd=npos+resolution,
                      distalStart=npos-length,
                      distalEnd=npos+length,
                      gene=gene)
    return(df)
  })
  file.to.py <- Reduce(rbind, file.to.py)
  temp <- tempfile(fileext = ".v4c")
  t.dir <- tempdir()
  print(t.dir)
  write.table(x = file.to.py, file = temp, sep = "\t", row.names = F, col.names = T, quote = F)
  #system(paste0("cp ", temp, " /bar/cfan/SIPG/micro-c/Dekker/test.v4c"))
  system(paste0("python ~/python_scripts/v4c_cool.py ", temp, " ", cool.file, " ", t.dir, " ", resolution))
  v4.list <- lapply(pos.list, function (x) {
    return(V4C(pos = x, length=length, resolution=resolution, plot=plot,return.plot = return.plot,
               indir=t.dir, ...))
  } )


}

# p <- V4C.cool(gene.name.to.location(B.cell.gene, "~/genomes/hg38/bins/hg38_1000_bins.prox"),
#               resolution=1000, prox.file="~/genomes/hg38/bins/hg38_1000_bins.prox", plot = T, cool.file = "~/SIPG/micro-c/Dekker/4DNFI9GMP2J8.mcool",
#               standard.x=F, rm.zero = T, find.enhancer = T, length = 30000,
#               find.loop = TRUE, loop.file = "~/SIPG/hic/V4/loops/GM_looplist.simple", enhancer.extend = 0,
#               loop.extend = 0, smooth.line = F, loop.anchor.as.lines = TRUE)

V4C <- function(pos, length, resolution, plot=T, ylim=NULL, normalization,
                indir, save.plot.dir=NULL, bedfile.dir=NULL, standard.x=F, rm.zero=T, find.enhancer=FALSE,
                find.loop=FALSE, loop.file, loop.anchor.as.lines=TRUE, n.ticks=10, return.plot=T,
                HMM.enhancer, enhancer.extend, loop.extend,
                smooth.line=FALSE, add.vline=NULL, title = "", pt.size = 0.5) {
  # note: normalization is not performed at this step. it's only used to plot ylab.
  m <- read.table(paste0(indir, "/", pos,".txt"), as.is = TRUE)
  colnames(m) <- c("prox",	"distal", "gene",	"value")
  npos <- as.integer(sub(".*:(.*)", "\\1", pos))
  if (npos-length < 0) {
    length <- floor(npos/resolution)*resolution
  }
  n.bin <- seq(npos-length, npos+length, resolution)
  full <- data.frame(prox=rep(pos, length(n.bin)),
                     distal=n.bin)
  full <- left_join(full, m)
  full[is.na(full)] <- 0
  # a position frame is used because in m, interactions with 0 values are omitted.
  #print(full)
  full$epi <- rep("N", nrow(full))
  if (find.enhancer==T) {
    df <- data.frame(chr=rep(sub("(.*):.*", "\\1", pos), nrow(full)),
                     left=full$distal)
    df$right <- df$left + resolution-1
    enhancer <- find.enhancer(df = df, HMM.enhancer = HMM.enhancer, extend = enhancer.extend)
    t.v.e <- rep("", nrow(full))
    t.v.e[full$distal %in% enhancer] <- ":enhancer"
    full$epi <- paste0(full$epi, t.v.e)
  }

  full$structure <- rep("N", nrow(full))
  if (find.loop==TRUE) {
    df <- data.frame(chr=rep(sub("(.*):.*", "\\1", pos), nrow(full)),
                     left=full$distal)
    df$right <- df$left + resolution-1
    l.a <- find.loop(df = df, loop.file = loop.file, extend = loop.extend)
    t.v.l <- rep("", nrow(full))
    t.v.l[full$distal %in% l.a] <- ":loop_anchor"
    full$structure <- paste0(full$structure, t.v.l)
  }


  if (standard.x==T) { npos <- 0; full$distal <- ((npos-length)/resolution):((npos+length)/resolution)}
  if (plot == T) {
    p.full <- full[,c("distal", "value", "epi", "structure")]
    colnames(p.full) <- c("x", "y","epi", "structure")
    if(rm.zero==T) p.full <- filter(p.full, y > 0)
    breaks <- break.generator(p.full$x, n.ticks)
    p <- ggplot(p.full, aes(x=x,y=y))

    if (!is.null(add.vline)) {
      p <- p+geom_vline(xintercept = add.vline %>% unlist(), color="green")
    }
    p <- p + geom_point(aes(color=epi, shape=structure), size = pt.size)+
      geom_line() +
      theme_classic() + xlab("distal position") + ylab(normalization) +
      ggtitle(paste0(pos, " ",m$gene[1], "\n", "bin width: ",resolution, "\n", title)) +
      geom_vline(xintercept = npos, color="black") +
      #geom_hline(yintercept = 1, color="grey") +
      scale_x_continuous(breaks=breaks) +
      theme(axis.text.x = element_text(angle = 45, hjust=1))
    if (sum(p.full$epi != "N")+sum(p.full$structure != "N") == 0)
      p <- p+theme(legend.position = "none")
    if (smooth.line) p <- p+geom_smooth(span=0.1, se=T)
    if (!is.null(ylim)) p <- p+ylim(0,ylim)
    if (loop.anchor.as.lines && find.loop) {
      x.i <- p.full$x[which(p.full$structure=="N:loop_anchor")]
      p <- p+geom_vline(xintercept = x.i, color="tan3", linetype="dashed")
    }

    if (grepl("oe", normalization))
      p <- p+geom_hline(yintercept = 1)
    if (!is.null(save.plot.dir)) {
      ggsave(save.plot.dir,p,device = "png", width = 10, height=3, units = "in",
             dpi = 100)
    }

    interval <- paste0(sub("(.*):.*", "\\1", pos),":",breaks)[c(1, length(breaks))]

    bed <- data.frame(chr=rep(sub("(.*):.*", "\\1", pos), length(breaks)), left=breaks, right=breaks+50, value=1)

    if (!is.null(bedfile.dir)) {
      write.table(bed, bedfile.dir, quote = FALSE, sep = "\t", row.names = F, col.names = F)
    }

    if (return.plot==T) return(p)
    else return(list(interval=interval, bed=bed))
  } else return(full)
}

break.generator <- function(x, n.ticks) {
  int <- (max(x)-min(x))/(2*n.ticks)
  return(seq(min(x), max(x), int))
}

find.loop <- function(df, loop.file, extend=0) {
  # loop.file should be a 3 column bed file containing the coordinates of loop anchors
  df[,2] <- df[,2] - extend
  df[,2][which(df[,2] < 1)] <- 1
  df[,3] <- df[,3] + extend
  #l.a means loop anchor
  l.a <- read.table(loop.file, as.is=TRUE) %>% arrange(V1, V2)
  df <- arrange(df, chr, left)
  o <- TEenrich::fanc.bedr.join.region(x = df, y = l.a, params = " -wao")
  colnames(o)[length(colnames(o))] <- "overlap"
  o <- filter(o, overlap > 0)
  return(o[,2])
}

find.enhancer <- function(df, HMM.enhancer, extend=1000) {
  df[,2] <- df[,2] - extend
  df[,2][which(df[,2] < 1)] <- 1
  df[,3] <- df[,3] + extend
  HMM.enhancer <- arrange(HMM.enhancer,V1, V2)
  df <- arrange(df, chr, left)
  o <- TEenrich::fanc.bedr.join.region(x = df, y = HMM.enhancer, params = " -wao -F 0.8")
  colnames(o)[length(colnames(o))] <- "overlap"
  o <- filter(o, overlap > 0)
  return(o[,2])
}



anchor2bed <- function(anchor, increment=5000) {
  df <- data.frame(chr=sub("(.*):(.*)", "\\1", anchor),
                   left=as.numeric(sub("(.*):(.*)", "\\2", anchor)),
                   right=as.numeric(sub("(.*):(.*)", "\\2", anchor)) + increment-1)
  df$id <- 1:nrow(df)
  return(df)
}

gene.name.to.location <- function (gene.name, prox.file, rm.nonCanonical=T, canonical.chr=NULL) {
  # takes in gene.name as a vector
  # returns a vector of loci
  prox <- read.table(prox.file, as.is = TRUE)
  prox <- separate_rows(prox, V4)
  gene.df <- data.frame(gene=gene.name)
  colnames(prox) <- c("chr", "left", "right", "gene")
  if (rm.nonCanonical==T)
  {
    if (is.null(canonical.chr)) {
      print("chr1-22, chrX, chrY")
      canonical.chr=paste0("chr", c(1:22, "X", "Y"))
    }
    prox <- filter(prox, chr %in% canonical.chr )
  }
  df <- left_join(gene.df, prox)
  locations <- paste0(df$chr, ":", df$left)
  return(locations)
}

resolution.convert <- function(x, out.res) {
  # x is a vector.
  x.converted <- floor(x/out.res)*out.res
  return(x.converted)
}

v4c2bdg <- function(v4c.list, bdg.dir, hic.file, length, root.name = NULL, genome = NULL,
                    resolution, add.vline=NULL, shift=0, liftover.params = NULL) {
  if (is.null(root.name)) {
    root.name <- basename(hic.file)
  }
  options(scipen = 20)
  system(paste0("mkdir -p ", bdg.dir))
  lapply(v4c.list, function(v4) {
    chr = sub(":.+", "", v4$prox [1])
    v4$prox <- (sub(".+:", "", v4$prox [1]) %>% as.numeric + shift) %>% paste0(chr, ":", .)
    v4$distal <- v4$distal + shift
    prox <- v4$prox [1]
    gene <- v4$gene[1]

    zero <- data.frame(chr = chr, left = 0, right = v4$distal[1], value = 0)
    bdg <- v4 %>% mutate(chr = chr, left = distal, right = distal + resolution) %>%
      dplyr::select(chr, left, right, value)
    bdg <- rbind(zero, bdg)

    bdg.name <- paste0(bdg.dir, "/", root.name, "_", genome, "_", prox, "_",resolution,"_", length,
     "_shift_",shift,  ".bdg")
    write.table(bdg, bdg.name, sep = "\t", quote = F, col.names = F, row.names = F)
    if (!is.null(liftover.params))
      bdg.name <- do.call(lift.over.core, c(liftover.params, list(in.file = bdg.name, out.file = paste0(bdg.name, ".lift"))))
    system(paste0("~/scripts/bed_browser_v2.sh ", bdg.name))

    if (!is.null(add.vline)) {
      bed <- as.data.frame(add.vline) %>% t()
      colnames(bed) <- c("left", "right")
      bed <- utilsFanc::add.column.fanc(df1 = bed, df2 = data.frame(chr = rep(chr, nrow(bed))), pos = 1)
      bed.name <- paste0(bdg.dir, "/", basename(hic.file), "_", prox, "_",resolution,"_", length, ".bed")
      write.table(bed, bed.name, sep = "\t", quote = F, col.names = F, row.names = F)
      if (!is.null(liftover.params))
        bed.name <- do.call(lift.over.core, c(liftover.params, list(in.file = bed.name, out.file = paste0(bed.name, ".lift"))))
      system(paste0("~/scripts/bed_browser_v2.sh ", bed.name))
    }
  })
}

lift.over.core <- function(in.file, out.file=NULL, in.genome, out.genome,
                           liftover = "/opt/apps/kentUCSC/334/liftOver", run = T) {
  if (!is.character(in.file)) {
    in.file.ori <- in.file
    in.file <- tempfile()
    write.table(in.file.ori, in.file, quote= F, col.names=F, row.names = F, sep = "\t")
  }
  out.file.ori <- out.file
  if (is.null(out.file)) {
    out.file <- tempfile()
  }

  chain <- paste0("/bar/cfan/genomes/liftover/", in.genome, "To", stringr::str_to_title(out.genome), ".over.chain.gz")
  cmd <- paste0(liftover, " ", in.file, " ", chain, " ", out.file, " /dev/null")
  print(cmd)
  if (run == T)
    system(cmd)
  if(!file.exists(out.file))
    stop(paste0(out.file, " was not successfully generated. liftover failed"))
  if (is.null(out.file.ori))
    out.file <- read.table(out.file, as.is=T, header = F)

  return(out.file)
}

V4C.loop <- function(hic.df, bdg.dir, length = 2000000, resolution = 5000, normalization = "VC_SQRT") {
  # fields required: file, genome, region, root.name, region.name, shift, use.1.6
  if (is.character(hic.df))
    hic.df <- read.table(hic.df, header = T)
  hic.df %>% split(., f = 1:nrow(hic.df)) %>% 
    lapply(function(x) {
      trash <- V4C.juicer(x$region, length = length, resolution=resolution, genome.build = x$genome, plot = F,
                                observed_or_oe = "observed", normalization = normalization, 
                          hic.file = x$file,
                                help.only = F, # add.vline=list(c(45279477, 45279842), c(44952775, 44957602)),
                                rm.zero=F, shift = x$shift, out.root.name = paste0(x$root.name, "_", x$region.name),
                                bdg.dir = bdg.dir, use.1.6 = x$use.1.6)
      return()
    })
  return()
}
