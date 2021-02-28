contact.stats <- function(prox, distal.vec.fg, distal.vec.bg=NULL, int.names=NULL, 
                          sample.info.tsv, 
                          resolution, genome, v4.length, normalization,
                          collapse.rep = F, pseudo.count = 0.01, threads = 6,
                          plot.out = NULL, return.df = F) {
  sample.info <- read.table(sample.info.tsv, as.is = T, sep = "\t", quote = "", header = T)
  
  fg.df <- distal.vec.fg %>% region.2.bins(include.left = T, include.right = T, resolution = resolution, 
                                           genome = genome, return.list = F, out.bed = NULL)
  if (!is.null(distal.vec.bg)) {
    if (distal.vec.bg == "flip") {
      prox.bin <- paste0(prox, as.numeric(sub("chr.+:", "-", prox)) + resolution -1)
      distal.vec.bg <- flip.bin(anchor = prox.bin, distal.vec = distal.vec.fg)
    }
    if (length(distal.vec.bg) == 1) {
      distal.vec.bg <- rep(distal.vec.bg, length(distal.vec.fg))
    }
    fg.2.bg <- distal.vec.bg
    names(fg.2.bg) <- distal.vec.fg
    bg.2.fg <- names(fg.2.bg)
    names(bg.2.fg) <- fg.2.bg
    bg.df <- distal.vec.bg %>% region.2.bins(include.left = T, include.right = T, resolution = resolution, 
                                             genome = genome, return.list = F, out.bed = NULL)
    bg.df$ori <- bg.2.fg[bg.df$ori]
    
  }
  
  int.df <- sample.info %>% split(f = 1:nrow(sample.info)) %>% 
    mclapply(function(s) {
      v4.df <- V4C.juicer( prox, length = v4.length, resolution=resolution, genome.build = genome, plot = F,
                           observed_or_oe = "observed", normalization = normalization, hic.file = s$hic.file,
                           help.only = F, use.1.6 = s$use.1.6, strip.chr = s$strip.chr,
                           rm.zero=F)
      v4.df <- v4.df[[1]][,c("prox", "distal", "value")]
      fg <- fg.df[,c("left", "ori")]
      colnames(fg) <- c("distal", "ori")
      int.df <- left_join(fg, v4.df)
      int.df <- int.df %>% group_by(ori) %>% summarise(int.freq.fg = sum(value))
      
      if (!is.null(distal.vec.bg)) {
        bg <- bg.df[,c("left", "ori")]
        colnames(bg) <- c("distal", "ori")
        int.bg <- left_join(bg, v4.df)
        int.bg <- int.bg %>% group_by(ori) %>% summarise(int.freq.bg = sum(value))
        int.df <- left_join(int.df, int.bg)
        int.df$bg.ori <- fg.2.bg[int.df$ori]
      } else {
        int.df$int.freq.bg <- 1
        int.df$bg.ori <- "no_bg"
      }
      int.df <- cbind(int.df, s[,c("sample.type", "rep", "hic.file")])
      return(int.df)
    }, mc.cores = threads, mc.cleanup = T) %>% Reduce(rbind,.)
  
  # browser()
  
  if (!is.null(int.names)) {
    names(int.names) <- distal.vec.fg
    int.df$int.name <- int.names[int.df$ori]
  } else {
    int.df$int.name <- paste0(int.df$ori, "/", int.df$bg.ori)
  }
  
  
  if (collapse.rep == T) {
    int.df <- int.df %>% group_by(int.name, sample.type, rep) %>% 
      summarise(int.freq.fg = sum(int.freq.fg),
                int.freq.bg = sum(int.freq.bg))
    if (is.null(distal.vec.bg))
      int.df$int.freq.bg <- 1
  }
  
  int.df$int <- int.df$int.freq.fg/(int.df$int.freq.bg + pseudo.count)

  p <- ggbarplot(int.df, x = "int.name", y = "int", color = "sample.type", position = position_dodge(),
            add = c("mean_se", "dotplot"))
  if (!is.null(plot.out)) {
    system(paste0("mkdir -p ", dirname(plot.out)))
    ggsave(plot.out, p, width = 6, height = 4, device = "png", units = "in", dpi = 150)
  }
  if (return.df == T)
    return(int.df)
  return(p)
}


flip.bin <- function(anchor.vec, distal.vec) {
  # anchor.vec format: chr1:5000-9999
  if (length(anchor.vec == 1))
    anchor.vec <- rep(anchor.vec, length(distal.vec))
  anchors <- get.pos.from.string(anchor.vec)
  distal <- get.pos.from.string(distal.vec)
  flipped <- sapply(seq_along(distal), function(i) {
    x <- distal[[i]]
    anchor <- anchors[[i]]
    dist <- x$left - anchor$right
    length <- x$right - x$left
    new <- x
    new$right <- anchor$left - dist
    new$left <- new$right - length
    new <- paste0(new$chr, ":", new$left, "-", new$right)
    return(new)
  })
  return(flipped)
}

get.pos.from.string <- function(pos, reduce = F) {
  # pos = c("chr1:5000-9999", "chr1:5000-9999", "chr1:5000-9999")
  res <- strsplit(pos, ":|-")
  res <- lapply(res, function(x) {
    l <- as.list(x)
    names(l) <- c("chr", "left", "right")
    l[2:3] <- as.numeric(l[2:3])
    return(l)
  })
  if (length(res) == 1 && reduce == T)
    res <- res[[1]]
  return(res)
}

region.2.bins <- function(region.vec, include.left = T, include.right = T, resolution, genome,
                          out.bed = NULL, return.list = F) {
  q.df <- data.frame(chr = sub(":.+$", "",region.vec),
                   left = sub("chr.+:", "",region.vec) %>% sub("-\\d+$", "", .) %>% as.numeric(),
                   right = sub("^.+-", "", region.vec) %>% as.numeric())
  bin.df <- read.table(paste0("/bar/cfan/genomes/", genome, "/bins/", genome, "_", resolution, ".bins"),
                       as.is = T, sep = "\t", quote = "")
  colnames(bin.df) <- c("chr", "left", "right")
  out.list <- q.df %>% split(f = 1:nrow(q.df)) %>% 
    lapply(function(x) {
      bin.df <- bin.df %>% filter(chr == x$chr) %>% mutate(., index = 1:nrow(.))
      index.1 <- bin.df %>% filter(left <= x$left, right >= x$left) %>% pull(index)
      if (include.left == F && bin.df[index.1, "left"] < x$left)
        index.1 <- index.1 + 1
      index.2 <- bin.df %>% filter(left <= x$right, right >= x$right) %>% pull(index)
      if (include.right == F && bin.df[index.2, "right"] > x$right)
        index.2 <- index.2 - 1
      bins <- bin.df[index.1:index.2,]
      bins$ori <- paste0(x$chr, ":", x$left,"-", x$right)
      bins$index <- NULL
      return(bins)
    })
  names(out.list) <- sapply(out.list, function(x)return(x$ori[1]))
  if (return.list == T)
    return(out.list)
  out.df <- Reduce(rbind, out.list)
  rownames(out.df) <- NULL
  if (!is.null(out.bed)) {
    system(paste0("mkdir -p ", dirname(out.bed)))
    write.table(out.df, out.bed, sep = "\t", row.names = F, col.names = F, quote = F)
    system(paste0("/bar/cfan/scripts/bed_browser_v2.sh ", out.bed))
  }
  return(out.df)
}

