.hic.smart.name <- function(.hic.vec) {
  file.names <- sapply(.hic.vec, function(file) {
    if (grepl("aligned", file) && grepl("inter", file)) {
      res <- file %>% sub("//", "/", .) %>% gsub("/*mega", "", .) %>% sub("/*aligned/", "_", .) %>% 
        basename()
      return(res)
    } else {
      return(basename(file))
    }
  }) %>% `names<-`(NULL)
  return(file.names)
}

cm.add.dimnames <- function(cm, simplify = T) {
  if (simplify == T) {
    dimnames(cm@matrix) <- list(cm@regions[cm@anchor1] %>% start(),
                                cm@regions[cm@anchor2] %>% start())
  } else {
    dimnames(cm@matrix) <- list(cm@regions[cm@anchor1] %>% utilsFanc::gr.get.loci(),
                                cm@regions[cm@anchor2] %>% utilsFanc::gr.get.loci())
  }

  return(cm)
}

circle.gen <- function(anchor.1, anchor.2, bin.size, extend.size = NULL,
                       n = 1) {
  a.1 <- utilsFanc::gr.fit2bin(anchor.1)
  a.2 <- utilsFanc::gr.fit2bin(anchor.2)
  
  
}

gr.full.interaction <- function(gr, allow.self = F, strict = T) {
  # make a GInteractions object out of gr. Each region in gr interacts with each other
  n <- length(gr)
  mat <- outer(1:n, 1:n, FUN = function(x, y) paste0(x, ";", y))
  if (allow.self == F)
    diag(mat) <- NA
  vec <- mat %>% t() %>%  as.vector() %>% .[!is.na(.)]
  a1 <- vec %>% sub(";.+", "", .) %>% as.numeric()
  a2 <- vec %>% sub(".+;", "", .) %>% as.numeric()
  gi <- GInteractions(anchor1 = a1, anchor2 = a2, regions = gr)
  if (strict == T) {
    gi <- unique(swapAnchors(gi))
  }
  return(gi)
}

intSet.add.coverage <- function(intSet, assay = "NONE") {
  # intSet <- hico.list$inter$res_5K
  
}

norm.vec.gen.m <- function(.hic.vec.list, 
                           genome, chr, norms = c("VC", "VC_SQRT"), resolution.vec,
                           scale.to = 100000,
                           threads = NULL,
                           use.1.6 = T, 
                           print.cmd = F, 
                           java.path = "/opt/apps/java/jdk1.8.0_92/bin/java",
                           juicer.path.legacy = "/bar/cfan/software/juicer/SLURM/scripts/juicer_tools.jar",
                           juicer.path.1.6 = "~/software/juicer-1.6/SLURM/scripts/juicer_tools.jar") {
  if (is.null(names(.hic.vec.list))) {
    stop("is.null(names(.hic.vec.list))")
  }
  
  normo <- lapply(.hic.vec.list, function(.hic.vec) {
    normo <- lapply(resolution.vec, function(resolution) {
      normo <- norm.vec.gen(.hic.vec = .hic.vec, genome = genome, chr = chr, 
                            norms = norms, resolution = resolution, 
                            scale.to = scale.to, threads = threads,
                            use.1.6 = use.1.6, print.cmd = print.cmd,
                            java.path = java.path, juicer.path.legacy = juicer.path.legacy,
                            juicer.path.1.6 = juicer.path.1.6)
      metadata(normo)$resolution <- resolution
      return(normo)
    })
    names(normo) <- paste0("res_", resolution.vec)
    return(normo)
  })
  return(normo)
}

norm.vec.gen <- function(.hic.vec, genome, chr, norms = c("VC", "VC_SQRT"), resolution,
                         scale.to = 100000,
                         threads = NULL,
                         use.1.6 = T, 
                         print.cmd = F, 
                         java.path = "/opt/apps/java/jdk1.8.0_92/bin/java",
                         juicer.path.legacy = "/bar/cfan/software/juicer/SLURM/scripts/juicer_tools.jar",
                         juicer.path.1.6 = "~/software/juicer-1.6/SLURM/scripts/juicer_tools.jar") {
  if (is.null(threads))
    threads <- min(length(.hic.vec), 12)
  if (is.null(names(.hic.vec))) {
    names(.hic.vec) <- .hic.smart.name(.hic.vec = .hic.vec)
  }
  if (length(use.1.6) == 1) {
    use.1.6 <- rep(use.1.6, length(.hic.vec))
  }
  if (is.null(names(use.1.6)))
    names(use.1.6) <- names(.hic.vec)
  
  chr.length <- read.table(paste0("~/genomes/", genome, "/", genome, ".chrom.sizes")) %>% 
    filter(V1 == chr) %>% pull(V2)
  if (length(chr.length) > 1) {
    stop("length(chr.length) > 1")
  }
  names(chr.length) <- paste0(chr)
  tiles <- tileGenome(chr.length, tilewidth = resolution, cut.last.tile.in.chrom = T) 
  assay.list <- utilsFanc::safelapply(norms, function(norm) {
    mat <- utilsFanc::safelapply(.hic.vec %>% names(), function(sample) {
      .hic <- .hic.vec[sample]
      use.1.6 <- use.1.6[sample]
      if (use.1.6 == T) {
        juicer.path <- juicer.path.1.6
      } else {
        juicer.path <- juicer.path.legacy
      }
      tempt.file <- tempfile()
      cmd <- paste0(java.path, " -jar ", juicer.path, " dump norm ", norm, 
                    " ", .hic, " ", chr, " BP ", resolution, " ", tempt.file)
      if (print.cmd == T) {
        print(cmd); system(cmd)
      } else {
        invisible(system(cmd, intern = T))
      }
      res <- readLines(tempt.file) %>% as.numeric()
      return(res)
    }, threads = threads) %>% Reduce(cbind, .)
    colnames(mat) <- names(.hic.vec)
    return(mat)
  })
  names(assay.list) <- norms
  assay.list.cpm <- lapply(assay.list, function(mat) {
    mat <- scale.to * edgeR::cpm(mat)/1000000
    return(mat)
  })
  names(assay.list.cpm) <- paste0(names(assay.list), "..cpm")
  se <- SummarizedExperiment(assays = c(assay.list, assay.list.cpm), rowRanges = tiles, 
                             metadata = list(resolution = resolution, 
                                             scale.to = scale.to))
  return(se)
}

assay.fanc <- function(se, assay.name) {
  gr <- rowRanges(se)
  df <- assay(se, assay.name) %>% as.data.frame()
  mcols(gr) <- cbind(mcols(gr), df)
  return(gr)
}

mat.add.by.dimnames <- function(mat.list, dim.names = NULL, threads = 1, 
                                parse.region.ids = T, sep = "_",
                                col.region.order = NULL,
                                row.region.order = NULL,
                                add.transpose = F) {
  if (is.null(dim.names)) {
    if (parse.region.ids == T) {
      if (is.null(names(mat.list))) {
        stop("is.null(names(mat.list))")
      }
      mat.list <- mapply(function(mat, name) {
        name <- name %>% sub("\\.\\..+$", "", .)
        dim.names <- strsplit(name, split = sep) %>% unlist()
        rownames(mat) <- paste0(dim.names[1], "..", rownames(mat))
        colnames(mat) <- paste0(dim.names[2], "..", colnames(mat))
        return(mat)
      }, mat.list, names(mat.list), SIMPLIFY = F)
    }
    if (add.transpose == T) {
      mat.list <- c(mat.list, lapply(mat.list, Matrix::t))
      # note: during transposation, names(mat.list) is not changed accordingly,
      # you shouldn't parse names(mat.list) anymore!
    }
    dim.names <- lapply(1:2, function(i) {
      if (i == 1) {
        region.order <- row.region.order
      } else {
        region.order <- col.region.order
      }
      res <- utilsFanc::safelapply(mat.list, function(mat) return(dimnames(mat)[[i]]),
                                   threads = threads) %>% 
        unlist() %>% unique() %>% utilsFanc::mixed.sort.by(y = region.order, sep = "\\.\\.", 
                                                           sep.out = "..")
      return(res)
    })
    names(dim.names) <- c("row", "col")
  }
  mat.list.transformed <- utilsFanc::safelapply(mat.list, function(mat) return(mat.expand.transform(mat = mat, dim.names = dim.names )),
                                                threads = threads)
  res <- Reduce(`+`, mat.list.transformed)
  return(res)
}

mat.expand.transform <- function(mat, dim.names) {
  if (! identical(sort(names(dim.names)), c("col", "row"))) {
    stop("dim.names must be a named list. names must be col and row")
  }
  dim.names <- dim.names[c("row", "col")]
  dim.names <- lapply(dim.names, function(x) return(as.character(x)))
  
  if (is.data.frame(mat))
    mat <- as.matrix(mat)
  
  mat <- mat %>% as(., "dgCMatrix")
  if (! "dgCMatrix" %in% class(mat)) {
    stop("mat is was not successfully converted to a diCMatrix")
  }
  mat.sum <- Matrix::summary(mat)
  
  indices <- lapply(c("row", "col"), function(margin) {
    if (margin == "row") {
      utilsFanc::check.intersect(x = rownames(mat), x.name = "rownames(mat)",
                                 y = dim.names$row, y.name = "dim.names$row",
                                 n.examples = 3)
      pointer <- 1
      
    }
    if (margin == "col") {
      utilsFanc::check.intersect(x = colnames(mat), x.name = "colnames(mat)",
                                 y = dim.names$col, y.name = "dim.names$col",
                                 n.examples = 3)
      pointer <- 2
    }
    
    frame.df <- data.frame(frame.name = dim.names[[margin]])
    frame.df$frame.id <- 1:nrow(frame.df)
    
    mat.df <- data.frame(frame.name = dimnames(mat)[[pointer]])
    mat.df$mat.id <- 1:nrow(mat.df)
    
    mapping <- suppressMessages(left_join(frame.df, mat.df))
    if (nrow(mapping) != nrow(frame.df)) {
      stop("nrow(mapping) != nrow(frame.df)")
    }
    
    query <- data.frame(mat.id = mat.sum[[pointer]])
    mapped <- suppressMessages(left_join(query, mapping)$frame.id)
    return(mapped)
  })
  names(indices) <- c("row", "col")
  
  res <- Matrix::sparseMatrix(i = indices$row, j = indices$col, x = mat.sum$x, 
                              dims = sapply(dim.names, length),
                              dimnames = dim.names)
  return(res)
}

gr.slide.2 <- function(gr, name.col, bin.size,
                       n.bins = NULL, n.bins.up = NULL, n.bins.down = NULL,
                       return.assign.vec = F) {
  if (!is.null(n.bins)) {
    n.bins.up <- n.bins
    n.bins.down <- n.bins
  }
  shift.vec <- 0
  if (!is.null(n.bins.up)) {
    shift.vec <- c(-1 * (n.bins.up:1), shift.vec)
  }
  if (!is.null(n.bins.down)) {
    shift.vec <- c(shift.vec, 1:n.bins.down)
  }
  shift.vec <- shift.vec * bin.size
  
  gr.shifted.list <- lapply(shift.vec, function(shift) {
    return(GenomicRanges::shift(gr, shift = shift))
  })
  gr.shifted <- Reduce(c, gr.shifted.list)
  if (return.assign.vec == T) {
    assign.df <- gr.shifted %>% as.data.frame() %>% 
      dplyr::select(start, !!as.name(name.col)) %>% 
      arrange(start) %>% 
      mutate(start = as.character(start))
    assign.vec <- assign.df[, name.col]
    names(assign.vec) <- assign.df$start
    return(assign.vec)
  }
  return(GenomicRanges::sort(gr.shifted))
}

gi.get.active.regions <- function(gi) {
  indices <- anchorIds(gi) %>% unlist() %>% unique() %>% sort()
  res <- regions(gi)[indices]
  return(res)
}

gr.use.ly49 <- function(gr, name.col = "forth") {
  mcols(gr)[, name.col][mcols(gr)[, name.col] == "Gm15854"] <- "Klra24"
  mcols(gr)[, name.col] <- mcols(gr)[, name.col] %>% gsub("[a-zA-Z]","",.) %>% 
    sub("\\-", "",.) %>% as.numeric() %>% letters[.] %>% paste0("Ly49", .) 
  return(gr)
}

make.unique.reduce.complexity <- function(vec) {
  df <- data.frame(ori = vec, 
                   root = sub("\\..+$", "", vec))
  df$suffix <- 0
  df$suffix[grepl("\\.\\d+$", df$ori)] <- sub(".+\\.", "", df$ori[grepl("\\.\\d+$", df$ori)]) %>% 
    as.numeric()
  df <- df %>% group_by(root) %>% summarise(min.suffix = min(suffix)) %>% ungroup() %>% 
    as.data.frame()
  df$min.suffix[df$min.suffix > 1] <- 1
  df$min.suffix[df$min.suffix == 0] <- ""
  out <- paste0(df$root, ".", df$min.suffix) %>% sub("\\.$","", .)
  return(out)
}



bw.pileup.browser <- function(gr, bw, out.file, output.each = F,
                              scale.to = NULL, take.mean = T, one.by.one = F,
                              chr = "chr6", start = 130120000) {
  if (length(unique(width(gr))) != 1) {
    stop("length(unique(width(gr))")
  }
  w <- width(gr)[1]
  frame <- data.frame(chr = chr, start = start, end = start + w -1) %>% 
    makeGRangesFromDataFrame()
  dir.create(dirname(out.file), recursive = T, showWarnings = F)
  signal <- lapply(seq_along(gr), function(i) {
    region <- gr[i]
    s <- rtracklayer::import(con = bw, which = region)
    if (!is.null(scale.to))
      s$score <- scale.to * s$score/sum(s$score)
    s <- as.data.frame(s) %>% mutate(seqnames = chr) %>% 
      makeGRangesFromDataFrame(keep.extra.columns = T)
    s <- shift(s, -1 * start(region))
    if (one.by.one == T) {
      start <- start + w * (i-1)
    }
    s <- shift(s, start)
    if (output.each == T) {
      out.file <- out.file %>% utilsFanc::insert.name.before.ext(i)
      rtracklayer::export(object = cov, con = out.file)
      system(paste0("~/scripts/bb.sh ", out.file))
    }
    return(s)
  })
  signal <- Reduce(c, signal)
  if (take.mean == T)
    signal$score <- signal$score/length(gr)
  cov <- coverage(signal, weight = "score") %>% as("GRanges")
  rtracklayer::export(object = cov, con = out.file)
  if (!grepl(".bw$", out.file)) {
    system(paste0("~/scripts/bb.sh ", out.file))
    out.file <- paste0(out.file, ".gz")
  }
  cat(utilsFanc::bash2ftp(filename = out.file))
  cat("\n")
  return(out.file)
}

gi.area <- function(gi, bin.size) {
  if (length(gi) > 1) {
    stop("only handle 1 interaction")
  }
  regions(gi) <- utilsFanc::gr.fit2bin(regions(gi), bin.size = bin.size)
  area <- (width(anchors(gi, "first"))/bin.size) * (width(anchors(gi, "second"))/bin.size)
  return(area)
}

formula.parse.fanc <- function(formula, df) {
  # format: value.col~colname.1:value1,value2,...,valuen;colname.n:value.1,value.n@colname.m:value.1,value.z
  # value.col: where numbers will be taken from.
  # to the left of @: x; to the right of @: y
  # if you want to use paired stats, add "pair.col>" before everything else
  if (grepl(">", formula)) {
    pair.col <- sub(">.+", "", formula)
    formula <- sub(".+>", "", formula)
    if (!pair.col %in% colnames(df))
      stop(paste0("pair.col: ", pair.col, " not found"))
  } else {
    pair.col <- NULL
  }
  value.col <- sub("~.+", "", formula)
  if (!value.col %in% colnames(df))
    stop(paste0("value.col: ", value.col, " not found"))
  res <- lapply(c("\\1", "\\2"), function(i) {
    filters <- sub(".+~(.+)@(.+)", i, formula) %>% strsplit(split=";") %>% unlist()
    if (length(filters) < 1) {
      stop("length(filters) < 1" %>% paste0(" for ", i))
    }
    for (filter in filters) {
      if (!grepl(":", filter)) {
        stop(paste0('!grepl(":", filter)', " for filter ", filter))
      }
      col <- sub(":.+", "", filter)
      include <- sub(".+:", "", filter)
      df <- df %>% .[.[, col] %in% include,]
    }
    if (!is.null(pair.col)) {
      df <- df[order(df[,pair.col]), ]
    }
    r <- df[, value.col]
    return(r)
  })
  names(res) <- c("x", "y")
  return(res)
}

stats.wrapper.ez <- function(vec.list, stat.fun, stat.params = NULL, slot = NULL) {
  # vec.list: 
  # $x
  # [1] 0.004307084 0.004511637
  # $y
  # [1] 0.002739732 0.002775481
  params <- vec.list
  names(params) <- NULL
  if (!is.null(stat.params)) {
    params <- c(params, stat.params)
  }
  res <- do.call(stat.fun, params)
  if (!is.null(slot)) {
    res <- res[[slot]]
  }
  return(res)
}

normo.2.bdg <- function(normo, assay = "VC..cpm", out.dir) {
  dir.create(out.dir, showWarnings = F, recursive = T)
  lapply(names(normo), function(mapq) {
    normo.mapq <- normo[[mapq]]
    lapply(names(normo.mapq), function(res) {
      normo.res <- normo.mapq[[res]]
      if (! "mean" %in% colnames(normo.res)) {
        mat <- assay(normo.res, assay) 
        mean <- mat %>% rowMeans()
        mat <- cbind(mat, mean)
        lapply(colnames(mat), function(sample) {
          bdg <- rowRanges(normo.res)
          bdg$score <- mat[, sample]
          bdg.file <- paste0(out.dir, "/", mapq, "..", res, "..", sample, ".bedGraph")
          rtracklayer::export.bedGraph(object = bdg, bdg.file)
          system(paste0("~/scripts/bb.sh ", bdg.file))
          return()
        })
      }
    })
  })
}