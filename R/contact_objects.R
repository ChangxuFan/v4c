read.hic.core <- function(.hic, region1, region2 = NULL,
                          observed.or.oe, normalization, resolution,
                          juicer.out.file = NULL,
                          java.path = "/opt/apps/java/jdk1.8.0_92/bin/java",
                          use.1.6 = T, return.mat = F, return.frame = F, print.cmd = F,
                          juicer.path.legacy = "/bar/cfan/software/juicer/SLURM/scripts/juicer_tools.jar",
                          juicer.path.1.6 = "~/software/juicer-1.6/SLURM/scripts/juicer_tools.jar") {
  # note: .hic could be a url
  if (use.1.6 == T) {
    juicer.path <- juicer.path.1.6
  } else {
    juicer.path <- juicer.path.legacy
  }

  if (is.null(region2)) {
    region2 <- region1
  } else {
    stop("currently region2 has to equal region1")
  }

  regions.in <- list(r1 = region1, r2 = region2)
  rm(region1); rm(region2)
  regions.in <- lapply(regions.in, function(r) {
    if (length(r) != 1) {
      stop("only single regions allowed in read.hic.core")
    }
    if ("GRanges" %in% class(r)) {
      res <- utilsFanc::gr.get.loci(gr = r)
    } else if (is.character(r)) {
      res <- r
    } else {
      stop("read.hic.core: only characters and Granges allowed as region input")
    }
    res <- res %>% sub("-", ":", .)
    return(res)
  })

  if (is.null(juicer.out.file))
    juicer.out.file <- tempfile()

  if (return.frame == F) {
    cmd <- paste0(java.path, " -jar ", juicer.path, " dump ", observed.or.oe,
                  " ", normalization, " ", .hic,
                  " ", regions.in$r1, " ", regions.in$r2, " BP ", resolution, " ", juicer.out.file)
    if (print.cmd == T) {
      print(cmd); system(cmd, intern = T)
    } else {
      invisible(system(cmd, intern = T))
    }
  }

  regions.out <- lapply(regions.in, function(r) {
    r.start.given <- sub("(.+):(.+):(.+)", "\\2", r) %>% as.numeric()
    r.start <- resolution * ceiling(r.start.given/resolution)
    r.end.given <- sub("(.+):(.+):(.+)", "\\3", r) %>% as.numeric()
    r.end <- resolution * floor(r.end.given/resolution)
    n.bins <- (r.end - r.start)/resolution + 1
    dim.name <- (0:(n.bins-1))*resolution + r.start
    return(dim.name)
  })
  if (return.frame == F) {
    df <- read.table(juicer.out.file, as.is = T)
    # a note a juicer dump's behaviour: it only reports the interaction between r1 and r2.
    # interactions with both anchors within r1 or r2 are not reported.
    # examples:
    # Browse[1]> df[! df$V1  %in% regions.out$r2 & ! df$V2  %in% regions.out$r2, ] %>% nrow()
    # [1] 0
    # Browse[1]> nrow(df)
    # [1] 267
    # Browse[1]> df[ df$V1  %in% regions.out$r2 |  df$V2  %in% regions.out$r2, ] %>% nrow()
    # [1] 267
    i <- (df$V1 - regions.out$r1[1])/resolution +1
    j <- (df$V2 - regions.out$r2[1])/resolution +1
    x <- df$V3
  } else {
    i <- 1; j <- 1; x <- 1
  }
  mat <- Matrix::sparseMatrix(i = i, j = j, x=x,
                              dims = sapply(regions.out, length),
                              dimnames = regions.out %>% `names<-`(NULL))
  mat <- mat + Matrix::t(mat) - Matrix::diag(Matrix::diag(mat))
  if (return.mat == T) {
    return(mat)
  }
  grs.out <- lapply(c("r1", "r2"), function(r) {
    chr <- regions.in[[r]] %>% gsub(":.+$", "",.)
    gr <- paste0(chr, ":", regions.out[[r]] + 1, "-", regions.out[[r]] + resolution) %>%
      utilsFanc::loci.2.df(loci.vec = ., return.gr = T)
    return(gr)
  }) %>% `names<-`(c("r1", "r2"))
  cm <- ContactMatrix(matrix = mat, anchor1 = grs.out$r1, anchor2 = grs.out$r2)
  if (return.frame == T) {
    int <- deflate(x = cm, use.zero = T)
    return(int@interactions)
  }
  return(cm)
}

read.hic <- function(.hic.vec,
                     region1, region2 = NULL,
                     observed.or.oe = "observed",
                     normalization = c("NONE", "VC_SQRT", "VC"),
                     resolution = c(1000, 5000, 10000, 25000),
                     use.1.6 = T, print.cmd = F,
                     threads.norm = 1, threads.res = 1, threads.sample = 1,
                     java.path = "/opt/apps/java/jdk1.8.0_92/bin/java",
                     juicer.path.legacy = "/bar/cfan/software/juicer/SLURM/scripts/juicer_tools.jar",
                     juicer.path.1.6 = "~/software/juicer-1.6/SLURM/scripts/juicer_tools.jar",
                     ...) {
  if (is.null(names(.hic.vec))) {
    names(.hic.vec) <- .hic.smart.name(.hic.vec = .hic.vec)
  }
  if (is.null(threads.norm))
    threads.norm <- length(normalization)
  if (is.null(threads.res))
    threads.res <- length(resolution)
  if (is.null(threads.sample))
    threads.sample <- length(.hic.vec)
  intSets <- utilsFanc::safelapply(resolution, function(res) {
    # prefix.res <- utilsFanc::so.formatter(x = res)
    mat.list <- utilsFanc::safelapply(normalization, function(norm) {
      # prefix.norm <- norm
      intSet <- utilsFanc::safelapply(.hic.vec, function(.hic) {
        # if (grepl("NK1.+neg", .hic)) {
        #   browser()
        # }
        cm <- read.hic.core(.hic = .hic, observed.or.oe = observed.or.oe, normalization = norm,
                            region1 = region1, region2 = region2,
                            resolution = res, return.mat = F,
                            java.path = java.path, print.cmd = print.cmd,
                            use.1.6 = use.1.6, juicer.path.legacy = juicer.path.legacy,
                            juicer.path.1.6 = juicer.path.1.6)
        intSet <- deflate(x = cm, use.zero = T)
        return(intSet)
      }, threads = threads.sample) %>% Reduce(cbind, .)
      mat <- assay(intSet, drop = F)
      colnames(mat) <- names(.hic.vec)
      return(mat)
    }, threads = threads.norm)
    names(mat.list) <- normalization
    gi <- read.hic.core(region1 = region1, region2 = region2, return.frame = T, resolution = res)
    intSet <- InteractionSet(assays = mat.list, interactions = gi, ...)
    metadata(intSet)$observed.or.oe <- observed.or.oe
    metadata(intSet)$resolution <- res
    return(intSet)
  }, threads = threads.res)
  names(intSets) <- paste0("res_", resolution)
  return(intSets)
}

contact.summary <- function(ins, gi, gi.is.list = F, gi.name.col,
                            assay, samples = NULL,
                            sum.fun = sum,
                            do.donut = F, donut.center, donut.peri,
                            donut.name = NULL, overwrite.peri = T,
                            bin.size,
                            coldata = NULL, stat.formula = NULL, stat.fun, stat.params = NULL) {
  if (gi.is.list == T) {
    stop("not developed yet")
  }
  if (!is.null(samples)) {
    bSamples <- colnames(ins) %in% samples
    if (sum(bSamples) < 1) {
      stop("bSamples < 1")
    }
    ins <- ins[, bSamples]
  }
  samples <- colnames(ins)
  gi.list <- lapply(1:length(gi), function(i) return(gi[i]))

  names(gi.list) <- sapply(gi, function(x) return(mcols(x)[, gi.name.col]))
  int.stat.list <- lapply(gi.list, function(gi) {
    int.name <- mcols(gi)[, gi.name.col]
    ins.sub <- subsetByOverlaps(x = ins, ranges = gi, ignore.strand = T)
    int.stat <- assays(ins.sub)[[assay]] %>% apply(2, FUN = sum.fun) %>%
      as.matrix() %>% t()
    rownames(int.stat) <- int.name
    expanded <- lapply(samples, function(sample) {
      cm <- inflate(ins.sub, rows = anchors(gi)$first, columns = anchors(gi)$second,
                    assay = assay, sample = sample, sparse = T) %>% cm.add.dimnames()
      return(cm)
    })
    names(expanded) <- samples
    res <- list(int.name = int.name, int.stat = int.stat,
                ins.sub = ins.sub, gi = gi, expanded = expanded)
    return(res)
  })
  int.stat <- lapply(int.stat.list, function(x) {
    return(x$int.stat)
  }) %>% Reduce(rbind, .)

  if (do.donut == T) {
    if (identical(all.equal(sum.fun, mean), TRUE)) {
      areas <- lapply(gi.list, gi.area, bin.size = bin.size)
      donut <- (int.stat[donut.peri, ] * areas[[donut.peri]] - int.stat[donut.center, ] * areas[[donut.center]]) /
        (areas[[donut.peri]] - areas[[donut.center]])
    } else {
      donut <- int.stat[donut.peri, ] - int.stat[donut.center, ]
    }
    donut <- as.matrix(donut) %>% t()
    if (is.null(donut.name))
      donut.name <- donut.peri %>% paste0("..donut")
    rownames(donut) <- donut.name
    int.stat <- rbind(int.stat, donut)
    donut.hole.ratio <- int.stat[donut.center, ]/donut
    rownames(donut.hole.ratio) <- "donut.hole.ratio"
    int.stat <- rbind(int.stat, donut.hole.ratio)
    if (overwrite.peri == T)
      int.stat <- int.stat %>% .[! rownames(.) == donut.peri, ]
  }
  int.stat.melt <- int.stat %>% as.data.frame() %>% mutate(., int.name = rownames(.)) %>%
    reshape2::melt(id.vars = "int.name", variable.name = "sample",
                   value.name = "int.freq")
  if (!is.null(coldata)) {
    if (! "sample" %in% colnames(coldata))
      coldata$sample <- rownames(coldata)
    int.stat.melt <- left_join(int.stat.melt, coldata)
  }

  if (!is.null(stat.formula)) {
    p.df <- lapply(stat.formula, function(formula) {
      vec.list <- formula.parse.fanc(formula = formula, df = int.stat.melt)
      p <- stats.wrapper.ez(vec.list = vec.list, stat.fun = stat.fun, stat.params = stat.params, slot = "p.value")
      res <- data.frame(comp = formula, p.value = p)
      return(res)
    }) %>% Reduce(rbind, .)
  } else {
    p.df <- NULL
  }
  res <- list(stat = int.stat, stat.melt = int.stat.melt,
              p.df = p.df,
              details = int.stat.list)
  return(res)
}

contact.plot.bar <- function(df, x, y, fill, sum.to.1 = F, filter.col = fill,
                             include = NULL, exclude = "donut.hole.ratio",
                             plot.out = NULL, title = NULL) {
  if (!is.null(filter.col)) {
    if (!is.null(include)) {
      df <- df[df[, filter.col] %in% include,]
    }
    if (!is.null(exclude)) {
      df <- df[!df[, filter.col] %in% exclude,]
    }
  }

  if (sum.to.1 == T) {
    df <- df %>% group_by(!!as.name(x)) %>%
      mutate(!!y := !!as.name(y)/sum(!!as.name(y))) %>%
      ungroup()
  }
  df[, y] <- df[, y] %>% round(digits = 5)
  p <- ggpubr::ggbarplot(df, x = x, y = y, fill = fill, label = T, color = fill, alpha = 0.75,
                    lab.pos = "in", title = title) +
    theme(aspect.ratio = 1)

  if (!is.null(plot.out)) {
    dir.create(dirname(plot.out), showWarnings = F, recursive = T)
    ggsave(plot.out, p, dpi = 100, width = 6 , height = 5, units = "in")
  }
  return(p)
}

contact.pipe <- function(anchor1, anchor2, gi = NULL, gi.name.col,
                         intSet, assay, samples = NULL, sum.fun = sum,
                         resolution,
                         sum.to.1 = T,
                         out.dir, root.name = NULL,
                         write.anchors = T, do.plot = T,
                         do.donut = F, donut.center, donut.peri,
                         coldata = NULL,
                         donut.name = NULL, overwrite.peri = T, bin.size,
                         stat.formula = NULL, stat.fun, stat.params = NULL) {
  if (!is.null(coldata)) {
    coldata$sample <- rownames(coldata)
  }
  if (is.null(root.name)) {
    root.name <- basename(out.dir)
  }
  system(paste0("mkdir -p ", out.dir))
  if (is.null(gi)) {
    if (length(anchor1) == 1 && length(anchor2) > 1) {
      anchor1 <- rep(list(anchor1), length(anchor2)) %>%
        Reduce(c, .)
    }
    gi <- GInteractions(anchor1, anchor2)
  }
  if (write.anchors == T) {
    utilsFanc::write.zip.fanc(anchors(gi, "first") %>% utilsFanc::gr.fit2bin(bin.size = resolution, expand = T),
                              out.file = paste0(out.dir, "/", root.name, "_anchor1.bed"))
    utilsFanc::write.zip.fanc(anchors(gi, "second") %>% utilsFanc::gr.fit2bin(bin.size = resolution, expand = T),
                              out.file = paste0(out.dir, "/", root.name, "_anchor2.bed"))
  }
  contact.stat <- contact.summary(ins = intSet, gi = gi, gi.is.list = F,
                          gi.name.col = gi.name.col, assay = assay,
                          samples = samples, sum.fun = sum.fun,
                          do.donut = do.donut, donut.center = donut.center,
                          donut.peri = donut.peri, donut.name = donut.name,
                          overwrite.peri = overwrite.peri, bin.size = bin.size,
                          coldata = coldata,
                          stat.formula = stat.formula, stat.fun = stat.fun,
                          stat.params = stat.params)

  if (!is.null(contact.stat$p.df)) {
    write.table(contact.stat$p.df, paste0(out.dir, "/", root.name, "_pvalues.tsv"),
                sep = "\t", quote = F, row.names = F, col.names = T)
  }
  if (do.plot == T) {
    p1 <- contact.plot.bar(df = contact.stat$stat.melt,
                           x = "sample", y = "int.freq", fill = "int.name",
                           sum.to.1 = sum.to.1, title = root.name,
                           plot.out = paste0(out.dir, "/", root.name, "_bar.pdf"))

    if (!is.null(coldata)) {
      contact.stat$stat.melt <- left_join(contact.stat$stat.melt, coldata) %>% filter(int.name != "donut.hole.ratio")
      p1.1 <- ggbarplot(contact.stat$stat.melt, x = "Ly49D", y = "int.freq", fill = "int.name",
                        position = position_dodge(), color = "black", alpha = 0.5, add = "mean_sd",
                        add.params = list(size = 1), size = 1) +
        geom_point(aes(fill = int.name, shape = bio.rep), position = position_dodge(width = 0.8),
                   size = 3.0, color = "grey50") +
        theme_pubr(base_size = 18, legend = "right") +
        theme(aspect.ratio = 1) +
        ggtitle(root.name) +
        ggsave(paste0(out.dir, "/", root.name, "_bar_coldata.pdf"), height = 4, width = 8,
               dpi = 100)
    } else {
      p1.1 <- NULL
    }
    if (sum.to.1 == T) {
      p2 <- contact.plot.bar(df = contact.stat$stat.melt,
                             x = "sample", y = "int.freq", fill = "int.name",
                             sum.to.1 = F, title = root.name,
                             plot.out = paste0(out.dir, "/", root.name, "_bar_ori.pdf"))
    } else {
      p2 <- NULL
    }

    if (do.donut == T && !is.null(coldata)) {
      if (!"sample" %in% col(coldata)) {
        coldata$sample <- rownames(coldata)
      }
      pd.df <- contact.stat$stat.melt %>% left_join(coldata) %>% filter(int.name == "donut.hole.ratio")
      # if (identical(all.equal(sum.fun, mean), T)) {
      #   pd.df$grouping <- paste0(pd.df$int.name, "..", pd.df$bio.rep)
      #   pd <- pd.df %>% ggplot(aes(x = Ly49D, y = int.freq, group = grouping, color = int.name)) +
      #     geom_line(aes(linetype = bio.rep)) +
      #     geom_point() +
      #     theme_pubr(legend = "right") +
      #     theme(aspect.ratio = 1.75) +
      #     ggsave(paste0(out.dir, "/", root.name, "_donut_cmp.pdf"), height = 4, width = 4,
      #            dpi = 100)
      # } else {
      #   pd <- ggline(data = pd.df, x = "Ly49D", y = "int.freq", shape = "bio.rep",
      #                color = "bio.rep") +
      #     ggsave(paste0(out.dir, "/", root.name, "_donut_ratio_line.pdf"), height = 4, width = 3,
      #            dpi = 100)
      # }
      pd <- ggline(data = pd.df, x = "Ly49D", y = "int.freq", shape = "bio.rep",
                   color = "bio.rep") +
        ggsave(paste0(out.dir, "/", root.name, "_donut_ratio_line.pdf"), height = 4, width = 3,
               dpi = 100)

    } else {
      pd <- NULL
    }
  } else {
    p1 <- p1.1 <- p2 <- pd <- NULL
  }
  p.list <- list(p1 = p1, p1.1 = p1.1, p2 = p2, pd = pd)

  res <- list(contact.stat = contact.stat, gi = gi, p.list = p.list)

  saveRDS(res, paste0(out.dir, "/", root.name, ".Rds"))
  return(res)

}

contact.pipe.m <- function(hico.by.mapq, norms,
                           anchor1 = NULL, anchor2 = NULL, gi = NULL, gi.name.col,
                           sum.fun, resolution,
                           do.donut = F, donut.center, donut.peri,
                           donut.name = NULL, overwrite.peri = T, coldata = NULL,
                           out.dir,
                           stat.formula = NULL, stat.fun, stat.params = NULL,
                           ...) {
  if (identical(all.equal(sum.fun, sum), TRUE)) {
    sum.to.1 <- T
  } else if (identical(all.equal(sum.fun, mean), TRUE)) {
    sum.to.1 <- F
  } else {
    stop("only accept sum or mean for sum.fun")
  }
  grid <-  lapply(names(hico.by.mapq), function(mapq) {
    # don't forget hico is used here as well!!
    hico <- hico.by.mapq[[mapq]]
    l <- lapply(norms, function(norm) {
      # c.o: contact object
      out.dir <- paste0(out.dir, "/", mapq, "..", norm, "/")
      if (!is.null(coldata)) {
        rownames(coldata) <- paste0(rownames(coldata), "_", mapq) %>%
          sub("inter30", "inter_30", .)
        if (!grepl(".hic$", rownames(coldata)[1])) {
          rownames(coldata) <- paste0(rownames(coldata), ".hic")
        }
      }

      c.o <- contact.pipe(anchor1 = anchor1, anchor2 = anchor2, gi = gi,
                          gi.name.col = gi.name.col,
                          intSet = hico[[paste0("res_", resolution)]],
                          assay = norm, sum.fun = sum.fun,
                          out.dir = out.dir, sum.to.1 = sum.to.1,
                          do.donut = do.donut, donut.center = donut.center, donut.peri = donut.peri,
                          overwrite.peri = overwrite.peri, coldata = coldata, bin.size = resolution,
                          write.anchors = F,
                          stat.formula = stat.formula, stat.fun = stat.fun,
                          stat.params = stat.params, ...)
      return(c.o)
    })
    names(l) <- paste0(mapq, "..", norms)
    return(l)
  }) %>% Reduce(c, .)
  saveRDS(grid, paste0(out.dir, "/grid.Rds"))
  lapply(c("p1", "p1.1", "p2", "pd"), function(p.name) {
    pl <- lapply(grid, function(x) return(x$p.list[[p.name]]))
    if (p.name == "p1.1") {
      width <- 8; height <- 4
    } else {
      width <- height <- 4
    }
    p <- scFanc::wrap.plots.fanc(plot.list = pl,
                                 plot.out = paste0(out.dir, "/summary_", p.name, ".png"),
                                 sub.height = height, sub.width = width)

    return(NULL)
  })
  if (!is.null(stat.formula)) {
    p.df <- lapply(names(grid), function(name) {
      x <- grid[[name]]
      df <- x$contact.stat$p.df %>% mutate(id = name)
      return(df)
    }) %>% Reduce(rbind, .)
    write.table(p.df, paste0(out.dir, "/pvalue.tsv"),
                sep = "\t", quote = F, row.names = F, col.names = T)
  }
  return(grid)
}

contact.map.viz <- function(mat, color.range = NULL,
                            max.quantile = NULL, color = "magenta4",
                            zero.as.min = T, min.quantile = NULL,
                            row.assign.vec = NULL, col.assign.vec = NULL,
                            cluster_rows = F, cluster_columns = F,
                            clustering_distance_rows = "euclidean",
                            clustering_distance_columns = "euclidean",
                            add.text = T, fontsize = 10,
                            digits = 3, magnify = NULL,
                            width = 700, height = 750, res = 100,
                            plot.out = NULL) {
  # column and row assign vecs: named vector. Names are col/row names
  mat <- as.matrix(mat)
  if (is.null(magnify))
    magnify <- 1
  mat <- mat * magnify
  if (is.null(color.range)) {
    if (is.null(max.quantile))
      max.quantile <- 1
    if (is.null(min.quantile))
      min.quantile <- 0

    color.max <- quantile(mat, max.quantile)
    color.min <- ifelse(zero.as.min, 0, quantile(mat, min.quantile) )
    color.range <- c(color.min, color.max)
  }

  col_fun = colorRamp2(color.range, c("white", color))
  if (add.text) {
    cell_fun <- function(j, i, x, y, width, height, fill) {
      grid.text(round(mat[i, j], digits = digits), x, y,
                gp = gpar(fontsize = fontsize))
    }
  } else {
    cell_fun <- NULL
  }

  if (!is.null(row.assign.vec)) {
    row_split <- row.assign.vec[rownames(mat)] %>%
      factor(., levels = unique(.))
  } else {
    row_split <- NULL
  }
  if (!is.null(col.assign.vec)) {
    col_split <- col.assign.vec[colnames(mat)] %>%
      factor(., levels = unique(.))
  } else {
    col_split <- NULL
  }
  hm <- Heatmap(mat, col = col_fun,
                cell_fun = cell_fun,
                cluster_rows = cluster_rows, cluster_columns = cluster_columns,
                clustering_distance_rows = clustering_distance_rows,
                clustering_distance_columns = clustering_distance_columns,
                row_split = row_split, column_split = col_split, border = T
  )
  if (!is.null(plot.out)) {
    scFanc::save.base.plot(p = hm, file = plot.out, width = width, height = height, res = res)
  }
  return(hm)
}

contact.match.bg.core <- function(fg.gi, bg.gi, fg.id.col = NULL, bg.id.col = NULL,
                                  n.bg.each = 10, bezle = 0,
                                  normo.se, distance.wobble,
                                  out.dir) {
  dir.create(out.dir, showWarnings = F, recursive = T)
  if (bezle > 0)
    bg.gi <- gi.remove.bezles(gi = bg.gi, n = bezle)
  bg.gi <- bg.gi %>% subsetByOverlaps(fg.gi, invert = T)
  if (is.null(fg.id.col)) {
    fg.id.col <- "id"
    fg.gi$id <- paste0("fg_", 1:length(fg.gi))
  }
  if (is.null(bg.id.col)) {
    bg.id.col <- "id"
    bg.gi$id <- paste0("bg_", 1:length(bg.gi))
  }

  if (! fg.id.col %in% colnames(mcols(fg.gi)))
    stop("! fg.id.col %in% colnames(mcols(fg.gi))")
  if (! bg.id.col %in% colnames(mcols(bg.gi)))
    stop("! bg.id.col %in% colnames(mcols(bg.gi))")
  if (any(duplicated(mcols(fg.gi)[, fg.id.col]))) {
    stop("some ids are duplicated in fg")
  }
  if (any(duplicated(mcols(bg.gi)[, bg.id.col]))) {
    stop("some ids are duplicated in bg")
  }
  gr.cov <- rowRanges(normo.se)
  gr.cov$cov <- assay(normo.se, "VC..cpm") %>% rowMeans()

  regions(fg.gi) <- regions(fg.gi) %>% plyranges::join_overlap_inner(gr.cov)
  regions(bg.gi) <- regions(bg.gi) %>% plyranges::join_overlap_inner(gr.cov)
  bg.gi.chosen.list <- list()
  for (i in 1:length(fg.gi)) {
    fg <- fg.gi[i]
    fg.id <- mcols(fg)[, fg.id.col]
    dist <- pairdist(fg)
    bg <- bg.gi %>% .[pairdist(.) >= dist - distance.wobble & pairdist(.) <= dist + distance.wobble]
    gi.list <- list(fg = fg, bg = bg)
    cov.mat <- lapply(names(gi.list), function(gi.name) {
      gi <- gi.list[[gi.name]]
      df <- data.frame(a1 = anchors(gi, "first")$cov, a2 = anchors(gi, "second")$cov)
      if (gi.name == "fg")
        rownames(df) <- mcols(gi)[, fg.id.col]
      else
        rownames(df) <- mcols(gi)[, bg.id.col]
      mat <- df %>% as.matrix()
      return(mat)
    }) %>% Reduce(rbind, .)
    bg.id.chosen <- scFanc::bg.gen.2(mat = cov.mat, fg.vec = fg.id, n.bg.each = n.bg.each, no.replace = T,
                          method = "euclidean")

    # assess how well the bg is chosen:
    fgbg <- c(fg, bg)
    df <- data.frame(a1 = anchors(fgbg, "first")$cov,
                     a2 = anchors(fgbg, "second")$cov,
                     fg.id = mcols(fgbg)[, fg.id.col],
                     bg.id = mcols(fgbg)[, bg.id.col])
    df$hl <- NA
    df$hl[1:length(fg)] <- "fg"
    df$hl[df$bg.id %in% bg.id.chosen] <- "bg"
    pl <- list()
    pl$fg <- scFanc::xy.plot(df = df, x = "a1", y = "a2",
                             highlight.var = "fg.id", highlight.ptsize = 0.2,
                             highlight.values = fg.id, title = "fg")

    pl$bg <- scFanc::xy.plot(df = df, x = "a1", y = "a2",
                             highlight.var = "bg.id", highlight.ptsize = 0.2,
                             highlight.values = bg.id.chosen, title = "bg")
    trash <- scFanc::wrap.plots.fanc(plot.list = pl,
                                     plot.out = paste0(out.dir, "/assess/", fg.id, ".png"))
    # bg.list[[fg.id]] <- data.frame(fg.id = fg.id, bg.id = bg.id.chosen)

    bg.gi.chosen.list[[fg.id]] <- bg.gi %>% .[mcols(.)[, bg.id.col] %in% bg.id.chosen]
    bg.gi.chosen.list[[fg.id]]$fg.id <- fg.id

    bg.gi <- bg.gi %>% .[! mcols(.)[, bg.id.col] %in% bg.id.chosen]
  }
  bg.gi.chosen <- Reduce(c, bg.gi.chosen.list)
  x <- bg.gi.chosen
  x$value <- 1
  utilsFanc::gint.2.lr(gint = x, value.col = "value",
                       out.file = paste0(out.dir, "/bg.lr"))
  saveRDS(bg.gi.chosen, paste0(out.dir, "/bg.chosen.Rds"))
  return(bg.gi.chosen)

}

contact.map.collapse <- function(cs, regions.include = NULL, samples.include = NULL,
                                 out.dir,
                                 plot.each = F, per.sample = T, per.region = T,
                                 sample.coldata = NULL, split.column, split.vec = NULL,
                                 normalize = T, take.mean = T,
                                 center.piece = c(2,2),
                                 return.mat.list = F, ...) {
  # cs: contact summary result
  regions <- names(cs$details)
  if (!is.null(regions.include)) {
    utilsFanc::check.intersect(x = regions.include, y = regions)
    regions <- regions.include
  }
  # assuming all regions have the same samples

  mat.list <- lapply(regions, function(region) {
    region.sum <- cs$details[[region]]
    samples <- names(region.sum$expanded)
    if (!is.null(samples.include)) {
      utilsFanc::check.intersect(x = samples.include, y = samples)
      samples <- samples.include
    }
    pf <- parent.frame(n = 2)
    pf$samples <- samples
    mat.list <- lapply(samples, function(sample) {
      mat <- region.sum$expanded[[sample]]@matrix
      if (plot.each == T) {
        plot.out <- paste0(out.dir, "/per_mat/", region, "..", sample, ".png")
        trash <- contact.map.viz(mat = mat, plot.out = plot.out, ...)
      }
      return(mat)
    })
    names(mat.list) <- paste0(region, "..", samples)
    return(mat.list)
  }) %>% Reduce(c, .)

  if (normalize == T) {
    mat.list <- mat.list %>% lapply(function(x) return(x/sum(x)))
  }
  ###  #######
  center.df <- data.frame(
    center.pct = sapply(mat.list, function(x) return(x[center.piece[1],center.piece[2]])),
    region = sub("\\.\\..+$", "", names(mat.list)),
    sample = sub("^.+\\.\\.", "", names(mat.list))
  )
  rownames(center.df) <- NULL


  if (return.mat.list == T) {
    return(mat.list)
  }
  mat.sums <- list(mat.list = mat.list, center.df = center.df)
  if (per.sample == T) {
    mat.sums$per.sample <- lapply(samples, function(sample) {
      mat.list.sub <- mat.list[grepl(paste0("\\.\\.", sample),names(mat.list))]
      mat.sum <- mat.list.collapse.core(mat.list = mat.list.sub, take.mean = T,
                                        plot.out = paste0(out.dir, "/per_sample/", sample, ".png"), ...)
      # saveRDS(mat.sum, paste0(out.dir, "/per_sample/", sample, ".Rds"))
      return(mat.sum)
    })
    names(mat.sums$per.sample) <- samples
    p <- ggbarplot(data = center.df, x = "sample", y = "center.pct", add = c("mean_se", "jitter")) +
      coord_flip() +
      ggsave(paste0(out.dir, "/per_sample/summary_bar.pdf"), width = 7, height = 5, units = "in")

  }

  if (per.region == T) {
    mat.sums$per.region <- lapply(regions, function(region) {
      mat.list.sub <- mat.list[grepl(paste0(region, "\\.\\."),names(mat.list))]
      mat.sum <- mat.list.collapse.core(mat.list = mat.list.sub, take.mean = T,
                                        plot.out = paste0(out.dir, "/per_region/", region, ".png"), ...)
      # saveRDS(mat.sum, paste0(out.dir, "/per_region/", region, ".Rds"))
      return(mat.sum)
    })
    names(mat.sums$per.region) <- regions
    p <- ggbarplot(data = center.df, x = "region", y = "center.pct", add = c("mean_se", "jitter")) +
      coord_flip() +
      ggsave(paste0(out.dir, "/per_region/summary_bar.pdf"), width = 7, height = 0.3*length(regions) + 1, units = "in")
  }

  mat.sums$all <- mat.list.collapse.core(mat.list = mat.list, take.mean = T,
                                    plot.out = paste0(out.dir, "/collapse_all.png"), ...)
  if (!is.null(sample.coldata)) {
    utilsFanc::check.intersect(x = samples, y = rownames(sample.coldata),
                               x.name =  "samples", y.name = "rownames(sample.coldata")
    split.vec <- sample.coldata %>% as.data.frame() %>% .[samples, split.column]
    names(split.vec) <- NULL
  }
  if (!is.null(split.vec)) {
    if (!is.null(names(split.vec))) {
      utilsFanc::check.intersect(x = samples, y = names(split.vec),
                                 x.name =  "samples", y.name = "names(split.vec)")
      split.vec <- split.vec[samples]
    }
    mat.sums$all.split <- mat.list %>% split(., f = split.vec) %>%
      mapply(function(mat.l, name) {
        mat.list.collapse.core(mat.list = mat.l, take.mean = T,
                               plot.out = paste0(out.dir, "/collapse_",name,".png"), ...) %>%
          return()
      }, ., names(.), SIMPLIFY = F)
  }
  return(mat.sums)
}

mat.list.collapse.core <- function(mat.list, take.mean = T, keep.dimnames = F,
                                   plot.out = NULL, ...) {
  mat.sum <- Reduce(`+`, mat.list)
  if (take.mean == T)
    mat.sum <- mat.sum/length(mat.list)
  if (keep.dimnames == F) {
    dimnames(mat.sum) <- list(NULL,NULL)
  }
  if (!is.null(plot.out)) {
    dir.create(dirname(plot.out), showWarnings = F, recursive = T)
    trash <- contact.map.viz(mat = mat.sum, plot.out = plot.out, ...)
  }
  return(mat.sum)
}


SIPG.pipe <- function(fg.gi, hico.list, normo.list,
                      fg.ext,fg.id.col,
                      n.bg.each = 10, distance.wobble,
                      resolutions, mapqs, norms,
                      center.piece = NULL, color.range = NULL,
                      sample.coldata = NULL, split.column = "Ly49D",
                      out.dir, plot.only = F, plot.bar = T,
                      plot.mat.grid = T, fg.region.name.col = "forth",
                      replot.collapse = F,
                      stat.formula = NULL, stat.fun, stat.params = NULL,
                      debug = F, ...) {

  if (!is.null(sample.coldata)) {
    sample.coldata$sample <- rownames(sample.coldata)
    sample.coldata <- as.data.frame(sample.coldata)
  }
  results <- utilsFanc::safelapply(mapqs, function(mapq) {
    utilsFanc::safelapply(resolutions, function(res) {
      if (is.null(center.piece)) {
        center.piece <- 1+ fg.ext/res
        if (grepl("\\.", as.character(center.piece))) {
          paste0("center piece (", center.piece, ") contains decimals!")
        }
        center.piece <- rep(center.piece, 2)
      }
      utilsFanc::safelapply(norms, function(norm) {
        sub.dir <- paste0(out.dir, "/", mapq, "..res_", res, "..", norm, "/")
        if (plot.only == F) {
          hico <- hico.list[[mapq]][[paste0("res_", res)]]
          normo.se <- normo.list[[mapq]][[paste0("res_", res)]]
          bg.gi <- contact.match.bg.core(fg.gi = fg.gi, bg.gi = hico@interactions, fg.id.col = fg.id.col,
                                         bezle = ceiling(fg.ext/res),
                                         n.bg.each = n.bg.each, normo.se = normo.se,
                                         distance.wobble = distance.wobble,
                                         out.dir = paste0(sub.dir, "/bg/match_bg/"))
          fg.gi.ext <- fg.gi
          regions(fg.gi.ext) <- regions(fg.gi.ext) + fg.ext
          bg.gi.ext <- bg.gi
          regions(bg.gi.ext) <- regions(bg.gi.ext) + fg.ext

          gi.list <- list(fg = fg.gi.ext, bg = bg.gi.ext)
          pipe.list <- utilsFanc::safelapply(names(gi.list), function(gi.name) {
            gi <- gi.list[[gi.name]]
            if (gi.name == "fg") {
              gi.name.col <- fg.id.col
              plot.each <- F # turned off because I will plot mat grid instead
              per.region <- F
            } else {
              gi.name.col <- "id"
              plot.each <- F
              per.region <- F
            }
            pipe.out <- contact.pipe(gi = gi, gi.name.col = gi.name.col, intSet = hico,
                                     assay = norm, sum.fun = sum, resolution = res,
                                     sum.to.1 = T, write.anchors = T, do.donut = F,
                                     out.dir = paste0(sub.dir, "/", gi.name,"/contact_pipe/"),
                                     do.plot = F)
            mat.sums <- contact.map.collapse(cs = pipe.out$contact.stat,
                                             out.dir = paste0(sub.dir, "/", gi.name,"/collapse/"),
                                             plot.each = plot.each, per.sample = T, per.region = per.region,
                                             normalize = T, take.mean = T, center.piece = center.piece,
                                             color.range = color.range, width = 550, height = 500,
                                             add.text = F, sample.coldata = sample.coldata,split.column = split.column,
                                             ...)
            res <- list(pipe.out = pipe.out, mat.sums = mat.sums)
          }, threads = 1)
          names(pipe.list) <- names(gi.list)
          center.df <- lapply(names(pipe.list), function(name) {
            res <- pipe.list[[name]]$mat.sums$center.df
            res$type <- name
            return(res)
          }) %>% Reduce(rbind, .)
          center.df.2 <- center.df %>% group_by(region, type) %>% summarise(center.pct = mean(center.pct)) %>%
            ungroup() %>% as.data.frame()
          if (!is.null(sample.coldata)) {
            center.df.3 <- left_join(center.df, sample.coldata[, c("sample", split.column)]) %>%
              group_by(region, type, !!as.name(split.column)) %>% summarise(center.pct = mean(center.pct)) %>%
              ungroup() %>% as.data.frame()
          }
          else {
            center.df.3 <- NULL
          }
          result <- list(fg.gi = fg.gi, fg.gi.ext = fg.gi.ext,
                         bg.gi = bg.gi, bg.gi.ext = bg.gi.ext,
                         pipe.list = pipe.list, center.df = center.df, center.df.2 = center.df.2,
                         center.df.3 = center.df.3)
          pl <- list()
        } else {
          result <- readRDS(paste0(sub.dir, "/result.Rds"))
          pl <- result$pl
        }
        if (plot.bar == T) {
          pl$per.sample <- ggbarplot(data = result$center.df, x = "sample", y = "center.pct",
                                     color = "type", add = c("mean_sd"),
                                     position = position_dodge()) +
            geom_jitter(size = 0.05, position = position_jitterdodge(jitter.width = 0.2), aes(color = type),
                        alpha = 0.5) +
            coord_flip() +
            ggsave(paste0(sub.dir, "/contrast_per_sample.pdf"), width = 6, height = 5, units = "in", dpi = 100)
          pl$mean <- ggbarplot(data = result$center.df.2, x = "type", y = "center.pct",
                               color = "black", fill = "type", add = c("mean_sd"),
                               add.params = list(size = 1.0), alpha = 0.5) +
            geom_jitter(size = 0.5, aes(fill = type), color = "grey50",
                        alpha = 1) +
            coord_flip() + theme(aspect.ratio = 1) + theme_pubr(base_size = 18) +
            ggsave(paste0(sub.dir, "/contrast_collapse.pdf"), width = 6, height = 5, units = "in", dpi = 100)
          if (!is.null(stat.formula)) {
            p.df <- lapply(stat.formula, function(formula) {
              vec.list <- formula.parse.fanc(formula = formula, df = result$center.df.2)
              p <- stats.wrapper.ez(vec.list = vec.list, stat.fun = stat.fun, stat.params = stat.params, slot = "p.value")
              if (p < 0.0001) {
                p <- formatC(p, format = "e", digits = 2)
              }
              res <- data.frame(comp = formula, p.value = p)
              return(res)
            }) %>% Reduce(rbind, .)
            write.table(p.df, paste0(sub.dir, "/contrast_collapse_pvalue.tsv"), sep = "\t", quote = F,
                        row.names = F, col.names = T)
          }

          if (!is.null(sample.coldata)) {
            if (is.null(result$center.df.3)) {
              result$center.df.3 <- left_join(result$center.df, sample.coldata[, c("sample", split.column)]) %>%
                group_by(region, type, !!as.name(split.column)) %>% summarise(center.pct = mean(center.pct)) %>%
                ungroup() %>% as.data.frame()
            }
            pl$mean_split <- ggbarplot(data = result$center.df.3, x = split.column, y = "center.pct",
                                 color = "black", fill = "type",  alpha = 0.5, size = 1,
                                 add = c("mean_sd"), position = position_dodge(),
                                 add.params = list(size = 1.0)) +
              geom_jitter(size = 0.8, aes(fill = type), color = "grey50",
                          alpha = 1, position = position_jitterdodge()) +
              coord_flip() +
              theme_pubr(base_size = 18) + theme(aspect.ratio = 1) +
              ggsave(paste0(sub.dir, "/contrast_collapse_split.pdf"), width = 6, height = 5, units = "in", dpi = 100)
          }
        }
        if (plot.mat.grid == T) {
          fg.gr <- gi.get.active.regions(fg.gi)
          region.order <- sort(fg.gr) %>% mcols() %>% .[, fg.region.name.col]
          fg.bin.assign <- gr.slide.2(fg.gr, name.col = fg.region.name.col,
                                      bin.size = res, n.bins = ceiling(fg.ext/res), return.assign.vec = T)
          names(fg.bin.assign) <- paste0(fg.bin.assign, "..", names(fg.bin.assign))
          if (1) {
            # is.null(result$pipe.list$fg$mat.sums$per.region)
            result$pipe.list$fg$mat.sums$per.region <- result$pipe.list$fg$mat.sums$mat.list %>%
              split(., f = factor(sub("\\.\\..+$", "", names(.)),
                                  levels = unique(sub("\\.\\..+$", "", names(.))))) %>%
              lapply(function(mats) {
                in.mat.list <- list(all = mats)
                if (!is.null(sample.coldata)) {
                  samples <- names(mats) %>% sub(".+\\.\\.", "", .)
                  split.id <- sample.coldata %>% as.data.frame() %>% .[samples, split.column]
                  in.mat.list <- c(in.mat.list,
                               mats %>% split(., f = split.id))
                }
                out.mat.list <- lapply(in.mat.list, function(mats) {
                  mat.list.collapse.core(mat.list = mats, take.mean = T, keep.dimnames = T)
                })

                region <- names(mats)[1] %>% sub("\\.\\..+$", "", .)
                names(out.mat.list) <- paste0(region, "..", "mean_", names(out.mat.list))
                return(out.mat.list)
              }) %>% Reduce(c, .)
          }


          mat.list.to.grid <- c(result$pipe.list$fg$mat.sums$mat.list,
                                result$pipe.list$fg$mat.sums$per.region)
          pl$mat.grid <-  mat.list.to.grid %>% split(., f = sub("^.+\\.\\.", "", names(.))) %>%
            lapply(function(mat.l) {
              sample <- names(mat.l)[1] %>% sub("^.+\\.\\.", "", .)
              plot.out <- paste0(sub.dir, "/mat_grid/", sample, ".png")
              mat.grid <- mat.l %>%
                mat.add.by.dimnames(parse.region.ids = T, sep = "_",
                                    col.region.order = region.order,
                                    row.region.order = region.order, add.transpose = T)

              p <- contact.map.viz(mat = mat.grid, row.assign.vec = fg.bin.assign,
                                   col.assign.vec = fg.bin.assign, add.text = F,
                                   width = 1600, height = 1500,
                                   plot.out = plot.out, ...)
              return(p)
            })

        }
        if (replot.collapse == T) {
          lapply(c("fg", "bg"), function(gi.name) {
            # per.region <- ifelse(gi.name == "fg", T, F)
            # plot.each <- ifelse(gi.name == "fg", T, F)
            # per.sample <- T
            per.region <- F
            plot.each <- F
            per.sample <- F
            trash <- contact.map.collapse(cs = result$pipe.list[[gi.name]]$pipe.out$contact.stat,
                                          out.dir = paste0(sub.dir, "/", gi.name,"/collapse/"),
                                          plot.each = plot.each, per.sample = per.sample, per.region = per.region,
                                          normalize = T, take.mean = T, center.piece = center.piece,
                                          color.range = color.range, width = 550, height = 500,
                                          add.text = F, sample.coldata = sample.coldata, split.column = "Ly49D",
                                          ...)
            return()
          })

        }
        result$pl <- pl
        saveRDS(result, paste0(sub.dir, "/result.Rds"))
        return(paste0(sub.dir, "/result.Rds"))
      }, threads = ifelse(debug == T, 1, length(norms))) %>% Reduce(c, .) %>% return()
    }, threads = ifelse(debug == T, 1, length(resolutions))) %>% Reduce(c, .) %>% return()
  }, threads = ifelse(debug == T, 1, length(mapqs))) %>% Reduce(c, .) %>% return()
  return(results)
}



t.f.parent.frame <- function(x) {
  lapply(1:x, function(i) {
    p <- parent.frame(2)
    p$miao.parent.frame <- i
    return()
  })
  print(miao.parent.frame)
}

gi.remove.bezles <- function(gi, n) {
  if (n < 1) {
    return(gi)
  }
  all.used.regions <- anchorIds(gi) %>% unlist() %>% unique() %>% sort()
  bezle.left <- all.used.regions[1:n]
  bezle.right <- all.used.regions %>% .[(length(.) - n + 1):length(.)]
  bezles <- c(bezle.left, bezle.right)
  gi <- gi[ ! anchorIds(gi, "first") %in% bezles &
            ! anchorIds(gi, "second") %in% bezles ]
  return(gi)
}

hico.list.gen <- function(inter.hic, inter30.hic, region, colData, save.rds = NULL) {
  hico.list <- lapply(list(inter.hic, inter30.hic), function(.hic.vec) {
    print(paste0("Processing these files: \n",
                 paste0(.hic.vec, collapse = "\n")))
    hico <- read.hic(.hic.vec = .hic.vec, region1 = region, observed.or.oe = "observed",
                     normalization = c("NONE", "VC_SQRT", "VC"), use.1.6 = T,
                     threads.norm = NULL, threads.res = 1, threads.sample = NULL,
                     colData = colData)
    return(hico)
  })
  names(hico.list) <- c("inter.hic", "inter_30.hic")
  if (!is.null(save.rds)) {
    dir.create(dirname(save.rds), showWarnings = F, recursive = T)
    saveRDS(hico.list, save.rds)
  }
  return(hico.list)
}

