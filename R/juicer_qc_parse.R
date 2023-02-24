juicer.qc.parse <- function(inter.txt.vec, out.file=NULL,
                            field.file = "~/R_packages/v4c/juicer.qc.fields.tsv", on.target.tsv = NULL) {
  if (is.null(names(inter.txt.vec)))
    stop("inter.txt.vec must be named")
  qc.df <- lapply(seq_along(inter.txt.vec), function(i) {
    table <- tempfile()
    cmd <- paste0("sed 's/general.*/general/' ", inter.txt.vec[i], " | sed 's/: */\t/g' > ", table)
    print(cmd); system(cmd)
    qc.vec <- read.table(table, header = F, sep = "\t", quote = "")$V2
    if (length(qc.vec) == 19) {
      qc.vec <- append(qc.vec, "mega", after = 0)
      qc.vec <- append(qc.vec, "NA", after = 6)
      qc.vec <- append(qc.vec, "NA", after = 11)
    }
    return(qc.vec) 
  }) %>% Reduce(cbind, .) %>% as.data.frame()
  colnames(qc.df) <- names(inter.txt.vec)
  field.df <- read.table(field.file, sep = "\t", quote = "")
  colnames(field.df) <- "qc"
  
  qc.df <- utilsFanc::add.column.fanc(df1 = qc.df, df2 = field.df, pos = 1)
  
  if (!is.null(on.target.tsv)) {
    # if (is.null(names(on.target.tsv)))
    #   stop("on.target.tsv must be a named list")
    
    samples <- names(inter.txt.vec)
    on.target.dfs <- lapply(seq_along(on.target.tsv), function(i) {
      tsv <- read.table(on.target.tsv[[i]], header = T)
      field <- names(on.target.tsv[i])
      tsv <- tsv[, c("sample.name", "space","start", "end", "frac.on.target")] %>% 
        mutate(pos = paste0(space , ":", start, "-", end)) %>% 
        mutate(space = NULL, start = NULL, end = NULL)
      df <- reshape2::acast(data = tsv, formula = pos~sample.name, value.var = "frac.on.target")
      df <- utilsFanc::add.column.fanc(df1 = df, df2 = data.frame(qc = rownames(df)), pos = 1)
      df <- df %>% mutate(qc = paste0(field, "::", qc))
      return(df)
    })
    qc.df <- Reduce(rbind, c(list(qc.df), on.target.dfs))
  }
  if (!is.null(out.file)) {
    write.table(qc.df, out.file, quote = F, col.names = T, row.names = F, sep = "\t")
    paste0("sed 's/([^()]*)//g' ", out.file, " > ", utilsFanc::insert.name.before.ext(name = out.file, insert = "ez", delim = "_"))
  }
    
  return(qc.df)
}


# t <- juicer.qc.parse(inter.txt.vec = inter.txt.vec, out.file = "qc.tsv",
#                      on.target.tsv = list(on.target.rate.raw = "capture_stats_raw.tsv",
#                                           on.target.rate.dedup = "capture_stats_dedup.tsv"))