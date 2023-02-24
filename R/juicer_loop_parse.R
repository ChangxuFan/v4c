loop.extract.motif <- function(loop.df, out.file) {
  df.list <- list()
  df.list$df1 <- loop.df[, c("chr1", "motif_x1", "motif_x2", "uniqueness_1", "orientation_1")]
  df.list$df2 <- loop.df[, c("chr2", "motif_y1", "motif_y2", "uniqueness_2", "orientation_2")]
  df <- df.list %>% lapply(function(df) {
    df <- df %>% na.omit()
    colnames(df) <- c("chr", "start", "end", "forth", "strand")
    df <- utilsFanc::add.column.fanc(df, data.frame(fifth = rep(".", nrow(df))), pos = 5)
    df$strand[df$strand == "p"] <- "+"
    df$strand[df$strand == "n"] <- "-"
    return(df)
  }) %>% Reduce(rbind, .)
  utilsFanc::write.zip.fanc(df = df, out.file = out.file, bed.shift = F)
  return(df)
}

loop.extract.anchor <- function(loop.df, out.file) {
  bMotif <- any(grepl("motif", colnames(loop.df)))
  df.list <- list()
  left.columns <- c("chr1", "x1", "x2", "uniqueness_1", "orientation_1")
  right.columns <- c("chr2", "y1", "y2", "uniqueness_2", "orientation_2")
  if (!bMotif) {
    df$uniqueness_1 <- "NA"
    df$uniqueness_2 <- "NA"
    df$orientation_1 <- "+"
    df$orientation_2 <- "-"
  }
  df.list$df1 <- loop.df[, left.columns]
  df.list$df2 <- loop.df[, right.columns]
  
  df <- df.list %>% lapply(function(df) {
    # df <- df %>% na.omit()
    colnames(df) <- c("chr", "start", "end", "forth", "strand")
    df <- utilsFanc::add.column.fanc(df, data.frame(fifth = rep(".", nrow(df))), pos = 5)
    df$strand[df$strand == "p"] <- "+"
    df$strand[df$strand == "n"] <- "-"
    df$strand[is.na(df$strand)] <- "."
    df$forth[is.na(df$forth)] <- "NA"
    return(df)
  }) %>% Reduce(rbind, .)
  utilsFanc::write.zip.fanc(df = df, out.file = out.file, bed.shift = F)
  return(df)
}

