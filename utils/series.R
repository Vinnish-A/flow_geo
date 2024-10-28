
# series --------------------------------------------------------------

library(GEOquery)

readSeries = function(dir_) {
  
  series_ = list.files(dir_, full.names = T)
  gse_ = getGEO(filename = series_, getGPL = F)
  
  expr_ = gse_ |> Biobase::exprs() |> as_tibble(rownames = 'ID')
  phen_ = pData(gse_) |> as_tibble(rownames = 'sample')  
  
  if (checkGPL(gse_@annotation)) {
    
    table_id_ = idmap(gse_@annotation, type = 'soft', destdir = 'tmp') |> 
      pull(symbol, ID)
    
  } else {
    
    ids_ = unique(expr_$ID)
    table_id_ = setNames(ids_, ids_)
    
  }
  
  expr_tidy_ = expr_ |> 
    mutate(symbol = table_id_[ID]) |> 
    filter(symbol != 'permuted_negative') |> 
    distinct(symbol, .keep_all = T) |> 
    column_to_rownames('symbol') |> 
    select(-ID) |> 
    t() |> 
    as_tibble(rownames = 'sample')
  
  phen_tidy_ = expr_tidy_[, 1] |> 
    left_join(phen_) |> 
    select(sample, title, contains(':ch')) |> 
    rename_with(~ .x |> str_sub(end = -5), contains(':ch'))
  
  return(list(expr = expr_tidy_, phen = phen_tidy_))
  
}

writeSeries = function(series_, dir_, RData = T, gse_ = NULL) {
  
  if (!dir.exists(dir_)) dir.create(dir_, recursive = T)
  
  if (!is.null(gse_)) names(series_) = paste(names(series_), gse_, sep = '_')
  
  cat('Writing expr...\n')
  write_csv(series_$expr, normalizePath(file.path(dir_, 'expr.csv'), '/', F))
  cat('Writing phen...\n')
  write_csv(series_$phen, normalizePath(file.path(dir_, 'phen.csv'), '/', F))
  
  if (RData) {
    cat('Writing RData...\n')
    env_ = list2env(series_)
    save(list = ls(envir = env_), file = normalizePath(file.path(dir_, 'all.RData'), '/', F), envir = env_)
  }
  
}
