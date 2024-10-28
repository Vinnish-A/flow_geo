
# download --------------------------------------------------------------

library(GEOquery)
library(rlang)

inTmpDir = function(expr_, dest_) {
  
  path_raw_ = getwd()
  
  if (!dir.exists(dest_)) dir.create(dest_, recursive = T)
  setwd(dest_)
  
  eval_tidy(expr_)
  
  setwd(path_raw_)
  
}

gzDir = function(dir_) {
  
  filenames_ = list.files(dir_, full.names = T, pattern = 'gz')
  if (length(filenames_) == 0) {
    cat('NULL')
    invisible(NULL)
  } else {
    lapply(filenames_, R.utils::gunzip, remove = T)
    invisible(NULL)
  }
  
}

