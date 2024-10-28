
# GSE13507 --------------------------------------------------------------

source('utils/download.R')

library(tidyverse)
library(AnnoProbe)
library(fromto)
library(rlang)
library(glue)

GSE = 'GSE13507'

## Illumina human-6 v2.0 expression beadchip ----

# Raw data missing
# deprecated

library(lumi)

data_raw = lumiR('data/GSE13507/raw/GSE13507_illumina_raw.txt', lib.mapping = 'Human-6-v2')
data_norm = lumiB(raw_data, method = 'bgAdjust')
data_norm = lumiN(norm_data, method = 'quantile')

## Series Mat ----

source('utils/series.R')

inTmpDir(geoget(GEO_IDs = GSE), glue('data/{ GSE }/series/'))
gzDir(glue('data/{ GSE }/series/'))

series_all = readSeries(glue('data/{ GSE }/series/'))
series_all$phen = series_all$phen |> 
  rename(OS = `overall survival`, OS_time = `survival month`) |> 
  mutate(OS = ifelse(str_detect(OS, 'survial'), 0, 1), 
         OS_time = as.numeric(OS_time), 
         OS_unit = 'month')

writeSeries(series_all, glue('data/{ GSE }/tidy/'))

## diff ----

source('utils/diff.R')

library(furrr)
library(transGI)

load(glue('data/{ GSE }/tidy/all.RData'))

multisession = function(expr_, n_ = 4) {
  
  if (future::nbrOfWorkers() > 1) plan('sequential')
  
  plan('multisession', workers = 6)
  tidy_eval(expr_)
  plan('sequential')
  
}

# uni-cox

lst_HR = future_imap(
  expr[, -1], 
  \(vec_, gene_) {
    
    res_ = list()
    
    res_[[gene_]] = phen |> 
      mutate(value = vec_) |> 
      select(value, OS_time, OS) |> 
      drop_na(value, OS_time, OS) |> 
      docal(cal_HR, value, OS_time, OS)
    
    return(res_)
    
  }
)
