
# normalized --------------------------------------------------------------


library(tidyverse)
library(GEOquery)
library(AnnoProbe)

gse = getGEO(filename = "data/GSE143272_series_matrix.txt", getGPL = F)

expr = gse |> exprs() |> as_tibble(rownames = 'ID')
phen = pData(gse) |> as_tibble(rownames = 'sample')  

table_id = idmap(gse@annotation, type = 'soft', destdir = 'tmp') |> 
  pull(symbol, ID)

expr_tidy = expr |> 
  mutate(symbol = table_id[ID]) |> 
  filter(symbol != 'permuted_negative') |> 
  distinct(symbol, .keep_all = T) |> 
  column_to_rownames('symbol') |> 
  select(-ID) |> 
  t() |> 
  as_tibble(rownames = 'sample')

phen_tidy = expr_tidy[, 1] |> 
  left_join(phen) |> 
  mutate(group = ifelse(str_detect(title, 'Healthy'), 'Healthy', 'Diseased'))

input = expr_tidy |> 
  column_to_rownames('sample') |> 
  t()
group = phen_tidy$group |> factor(levels = c('Healthy', 'Diseased'))

## diff ----

diff_limma = function(log2_tpm, group) {
  
  library(limma)
  
  design = model.matrix(~0 + group)
  colnames(design) = c('group1', 'group2')
  
  contrast = makeContrasts(group2 - group1, levels = design)
  
  fit = lmFit(log2_tpm, design)
  fit = contrasts.fit(fit, contrast)
  fit = eBayes(fit)
  
  
  result = topTable(fit, adjust.method = "BH", number = Inf)
  
  return(as_tibble(result, rownames = 'symbol'))
  
}


diff_deseq2 = function(count_matrix, group) {
  
  library(DESeq2)
  
  colData = data.frame(group = factor(group))
  dds = DESeqDataSetFromMatrix(
    countData = count_matrix,
    colData = colData,
    design = ~ group
  )
  
  dds = DESeq(dds)
  res = results(dds)
  
  res_df = as.data.frame(res)
  res_df = res_df[order(res_df$padj), ]
  
  return(as_tibble(res_df, rownames = 'symbol'))
  
}

diff_edger = function(count_matrix, group) {
  
  library(edgeR)
  
  group = factor(group)
  y = DGEList(counts = count_matrix, group = group)
  
  keep = filterByExpr(y)
  y = y[keep, , keep.lib.sizes = FALSE]
  
  y = calcNormFactors(y)
  design = model.matrix(~ group)
  y = estimateDisp(y, design)
  
  fit = glmFit(y, design)
  lrt = glmLRT(fit, coef = 2)  
  res = topTags(lrt, n = Inf)
  
  res_df = as.data.frame(res)
  
  return(as_tibble(res_df, rownames = 'symbol'))
  
}

deg_limma = diff_limma(input, group)

deg_limma |> 
  filter(adj.P.Val < 0.05) |> 
  nrow()
