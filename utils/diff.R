
# diff --------------------------------------------------------------


diff_limma = function(normal_matrix, group) {
  
  library(limma)
  
  design = model.matrix(~0 + group)
  colnames(design) = c('group1', 'group2')
  
  contrast = makeContrasts(group2 - group1, levels = design)
  
  fit = lmFit(normal_matrix, design)
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
