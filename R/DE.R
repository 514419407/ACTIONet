identify.DE.genes.wilcox <- function(ACTIONet.out, sce, Labels, gene.per.class = 10, pval.threshold = 1e-100, AUC.threshold = 0.5) {

  IDX = split(1:length(Labels), Labels)
  
  C = sapply(IDX, function(idx) as.numeric(Matrix::sparseVector(x = 1/length(idx), i = idx, length = length(Labels))))
  
  
	reduced.profile = (t(sce@reducedDims[["S_r"]]) %*% C)
	
	Annot = rowData(sce)
	cnames = colnames(Annot)
	idx = grep('^PC', cnames)
	V = as.matrix(Annot[, idx])
	
	perm = order(sapply(cnames[idx], function(str) as.numeric(stringr::str_replace(str, "PC", ""))))
	V = V[, perm]
	
	gene.profile = V %*% reduced.profile 
	rownames(gene.profile) = rownames(sce)
  
	sorted.genes = apply(gene.profile, 2, function(x) {
	 return(rownames(sce)[order(x, decreasing = TRUE)])
	})
	
	sce$groups = Labels
	wilc.out = presto::wilcoxauc(sce, "groups")
	
	# marker.genes = presto::top_markers(wilc.out, n = gene.per.class, auc_min = AUC.threshold, padj_max = pval.threshold)[-1]
	# 
	# wilc.out = wilc.out[(wilc.out$auc > AUC.threshold) & (wilc.out$padj < pval.threshold), ]

	ii = match(wilc.out$feature, rownames(sce))
  jj = match(wilc.out$group, levels(Labels))
  vv = -log10(wilc.out$padj)
  vv[wilc.out$auc < AUC.threshold] = 0
  logPval.mat = as.matrix(Matrix::sparseMatrix(i = ii, j = jj, x = vv, dims = c(nrow(sce), length(levels(Labels)))))
  rownames(logPval.mat) = rownames(sce)
  colnames(logPval.mat) = levels(Labels)
  
  counts = Matrix::rowSums(logPval.mat > -log10(pval.threshold))
  
  logPval.mat = logPval.mat[counts == 1, ]

  sig.genes = split(rownames(logPval.mat), apply(logPval.mat, 1, which.max))
  names(sig.genes) = colnames(logPval.mat)
  
	marker.genes = sapply(colnames(gene.profile), function(celltype) {
	  x = gene.profile[, celltype]
	  sorted.genes = rownames(sce)[order(x, decreasing = TRUE)]
	  
	  sorted.genes = sorted.genes[(sorted.genes %in% sig.genes[[celltype]])]

	  return(sorted.genes[1:min(gene.per.class, length(sorted.genes))])
	})
	
	return(marker.genes)
}


identify.gene.modules <- function(ACTIONet.out, sce, Labels, interactome = NA) {
  GCN = cor(t(ACTIONet.out$signature.profile))
  rownames(GCN) = colnames(GCN) = rownames(sce)
  diag(GCN) = 0
  
  if(!is.na(interactome)) {
    common.genes = intersect(rownames(sce), rownames(interactome))
    GCN = GCN[common.genes, common.genes]
    interactome = interactome[common.genes, common.genes]
    
    interactome = as(interactome, 'dgTMatrix')
    idx = (interactome@i)*length(common.genes) + (interactome@j + 1)
    w = GCN[idx]
    
    Adj = Matrix::sparseMatrix(i = interactome@i + 1, j = interactome@j + 1, w, dims = dim(interactome))
    rownames(Adj) = colnames(Adj) = common.genes
  } else {
    Adj = as(GCN, 'sparseMatrix')
  }
  
  gene.clusters = ACTIONet::signed_cluster(Adj, seed = 0)
  
  IDX = split(1:length(gene.clusters), gene.clusters)
  IDX = IDX[sapply(IDX, length) >= 30]
  
  module.genes = sapply(IDX, function(idx) rownames(Adj)[idx])
  
  
  IDX = split(1:length(Labels), Labels)
  
  C = sapply(IDX, function(idx) as.numeric(Matrix::sparseVector(x = 1/length(idx), i = idx, length = length(Labels))))
  
  
	reduced.profile = (t(sce@reducedDims[["S_r"]]) %*% C)
	
	Annot = rowData(sce)
	cnames = colnames(Annot)
	idx = grep('^PC', cnames)
	V = as.matrix(Annot[, idx])
	
	perm = order(sapply(cnames[idx], function(str) as.numeric(stringr::str_replace(str, "PC", ""))))
	V = V[, perm]
	
	gene.profile = V %*% reduced.profile 
	rownames(gene.profile) = rownames(sce)

	
  spec.scores = t(scale(sapply(module.genes, function(genes) {
    scores = Matrix::colMeans(gene.profile[genes, ])
  })))

  	
	return(list(Modules = module.genes, Specificity = spec.scores))
}



identify.phenotype.associated.genes.continuous <- function(A, sce, rand.no = 100, seed = 0) {
	
	Sr = t(sce@reducedDims[["S_r"]])

	Annot = rowData(sce)
	cnames = colnames(Annot)
	idx = grep('^PC', cnames)
	V = as.matrix(Annot[, idx])
	perm = order(sapply(cnames[idx], function(str) as.numeric(stringr::str_replace(str, "PC", ""))))
	V = V[, perm]


	A = scale(A)  

	N = nrow(A)
	O = as.matrix(array(1, N))

	SrSrt = Sr %*% t(Sr)
	SrO = Sr %*% O
	Core = SrSrt - (SrO %*% t(SrO)) / N
	coVar = (V %*% Core %*% t(V)) / N
	gene.sd = sqrt(diag(coVar))


	Stats = V %*% (Sr %*% A)
	Stats.scaled = apply(Stats, 2, function(v) return(v / gene.sd)) / N

	
	set.seed(seed)
  Enrichment.Z = sapply(1:ncol(Stats), function(k) {
		print(k)

		a = A[, k]

		R = replicate(rand.no, sample(a))

		randStats = as.matrix(V %*% (Sr %*% R))
		randStats.scaled = apply(randStats, 2, function(v) return(v / gene.sd)) / N


		Mu = apply(randStats.scaled, 1, mean)
		Std = (apply(randStats.scaled, 1, sd))

		z = as.matrix((Stats.scaled[, k] - Mu) / Std)

		return(z)
  })
  rownames(Enrichment.Z) = rownames(sce)
  # colnames(Enrichment.Z) = arch.Labels.human[ACTIONet.out$core.out$core.archs]
	
  return(Stats)
}

identify.phenotype.associated.genes.discrete <- function(Labels, sce, rand.no = 100, seed = 0) {
	if(!is.factor(Labels)) {
		Labels = factor(as.character(Labels), as.character(sort(unique(Labels))))
	}

	A = sapply(levels(Labels), function(L) as.numeric(Labels == L))

	Enrichment.Z = identify.phenotype.associated.genes.continuous(A = A, sce = sce, rand.no = rand.no, seed = seed)


	return(Enrichment.Z)
}

identify.ACTIONet.associated.genes <- function(ACTIONet.out, sce) {

  A = as(ACTIONet.out$build.out$ACTIONet, 'dgTMatrix')
  
  # Network enhancement
  rs = Matrix::rowSums(A)
  rs[rs == 0] = 1
  # P = t(scale(t(A), center = FALSE, scale = rs))
  P = as(sparseMatrix(i = A@i + 1, j = A@j + 1, x = A@x / rs[A@i + 1], dims = dim(A)), 'dgTMatrix')

  eps = 1e-16
  w = sqrt(Matrix::colSums(P)+eps);
  W = sparseMatrix(i = P@i + 1, j = P@j + 1, x = P@x / w[A@j + 1], dims = dim(A))
  P.NE = W %*% Matrix::t(W)
  
  X = sce@assays[["logcounts"]]
  

	Sr_t = sce@reducedDims[["S_r"]]

	Annot = rowData(sce)
	cnames = colnames(Annot)
	idx = grep('^PC', cnames)
	V = as.matrix(Annot[, idx])
	perm = order(sapply(cnames[idx], function(str) as.numeric(stringr::str_replace(str, "PC", ""))))
	V = V[, perm]


	kxk = (Matrix::t(Sr_t) %*% P.NE) %*% Sr_t
	gxg = (V %*% kxk) %*% Matrix::t(V) # x'Lx
  	
	gene.scores = diag(gxg)

	names(gene.scores) = rownames(sce)

	return(gene.scores)
}


# identify.ACTIONet.associated.genes <- function(ACTIONet.out, sce) {
# 	Sr_t = sce@reducedDims[["S_r"]]
# 	
# 	Annot = rowData(sce)
# 	cnames = colnames(Annot)
# 	idx = grep('^PC', cnames)
# 	V = as.matrix(Annot[, idx])
# 	perm = order(sapply(cnames[idx], function(str) as.numeric(stringr::str_replace(str, "PC", ""))))
# 	V = V[, perm]
# 
# 	
# 	kxk = (t(Sr_t) %*% L) %*% Sr_t
# 	gxg = (V %*% kxk) %*% t(V) # x'Lx
# 	
# 	num = as.numeric(diag(gxg))
# 
# 	# Normalize
# 	N = nrow(Sr_t)
# 	Wtot = sum(ACTIONet.out$build.out$ACTIONet)
# 	
# 	X = (t(Sr_t) %*% Sr_t) %*% t(V)
# 	
# 	o = as.matrix(array(1, N))
# 	s = t(o) %*% Sr_t
# 	Y = t(s) %*% (s  %*% t(V)) / N
# 	
# 	Delta = X - Y
# 	
# 	denom = as.numeric(diag(V %*% (X - Y)))
# 
# 	gene.autocorrelation.C = ((N-1) /  Wtot) * (num / denom)
# 
# # ## Test::
# # 	i = 110
# # 	x = V[i, ] %*% t(Sr_t)
# # 
# # 	
# # 	mu = mean(x)
# #   x_c = (x-mu)
# # 
# #   x%*%L%*%t(x)
# #   num[i]
# #   
# #   sum( x_c^2 )
# #   denom[i]
# # 
# #   ((N-1) / (Wtot)) * (num[i] / sum( (x_c)^2 ))
# # 	gene.autocorrelation.C[i]
# 	
# 	perm0  = order(num, decreasing = TRUE)
#   rownames(sce)[perm0[1:10]]
#   num[perm0[1:10]]
#   denom[perm0[1:10]]
#   
#   Matrix::rowMeans(gxg[perm0[1:10], ])
#   
# 
# 	perm = order(gene.autocorrelation.C, decreasing = TRUE)
#   rownames(sce)[perm[1:10]]
#   num[perm[1:10]]
#   denom[perm[1:10]]
#   
# 
# 	names(gene.autocorrelation.C) = rownames(sce)
# 
# 	gene.autocorrelation.C = as.matrix(gene.autocorrelation.C)
# 	
#   return(gene.autocorrelation.C)
# }


visualize.expression.heatmap <- function(ACTIONet.out, sce, genes, Labels, scale.rows = TRUE, filter.duplicates = FALSE, reorder.genes = FALSE, prune.expression = FALSE, CPal = "d3") {
require(ComplexHeatmap)

  if(!is.factor(Labels)) {
    Labels = factor(as.character(Labels), levels = sort(unique(as.character(Labels))))
  }
    
  imputed.genes = impute.genes.using.ACTIONet(ACTIONet.out, sce, genes, prune = prune.expression, rescale = TRUE)

  if(filter.duplicates) {
   S = t(imputed.genes[order(Labels), order(match(colnames(imputed.genes), genes))]) 
   rownames(S) = colnames(imputed.genes)
  } else {
    S = t(imputed.genes[order(Labels), genes])
    rownames(S) = genes
  }
  colnames(S) = sort(Labels)

  # Filter outliers
  Z = (S - median(S))/mad(S)


  
  # Normalize each genes across samples
  if(scale.rows == TRUE) {
  #   
  #   # Dr = Matrix::Diagonal(nrow(S), 1 / sqrt(as.numeric(Matrix::rowSums(S))))
  #   # Dc = Matrix::Diagonal(ncol(S), 1 / sqrt(as.numeric(Matrix::colSums(S))))
  #   # Z = as.matrix(Dr %*% S %*% Dc)
  # 
    ss = apply(Z, 1, sd)
    ss[ss == 0] = 1
    Z = t(scale(t(Z), center = TRUE, scale = ss))

   z.threshold = 3.5
   Z[Z > z.threshold] = z.threshold
   Z[Z < -z.threshold] = -z.threshold
    
  }
  # else {
  #   Z = S
  # }

  
  if(reorder.genes == TRUE) {
    set.seed(0)
    
    IDX = split(1:ncol(Z), colnames(Z))
    gene.avg = sapply(IDX, function(idx) as.numeric(Matrix::rowMeans(Z[, idx])))
    rownames(gene.avg) = rownames(Z)

    max.col = apply(gene.avg, 1, which.max)
    max.val = apply(gene.avg, 1, max)
    IDX = split(1:nrow(Z), max.col)
    row.perm = as.numeric(unlist(sapply(IDX, function(idx) {
      return(idx[order(max.val[idx], decreasing = TRUE)])
    })))
    Z = Z[row.perm, ]
  }  
  

	# Visualize
	if(is.list(CPal)) {
		Pal = CPal[1:length(levels(Labels))]
	} else {
		Pal = ggpubr::get_palette(CPal, length(levels(Labels)))
	}
	names(Pal) = levels(Labels)

	ha_column = HeatmapAnnotation(df = list(Celltype = colnames(Z)), col = list(Celltype = Pal), annotation_legend_param = list(Celltype=list(title_gp = gpar(fontsize = 8), labels_gp = gpar(fontsize = 5))), which = "column")

	gradPal = grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(n = 9, name = "RdYlBu")))(100)

	ht = ComplexHeatmap::Heatmap(Z, row_names_gp = gpar(fontsize = 7), column_names_gp = gpar(fontsize = 0), top_annotation = ha_column, name = "z-score", cluster_rows = FALSE, cluster_columns = FALSE, col = gradPal)

	return(ht)
}


visualize.gene.summaries <- function(ACTIONet.out, sce, genes, Labels, alpha_val = 0.9, thread_no = 8, prune = FALSE, rescale = TRUE, expr.slot = "logcounts") {
  imputed.expression = impute.genes.using.ACTIONet(ACTIONet.out, sce, genes, alpha_val = alpha_val, thread_no = thread_no, prune = prune, rescale = rescale, expr.slot = "logcounts")
  
  Z = scale(imputed.expression)
  
  IDX = split(1:length(Labels), as.character(Labels))
  
  expr.sum = t(sapply(IDX, function(idx) {
    vals = as.numeric(Matrix::colMeans(Z[idx, ]))
  }))
  colnames(expr.sum) = colnames(imputed.expression)
  
  expr.sum[expr.sum < 0] = 0
  
  
  Heatmap(expr.sum, row_names_gp = gpar(fontsize =8), column_names_gp = gpar(fontsize = 8), cluster_rows = TRUE, cluster_columns = TRUE, col = blues9, raster_quality = 100)
}


filter.shared.genes <- function(sce, genes, Labels, AUC.threshold = 0.55, pval.threshold = 1e-10) {
	sce$groups = Labels
	wilc.out = presto::wilcoxauc(sce[genes, ], "groups") 
	
	
  ii = match(wilc.out$feature, genes)
  jj = match(wilc.out$group, levels(Labels))
  vv = -log10(wilc.out$padj)
  vv[wilc.out$auc < AUC.threshold] = 0
  logPval.mat = as.matrix(Matrix::sparseMatrix(i = ii, j = jj, x = vv, dims = c(length(genes), length(levels(Labels)))))
  rownames(logPval.mat) = genes
  colnames(logPval.mat) = levels(Labels)
  
  counts = Matrix::rowSums(logPval.mat > -log10(pval.threshold))
	
  selected.genes = genes[which(counts == 1)]
  
  perm = order(apply(logPval.mat[selected.genes, ], 1, which.max))
    
  return(selected.genes[perm])
}

