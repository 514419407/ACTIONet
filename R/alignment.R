alignDatasets <- function(ACTIONet.out.ds1, sce.ds1, ACTIONet.out.ds2, sce.ds2, row.mappings, core.only = FALSE) {
	
	if(is.character(row.mappings[, 1])) {			
		X1 = ACTIONet.out.ds1$reconstruct.out$archetype_profile
		X2 = ACTIONet.out.ds2$reconstruct.out$archetype_profile
		
		rownames(X1) = toupper(rownames(X1))
		rownames(X2) = toupper(rownames(X2))
	
		row.mappings = toupper(row.mappings)
		
		mask = (row.mappings[, 1] %in% rownames(X1)) & (row.mappings[, 2] %in% rownames(X2))
		
		row.mappings = row.mappings[mask, ]

		ii = match(row.mappings[, 1], rownames(X1))
		jj = match(row.mappings[, 2], rownames(X2))
		
		if(dim(row.mappings)[2] == 3) { # Many-to-many mappings: Perform MWM
		}
	} else if(is.numeric(row.mappings[, 1])) {
		ii = row.mappings[, 1]
		jj = row.mappings[, 2]
	}
	
	# Construct raw archetype profiles	
	profile.ds1 = ACTIONet.out.ds1$reconstruct.out$archetype_profile[ii, ]
	profile.ds2 = ACTIONet.out.ds2$reconstruct.out$archetype_profile[jj, ]

	if(core.only == TRUE) {
		profile.ds1 = profile.ds1[, ACTIONet.out.ds1$core.out$core.archs]
		profile.ds2 = profile.ds2[, ACTIONet.out.ds2$core.out$core.archs]
	}
	
	# Recover V (PCA projection matrix)
	## Dataset 1
	Annot = rowData(sce.ds1)
	cnames = colnames(Annot)
	idx = grep('^PC', cnames)
	V = as.matrix(Annot[, idx])

	perm = order(sapply(cnames[idx], function(str) as.numeric(stringr::str_replace(str, "PC", ""))))
	V = V[, perm]
	V.ds1 = V[ii, ]

	## Dataset 2     
	Annot = rowData(sce.ds2)
	cnames = colnames(Annot)
	idx = grep('^PC', cnames)
	V = as.matrix(Annot[, idx])

	perm = order(sapply(cnames[idx], function(str) as.numeric(stringr::str_replace(str, "PC", ""))))
	V = V[, perm]
	V.ds2 = V[jj, ]

	V.alignment = V.ds1 %*% t(V.ds1) %*% V.ds1

	# Align
	## Center profiles
	profile.ds1_centered = scale(profile.ds1, scale = FALSE)
	profile.ds2_centered = scale(profile.ds2, scale = FALSE)

	## Project
	S_r.ds1 = t(V.ds1) %*% profile.ds1_centered
	S_r.ds2 = t(V.alignment) %*% profile.ds2_centered

	# Compute significance by creating pseudo-archetypes
	set.seed(0)

	C1 = ACTIONet.out.ds1$reconstruct.out$C_stacked
	IC1 = which(Matrix::rowSums(C1) > 0) ## Influential cells
	counts1 = sample(Matrix::colSums(C1 > 0), ncol(C1), replace = TRUE)
	C1.rand = as(sapply(counts1, function(c) {
	  ii = sample(IC1, size = c)
	  return(as.numeric(Matrix::sparseVector(i = ii, x = 1/length(c), length = nrow(C1))))
	}), 'sparseMatrix')
	profile1.rand = sce.ds1@assays[["logcounts"]] %*% C1.rand
	S_r.ds1.rand = t(V.ds1) %*% scale(profile1.rand, scale = FALSE)


	C2 = ACTIONet.out.ds2$reconstruct.out$C_stacked
	IC2 = which(Matrix::rowSums(C2) > 0) ## Influential cells
	counts2 = sample(Matrix::colSums(C2 > 0), ncol(C2), replace = TRUE)
	C2.rand = as(sapply(counts2, function(c) {
	  ii = sample(IC2, size = c)
	  return(as.numeric(Matrix::sparseVector(i = ii, x = 1/length(c), length = nrow(C2))))
	}), 'sparseMatrix')
	profile2.rand = sce.ds2@assays[["logcounts"]] %*% C2.rand
	S_r.ds2.rand = t(V.alignment) %*% scale(profile2.rand, scale = FALSE)

	# Compute similarities
	pairwise.similarity = cor(S_r.ds1, S_r.ds2)	
	pairwise.similarity.rand = cor(S_r.ds1.rand, S_r.ds2.rand)
	
	pairwise.similarity.z = (pairwise.similarity - median(pairwise.similarity.rand)) / mad(pairwise.similarity.rand)	


	A = cor(S_r.ds1)	
	A.rand = cor(S_r.ds1.rand)
	A.z = (A - median(A.rand)) / mad(A.rand)

	B = cor(S_r.ds2)	
	B.rand = cor(S_r.ds2.rand)
	B.z = (B - median(B.rand)) / mad(B.rand)

	
	# Output list
	out.list = list(S_r.ds1 = S_r.ds1, V.ds1 = V.ds1, S_r.ds2 = S_r.ds2, V.ds2 = V.alignment, pairwise.similarity = pairwise.similarity, pairwise.similarity.z = pairwise.similarity.z, A = A, A.z = A.z, B = B, B.z = B.z)

	return(out.list)
}

netAlign.backbone <- function(ACTIONet.out.ds1, sce.ds1, ACTIONet.out.ds2, sce.ds2, row.mappings, core.only = FALSE) {

	# Recover V (PCA projection matrix)
	## Dataset 1
	Annot = rowData(sce.ds1)
	cnames = colnames(Annot)
	idx = grep('^PC', cnames)
	V = as.matrix(Annot[, idx])

	perm = order(sapply(cnames[idx], function(str) as.numeric(stringr::str_replace(str, "PC", ""))))
	V = V[, perm]
	V.ds1 = V[ii, ]

	## Dataset 2     
	Annot = rowData(sce.ds2)
	cnames = colnames(Annot)
	idx = grep('^PC', cnames)
	V = as.matrix(Annot[, idx])

	perm = order(sapply(cnames[idx], function(str) as.numeric(stringr::str_replace(str, "PC", ""))))
	V = V[, perm]
	V.ds2 = V[jj, ]

	V.alignment = V.ds1 %*% t(V.ds1) %*% V.ds1

	# Align
	## Center profiles
	profile.ds1_centered = scale(profile.ds1, scale = FALSE)
	profile.ds2_centered = scale(profile.ds2, scale = FALSE)

	## Project
	S_r.ds1 = t(V.ds1) %*% profile.ds1_centered
	S_r.ds2 = t(V.alignment) %*% profile.ds2_centered

	# Compute similarities
	pairwise.similarity = cor(S_r.ds1, S_r.ds2)

	out.list = list(S_r.ds1 = S_r.ds1, V.ds1 = V.ds1, S_r.ds2 = S_r.ds2, V.ds2 = V.alignment, pairwise.similarity = pairwise.similarity)
	
}

visualize.pairwise.mapping.heatmap <-function(pairwise.mapping, ds1.arch.labels = NA, ds2.arch.labels = NA, CPal = "d3", ds1.name = "Dataset 1", ds2.name = "Dataset 2") {
	require(ComplexHeatmap)
	require(seriation)
	
	Sim = pairwise.mapping$pairwise.similarity
	Sim[is.na(Sim)] = 0

	
	#Dr = Matrix::Diagonal(nrow(Sim), 1 / sqrt(as.numeric(Matrix::rowSums(Sim))))
	#Dc = Matrix::Diagonal(ncol(Sim), 1 / sqrt(as.numeric(Matrix::colSums(Sim))))
	#Sim = as.matrix(Dr %*% Sim %*% Dc)
	
	#MED = median(Sim)
	#MAD = mad(as.numeric(Sim))
	
	#Sim = (Sim - MED) / MAD
	#Sim = 1 / (1+exp(-0.5*Sim))
	
	
	if(!is.factor(ds1.arch.labels)) {
	  ds1.arch.labels = factor(ds1.arch.labels, levels = sort(unique(ds1.arch.labels)))
	}

	if(!is.factor(ds2.arch.labels)) {
	  ds2.arch.labels = factor(ds2.arch.labels, levels = sort(unique(ds2.arch.labels)))
	}

	rownames(Sim) = ds1.arch.labels
	colnames(Sim) = ds2.arch.labels
	Sim = Sim[order(ds1.arch.labels), order(ds2.arch.labels)]


	#perm = seriation::get_order(seriation::seriate(dist(t(Sim))))
	#Sim = Sim[, perm]
	
	#Sim = Sim[Matrix::rowSums(Sim != 0) > 0, Matrix::colSums(Sim != 0) > 0]
	
	all.levels = union(levels(ds1.arch.labels), levels(ds2.arch.labels))	
	Pal = ggpubr::get_palette(CPal, length(all.levels))
	names(Pal) = all.levels
	
	ha_row = rowAnnotation(df = data.frame(ds1.annotations = rownames(Sim)), col = list(ds1.annotations = Pal), annotation_legend_param = list(ds1.annotations=list(title_gp = gpar(fontsize = 8), labels_gp = gpar(fontsize = 5))))
	ha_column = HeatmapAnnotation(df = data.frame(ds2.annotations = colnames(Sim)), col = list(ds2.annotations = Pal), annotation_legend_param = list(ds2.annotations=list(title_gp = gpar(fontsize = 8), labels_gp = gpar(fontsize = 5))))

	gradPal = grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(n = 9, name = "RdYlBu")))(100)
	ht_list = ComplexHeatmap::Heatmap(Sim, row_names_gp = gpar(fontsize =0), column_names_gp = gpar(fontsize = 0), top_annotation = ha_column, name = 'Mapping score', cluster_rows = FALSE, cluster_columns = FALSE, left_annotation = ha_row, col = gradPal, row_title = ds1.name, column_title = ds2.name, raster_quality = 100)

	return(ht_list)	
}


visualize.homologous.cell.states <-function(pairwise.mapping, AG.graph, ds1.annotations = NA, ds2.annotations = NA, CPal = "d3", ds1.name = "Dataset 1", ds2.name = "Dataset 2") {
	require(ComplexHeatmap)
	require(seriation)
	
	Sim = pairwise.mapping$pairwise.similarity
	Sim[is.na(Sim)] = 0

	
	#Dr = Matrix::Diagonal(nrow(Sim), 1 / sqrt(as.numeric(Matrix::rowSums(Sim))))
	#Dc = Matrix::Diagonal(ncol(Sim), 1 / sqrt(as.numeric(Matrix::colSums(Sim))))
	#Sim = as.matrix(Dr %*% Sim %*% Dc)
	
	#MED = median(Sim)
	#MAD = mad(as.numeric(Sim))
	
	#Sim = (Sim - MED) / MAD
	#Sim = 1 / (1+exp(-0.5*Sim))
	
	
	if(!is.factor(ds1.arch.labels)) {
	  ds1.arch.labels = factor(ds1.arch.labels, levels = sort(unique(ds1.arch.labels)))
	}

	if(!is.factor(ds2.arch.labels)) {
	  ds2.arch.labels = factor(ds2.arch.labels, levels = sort(unique(ds2.arch.labels)))
	}

	rownames(Sim) = ds1.arch.labels
	colnames(Sim) = ds2.arch.labels
	Sim = Sim[order(ds1.arch.labels), order(ds2.arch.labels)]


	#perm = seriation::get_order(seriation::seriate(dist(t(Sim))))
	#Sim = Sim[, perm]
	
	#Sim = Sim[Matrix::rowSums(Sim != 0) > 0, Matrix::colSums(Sim != 0) > 0]
	
	all.levels = union(levels(ds1.arch.labels), levels(ds2.arch.labels))	
	Pal = ggpubr::get_palette(CPal, length(all.levels))
	names(Pal) = all.levels
	
	ha_row = rowAnnotation(df = data.frame(ds1.annotations = rownames(Sim)), col = list(ds1.annotations = Pal), annotation_legend_param = list(ds1.annotations=list(title_gp = gpar(fontsize = 8), labels_gp = gpar(fontsize = 5))))
	ha_column = HeatmapAnnotation(df = data.frame(ds2.annotations = colnames(Sim)), col = list(ds2.annotations = Pal), annotation_legend_param = list(ds2.annotations=list(title_gp = gpar(fontsize = 8), labels_gp = gpar(fontsize = 5))))

	gradPal = grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(n = 9, name = "RdYlBu")))(100)
	ht_list = ComplexHeatmap::Heatmap(Sim, row_names_gp = gpar(fontsize =0), column_names_gp = gpar(fontsize = 0), top_annotation = ha_column, name = 'Mapping score', cluster_rows = FALSE, cluster_columns = FALSE, left_annotation = ha_row, col = gradPal, row_title = ds1.name, column_title = ds2.name, raster_quality = 100)

	return(ht_list)	
}


match.pairwise.mapping <- function(pairwise.mapping, method = "MWM", z.threshold = 1) {
	W = pairwise.mapping$pairwise.similarity.z
	W[W < z.threshold] = 0
	
	if(method == "MWM") {
		W.matched = MWM(W)
	} 
	else if(method == "MNN") {
		row.max.idx = apply(W, 1, which.max)
		W.forward.match = sparseMatrix(1:dim(W)[1], row.max.idx, dims = dim(W))

		col.max.idx = apply(W, 2, which.max)
		W.backward.match = sparseMatrix(col.max.idx, 1:dim(W)[2], dims = dim(W))

		W.matched = W * W.forward.match * W.backward.match		
	}
	
	W.matched = as(W.matched, 'dgTMatrix')
	df = data.frame(ii = W.matched@i + 1, jj = W.matched@j + 1, score = W.matched@x)		
	
	return(df)
}

visualize.pairwise.mapping.alluvial <-function(pairwise.mapping, ds1.arch.labels, ds2.arch.labels, CPal = "d3", ds1.name = "Dataset 1", ds2.name = "Dataset 2", matching.method = "MWM") {
	require(alluvial)

	matches = match.pairwise.mapping(pairwise.mapping, matching.method)
	df = data.frame(L1 = ds1.arch.labels[matches$ii], L2 = ds2.arch.labels[matches$jj], freq = 100*matches$score/sum(matches$score), stringsAsFactors = FALSE);
	
	all.levels = union(levels(ds1.arch.labels), levels(ds2.arch.labels))	
	Pal = ggpubr::get_palette(CPal, length(all.levels))
	names(Pal) = all.levels

	cols = rgb(t((col2rgb(Pal[df$L1]) + col2rgb(Pal[df$L2]))/512))
	
	alluvial(df[,1:2], freq=df[, 3], axis_labels = c(ds1.name, ds2.name), col = cols, gap.width=0.1)
}


compute.pairwise.gene.diff <-function(pairwise.mapping, ds1.arch.labels = NA, ds2.arch.labels = NA, matching.method = "MNN", prune = TRUE, nonparametric = FALSE, score.threshold = 0.5) {
	matches = match.pairwise.mapping(pairwise.mapping, matching.method)
	if(length(ds1.arch.labels) == 1 || length(ds2.arch.labels) == 1) {
		prune = FALSE
	}
	if( prune == TRUE ) {
		df = data.frame(L1 = ds1.arch.labels[matches$ii], ii = matches$ii, L2 = ds2.arch.labels[matches$jj], jj = matches$jj, score = matches$score, stringsAsFactors = FALSE)
		df = df[df$L1 == df$L2, ]
		Labels = sapply(1:dim(df)[1], function(x) sprintf('%s', df$L1[x]))
	}
	else {
		df = data.frame(ii = matches$ii, jj = matches$jj, score = matches$score)
		Labels = sapply(1:dim(df)[1], function(x) sprintf('Pair %d: %d <-> %d', x, df$ii[x], df$jj[x]))
	}
	df = df[df$score > score.threshold, ]
	
	S1 = CTL_vs_SCZ_mapping$V.ds1 %*% CTL_vs_SCZ_mapping$S_r.ds1[, df$ii]
	S2 = CTL_vs_SCZ_mapping$V.ds2 %*% CTL_vs_SCZ_mapping$S_r.ds2[, df$ii]

	if(nonparametric == TRUE) {
		S1 = apply(S1, 2, function(v) rank(v)/length(v))
		S2 = apply(S2, 2, function(v) rank(v)/length(v))
	}
	
	delta.r = S1 - S2
	colnames(delta.r) = Labels
	
	return(delta.r)
}


visualize.alignment.graph <- function(pairwise.mapping, ds1.arch.labels = NA, ds2.arch.labels = NA, CPal = "d3", ds1.name = "Dataset 1", ds2.name = "Dataset 2") {
	require(ComplexHeatmap)
	X1 = orthoProject(pairwise.mapping$S_r.ds1, rowMeans(pairwise.mapping$S_r.ds1))
	X2 = orthoProject(pairwise.mapping$S_r.ds2, rowMeans(pairwise.mapping$S_r.ds2))

	ds1.W = cor(X1); ds1.W = ds1.W / max(ds1.W)
	ds2.W = cor(X2); ds2.W = ds2.W / max(ds2.W)
	ds12.W = pairwise.mapping$pairwise.similarity; ds12.W = ds12.W / max(ds12.W)
	W1 = cbind(ds1.W, ds12.W)
	W2 = cbind(t(ds12.W), ds2.W)
	W = rbind(W1, W2)

	Labels = c(as.character(ds1.arch.labels), as.character(ds2.arch.labels))
	dsLabels = c(rep(ds1.name, length(ds1.arch.labels)), rep(ds2.name, length(ds2.arch.labels)))  


	all.levels = sort(unique(Labels))
	Pal = ggpubr::get_palette(CPal, length(all.levels))
	names(Pal) = all.levels

	dsPal = toupper(c("#d95f02", "#7570b3"))
	names(dsPal) = c(ds1.name, ds2.name)

	ha_row = rowAnnotation(df = data.frame(Type = Labels, Dataset = dsLabels), col = list(Type = Pal, Dataset = dsPal), annotation_legend_param = list(Type=list(title_gp = gpar(fontsize = 8), labels_gp = gpar(fontsize = 5)), Dataset=list(title_gp = gpar(fontsize = 8), labels_gp = gpar(fontsize = 5)) ))


	ha_column = HeatmapAnnotation(df = data.frame(Type = Labels, Dataset = dsLabels), col = list(Type = Pal, Dataset = dsPal), annotation_legend_param = list(Type=list(title_gp = gpar(fontsize = 8), labels_gp = gpar(fontsize = 5)), Dataset=list(title_gp = gpar(fontsize = 8), labels_gp = gpar(fontsize = 5)) ))


	ht_list = ComplexHeatmap::Heatmap(W, row_names_gp = gpar(fontsize =0), column_names_gp = gpar(fontsize = 0), top_annotation = ha_column, name = 'Mapping score', cluster_rows = TRUE, cluster_columns = TRUE, left_annotation = ha_row, col = RColorBrewer::brewer.pal(11, "RdYlBu")[seq(11, 1, by = -1)], row_title = ds1.name, column_title = ds2.name, raster_quality = 100)

	return(ht_list)
}	


constructAlignmentGraph <- function(pairwise.mapping, method = "MWM", z.threshold = 1) {
	require(seriation)
	require(ACTIONet)
	
	alignment.df = match.pairwise.mapping(pairwise.mapping, method = method)
	alignment.df = alignment.df[alignment.df$score > 0, ]

	A = pairwise.mapping$A.z
	B = pairwise.mapping$B.z

	X = A[alignment.df$ii, alignment.df$ii]
	Y = B[alignment.df$jj, alignment.df$jj]
	
	AG = (X + Y) / sqrt(2) # Should be still "normal"
	
	diag(AG) = 0

	# stats = t(scale(X)) %*% Diagonal(nrow(AG), alignment.df$score) %*% scale(Y)
	# 
	# stats = as.numeric(diag(t(scale(X)) %*% Diagonal(nrow(AG), alignment.df$score) %*% scale(Y)))
	# rand.stats = sapply(1:1000, function(i) {
	# 	rand.stat = as.numeric(diag(t(scale(X)) %*% Diagonal(nrow(AG), sample(alignment.df$score)) %*% scale(Y)))    
	# })	
	# Mu = apply(rand.stats, 1, mean)
	# Sigma = apply(rand.stats, 1, sd)
	# z = (stats - Mu) / Sigma  
  
	# mask = z.threshold < z
	# 
	# AG = AG[mask, mask]
	# alignment.df = alignment.df[mask, ]
	# z = z[mask]


	aligned.modules = signed_cluster(as(AG, 'sparseMatrix'), resolution_parameter = 1.8, seed = 0)

	AG.graph = graph_from_adjacency_matrix(AG, weighted = TRUE, mode = "undirected")
	V(AG.graph)$module = aligned.modules
	
	V(AG.graph)$ds1.arch = alignment.df$ii
	V(AG.graph)$ds2.arch = alignment.df$jj
	V(AG.graph)$alignment.score = alignment.df$score
	#V(AG.graph)$neighborhood.alignment.score = z
  
	# Arrange vertices to focus modularity of cell states
	IDX = split(1:length(aligned.modules), aligned.modules)

	# Order within each module
	pruned.modules = sapply(IDX, function(idx) {
	  if(length(idx) == 1) {
		return(idx)
	  } else {
		X = AG[idx, idx]
		pp = get_order(seriate(dist(X)))

		return(idx[pp])
	  }
	})

	# Order modules themselves
	cent = sapply(pruned.modules, function(idx) {
	  if(length(idx) == 1) {
		return(AG[idx, ])
	  } else {
		return(Matrix::colMeans(AG[idx, ]))
	  }
	})
	module.order = get_order(seriate(dist(t(cent))))

	perm = unlist(pruned.modules[module.order])

	V(AG.graph)$order = perm
	
	return(AG.graph)  
}

visualizeAlignmentGraph.heatmap <- function(AG.graph, annotations.df, CPal = "d3") {
	require(ComplexHeatmap)
	require(RColorBrewer)
	
	W = as.matrix(get.adjacency(AG.graph, attr = "weight"))
	perm = V(AG.graph)$order 

	W = W[perm, perm]
	annotations.df = annotations.df[perm, ]

	annotations = colnames(annotations.df)
	Pals = lapply(annotations, function(annotation){
		L = annotations.df[, annotation]
		if(!is.factor(L)) {
		  L = factor(L, levels = sort(unique(L)))
		}

		all.levels = levels(L)
		Pal = ggpubr::get_palette(CPal, length(levels(L)))
		names(Pal) = levels(L)
		return(Pal)
	})
	names(Pals) = annotations

	ha_row = HeatmapAnnotation(df = annotations.df, col = Pals, show_legend = FALSE, show_annotation_name = FALSE, which = "row")

	ha_column = HeatmapAnnotation(df = annotations.df, col = Pals, show_legend = TRUE, show_annotation_name = TRUE, which = "column")

	gradPal = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)

	ht = ComplexHeatmap::Heatmap(W, row_names_gp = gpar(fontsize =0), column_names_gp = gpar(fontsize = 0), top_annotation = ha_column, name = "z-score", cluster_rows = FALSE, cluster_columns = FALSE, left_annotation = ha_row, col = gradPal, raster_quality = 100)

	return(ht)
}


visualizeAlignmentGraph.alluvial <- function(AG.graph, ds1.arch.labels, ds2.arch.labels, CPal = "d3", ds1.name = "Dataset 1", ds2.name = "Dataset 2") {
	require(alluvial)

	A = as(get.adjacency(AG, attr = "weight"), 'dgTMatrix')

	Sim = pairwise.mapping$pairwise.similarity
	ii = as.numeric(V(AG)$ds1.arch)
	jj = as.numeric(V(AG)$ds2.arch)
	idx = ii + (jj-1)*nrow(Sim)
	vv = Sim[idx]

	df = data.frame(src = ds1.arch.labels[ii], dst = ds2.arch.labels[jj], weight = vv, stringsAsFactors = FALSE)

	Pal1 = ggpubr::get_palette(CPal, length(levels(ds1.arch.labels)))
	names(Pal1) = levels(ds1.arch.labels)

	Pal2 = ggpubr::get_palette(CPal, length(levels(ds2.arch.labels)))
	names(Pal2) = levels(ds2.arch.labels)


	alluvial(df[,1:2], freq=df[, 3], axis_labels = c(ds1.name, ds2.name), col = Pal1[df$src], border = Pal2[df$dst], gap.width=0.1)
}




visualizeAlignmentGraph.Sankey <- function(AG.graph, pairwise.mapping, ds1.arch.labels, ds2.arch.labels) {
	library(networkD3)

	A = as(get.adjacency(AG.graph, attr = "weight"), 'dgTMatrix')

	Sim = pairwise.mapping$pairwise.similarity
	ii = as.numeric(V(AG.graph)$ds1.arch)
	jj = as.numeric(V(AG.graph)$ds2.arch)
	idx = ii + (jj-1)*nrow(Sim)
	vv = Sim[idx]

	Nodes = data.frame(Name = c(levels(ds1.arch.labels), levels(ds2.arch.labels)))
	
	Links = data.frame(Source = as.numeric(ds1.arch.labels[ii]) - 1, Target = length(levels(ds1.arch.labels)) + as.numeric(ds2.arch.labels[jj]) - 1, Value = vv)
	
	
	Adj = Matrix::sparseMatrix(i = Links$Source+1, j = Links$Target+1, x = 1, dims = c(nrow(Nodes), nrow(Nodes))) 
	
	filter.mask = (Matrix::rowSums(Adj) + Matrix::colSums(Adj)) == 0  
  
	Nodes = data.frame(Name = as.character(Nodes$Name[!filter.mask]), stringsAsFactors = FALSE)
	Adj = as(Adj[!filter.mask, !filter.mask], 'dgTMatrix')
	
	Links = data.frame(Source = Adj@i, Target = Adj@j, Value = Adj@x)
	

  
	sn = sankeyNetwork(Links = Links, Nodes = Nodes, Source = 'Source',
             Target = 'Target', Value = 'Value', NodeID = 'Name',
             units = 'TWh', fontSize = 12, nodeWidth = 30)  
  
	return(sn)
}


visualize.pairwise.alignment <- function(pairwise.mapping, AG.graph, ds1.annotations = NULL, ds2.annotations = NULL, ds1.name = "Dataset 1", ds2.name = "Dataset 2", CPal = "d3") {
  if(! (is.null(ds1.annotations)) )
    colnames(ds1.annotations) = sapply(colnames(ds1.annotations), function(x) sprintf('%s_%s', ds1.name, x))

  if(! (is.null(ds2.annotations)) )
    colnames(ds2.annotations) = sapply(colnames(ds2.annotations), function(x) sprintf('%s_%s', ds2.name, x))
  
  ii = V(AG.graph)$ds1.arch
  jj = V(AG.graph)$ds2.arch
  Sim = pairwise.mapping$pairwise.similarity.z[ii, jj]
  
  ds1.annotations = ds1.annotations[ii, , drop = FALSE]
  ds2.annotations = ds2.annotations[jj, , drop = FALSE]
  
  if(! (is.null(ds1.annotations)) ) {
    ds1.annotations$Module = V(AG.graph)$module
  } else
  ds1.annotations = data.frame(Module = V(AG.graph)$module)
  
  #ds1.annotations$Consistency = V(AG.graph)$alignment.score #1 + floor(9.9 / (1 + exp(-2*scale(V(AG.graph)$alignment.score))))
  ds1.annotations$Consistency = 1 / (1 + exp(-3*scale(V(AG.graph)$alignment.score)))
  
  
  perm = V(AG.graph)$order 
  Sim = Sim[perm, perm]
  
  ds1.annotations = ds1.annotations[perm, , drop = FALSE]
  ds2.annotations = ds2.annotations[perm, , drop = FALSE]
  ii = ii[perm]
  jj = jj[perm]
  
  
  
  annotations = colnames(ds1.annotations)
  Pals.ds1 = lapply(annotations, function(annotation){
	  L = ds1.annotations[, annotation]
	  if(!is.factor(L)) {
		L = factor(L, levels = sort(unique(L)))
	  }
	  
	  all.levels = levels(L)
	  Pal = ggpubr::get_palette(CPal, length(levels(L)))
	  names(Pal) = levels(L)
	  return(Pal)
  })
  names(Pals.ds1) = annotations
  
  Pal = ggpubr::get_palette('igv', length(sort(unique(ds1.annotations$Module))))
  names(Pal) = sort(unique(ds1.annotations$Module))
  Pals.ds1$Module = Pal
  

  Consistency_col_fun = circlize::colorRamp2(seq(0, 1, length.out = 9), RColorBrewer::brewer.pal(n = 9, name = "Blues"))  
  Pals.ds1$Consistency = Consistency_col_fun
  
  
  annotations = colnames(ds2.annotations)
  Pals.ds2 = lapply(annotations, function(annotation){
	  L = ds2.annotations[, annotation]
	  if(!is.factor(L)) {
		L = factor(L, levels = sort(unique(L)))
	  }
	  
	  all.levels = levels(L)
	  Pal = ggpubr::get_palette(CPal, length(levels(L)))
	  names(Pal) = levels(L)
	  return(Pal)
  })
  names(Pals.ds2) = annotations
  

  
  ha_row = HeatmapAnnotation(df = ds1.annotations, col = Pals.ds1, show_legend = TRUE, show_annotation_name = TRUE, which = "row")
  
  if(! is.null(ds2.annotations))
    ha_column = HeatmapAnnotation(df = ds2.annotations, col = Pals.ds2, show_legend = TRUE, show_annotation_name = TRUE, annotation_name_side = "left", which = "column")
  
  gradPal = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)
  
  if(! (is.null(ds2.annotations))) {
    ht = ComplexHeatmap::Heatmap(Sim, row_names_gp = gpar(fontsize = 0), column_names_gp = gpar(fontsize = 0), top_annotation = ha_column, name = "Score", cluster_rows = FALSE, cluster_columns = FALSE, left_annotation = ha_row, col = gradPal, raster_quality = 100, row_title = ds1.name, column_title = ds2.name)
  } else {
    ht = ComplexHeatmap::Heatmap(Sim, row_names_gp = gpar(fontsize =0), column_names_gp = gpar(fontsize = 0), name = "Score", cluster_rows = FALSE, cluster_columns = FALSE, left_annotation = ha_row, col = gradPal, raster_quality = 100, row_title = ds1.name, column_title = ds2.name)
    
  }
  return(ht)
}
