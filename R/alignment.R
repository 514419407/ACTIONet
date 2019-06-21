alignDatasets <- function(ACTIONet.out.ds1, sce.ds1, ACTIONet.out.ds2, sce.ds2, row.mappings, core.only = FALSE) {
	# Construct raw archetype profiles
	profile.ds1 = ACTIONet.out.ds1$reconstruct.out$archetype_profile[row.mappings[, 1], ]
	profile.ds2 = ACTIONet.out.ds2$reconstruct.out$archetype_profile[row.mappings[, 2], ]

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
	V.ds1 = V[row.mappings[, 1], ]

	## Dataset 2     
	Annot = rowData(sce.ds2)
	cnames = colnames(Annot)
	idx = grep('^PC', cnames)
	V = as.matrix(Annot[, idx])

	perm = order(sapply(cnames[idx], function(str) as.numeric(stringr::str_replace(str, "PC", ""))))
	V = V[, perm]
	V.ds2 = V[row.mappings[, 2], ]

	V.alignment = V.ds1 %*% t(V.ds1) %*% V.ds1

	# Align
	## Center profiles
	profile.ds1_centered = scale(profile.ds1)
	profile.ds2_centered = scale(profile.ds2)

	## Project
	S_r.ds1 = t(V.ds1) %*% profile.ds1_centered
	S_r.ds2 = t(V.alignment) %*% profile.ds2_centered

	# Compute similarities
	pairwise.similarity = cor(S_r.ds1, S_r.ds2)

	out.list = list(S_r.ds1 = S_r.ds1, V.ds1 = V.ds1, S_r.ds2 = S_r.ds2, V.ds2 = V.alignment, pairwise.similarity = pairwise.similarity)

	return(out.list)
}


visualize.pairwise.mapping.heatmap <-function(pairwise.mapping, ds1.arch.labels, ds2.arch.labels, CPal = "d3", ds1.name = "Dataset 1", ds2.name = "Dataset 2") {
	require(ComplexHeatmap)
	
	Sim = pairwise.mapping$pairwise.similarity
		
	if(!is.factor(ds1.arch.labels)) {
	  ds1.arch.labels = factor(ds1.arch.labels, levels = sort(unique(ds1.arch.labels)))
	}

	if(!is.factor(ds2.arch.labels)) {
	  ds2.arch.labels = factor(ds2.arch.labels, levels = sort(unique(ds2.arch.labels)))
	}

	rownames(Sim) = ds1.arch.labels
	colnames(Sim) = ds2.arch.labels
	Sim = Sim[order(ds1.arch.labels), order(ds2.arch.labels)]

	
	Sim = Sim[Matrix::rowSums(Sim != 0) > 0, Matrix::colSums(Sim != 0) > 0]
	
	all.levels = union(levels(ds1.arch.labels), levels(ds2.arch.labels))	
	Pal = ggpubr::get_palette(CPal, length(all.levels))
	names(Pal) = all.levels
	
	ha_row = rowAnnotation(df = data.frame(ds1.annotations = rownames(Sim)), col = list(ds1.annotations = Pal), annotation_legend_param = list(ds1.annotations=list(title_gp = gpar(fontsize = 8), labels_gp = gpar(fontsize = 5))))
	ha_column = HeatmapAnnotation(df = data.frame(ds2.annotations = colnames(Sim)), col = list(ds2.annotations = Pal), annotation_legend_param = list(ds2.annotations=list(title_gp = gpar(fontsize = 8), labels_gp = gpar(fontsize = 5))))

	ht_list = ComplexHeatmap::Heatmap(Sim, row_names_gp = gpar(fontsize =0), column_names_gp = gpar(fontsize = 0), top_annotation = ha_column, name = 'Mapping score', cluster_rows = FALSE, cluster_columns = FALSE, left_annotation = ha_row, col = RColorBrewer::brewer.pal(11, "RdYlBu")[seq(11, 1, by = -1)], row_title = ds1.name, column_title = ds2.name, raster_quality = 100)

	return(ht_list)	
}

match.pairwise.mapping <- function(pairwise.mapping, method = "MWM") {
	W = pairwise.mapping$pairwise.similarity
	W[W < 0] = 0
	
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
	
	alluvial(df[,1:2], freq=df[, 3], axis_labels = c(ds1.name, ds2.name), col = cols)
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
