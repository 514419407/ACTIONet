reduce.sce <- function(sce, reduced_dim = 50, max.iter = 5) {
	require(scran)
	require(scater)
	require(ACTIONet)

	if( !("logcounts" %in% names(sce@assays)) ) {
		sce.norm = sce    
		A = as(sce@assays[["counts"]], 'dgTMatrix')
		cs = Matrix::colSums(A)    
		B = Matrix::sparseMatrix(i = A@i+1, j = A@j+1, x = log(1 + median(cs)*(A@x / cs[A@j + 1])), dims = dim(A))
		sce.norm@assays[["logcounts"]] = B   
	} else {
		sce.norm = sce
	}
	
	set.seed(0)
	reduction.out = reduceGeneExpression(as(sce.norm@assays[["logcounts"]], 'sparseMatrix'), reduced_dim = reduced_dim, method = 1, iters = max.iter)    

	SingleCellExperiment::reducedDim(sce.norm, "S_r") <- t(reduction.out$S_r)

	V = reduction.out$V
	colnames(V) = sapply(1:dim(V)[2], function(i) sprintf('PC%d', i))

	X = rowData(sce.norm)
	PC.idx = -grep('^PC', colnames(X))
	if(length(PC.idx) > 0)
	X = X[, PC.idx]	
	rowData(sce.norm) = cbind(V, X)

	return(sce.norm)  
}

batch.correct.sce.Harmony <- function(sce, batch.attr) {
	require(harmony)
	sce@reducedDims$S_r = harmony::HarmonyMatrix(sce@reducedDims$S_r, batch.attr, do_pca=FALSE)	
	return(sce)
}


reduce.and.batch.correct.sce.Harmony <- function(sce, batch.attr=NA, reduced_dim = 50) {
	if(is.na(batch.attr)) {
		print("You need to provide the batch vector/attr");
		return(sce)
	}

	set.seed(0)
	sce = reduce.sce(sce, reduced_dim = reduced_dim)
	batch.correct.sce.Harmony(sce, batch.attr)
}

reduce.and.batch.correct.sce.MNN <- function(sce, batch.attr=NA, reduced_dim = 50, k=20) {
	require(scran)
	require(scater)
	require(ACTIONet)

	if(length(batch.attr) ==  1 && is.na(batch.attr)) {
		print("You need to provide the batch vector/attr");
		return(sce)
	}

	sce = clearSpikes(sce)
	sce@assays[["counts"]] = as(sce@assays[["counts"]], 'sparseMatrix')

	if( length(batch.attr) ==  1) {
		IDX = split(1:dim(sce)[2], colData(sce)[[batch.attr]])	  
	}
	else {
		IDX = split(1:dim(sce)[2], batch.attr)	  
	}

	sce.list = lapply(IDX, function(idx) suppressWarnings(scater::normalize(sce[, idx])))
	sce.list.norm = do.call(scran::multiBatchNorm, sce.list)

	# Sort based on "complexity"
	perm = order(sapply(sce.list.norm, function(sce) dim(sce)[2]), decreasing = TRUE)    
	sce.list.norm = sce.list.norm[perm]

	sce.norm = do.call(SingleCellExperiment::cbind, sce.list.norm)
	sce.norm@assays[["logcounts"]] = as(sce.norm@assays[["logcounts"]], 'sparseMatrix')

	set.seed(0)
	mnn.out <- do.call(scran::fastMNN, c(sce.list.norm, list(k=20, d=reduced_dim, auto.order = FALSE, approximate=TRUE, cos.norm = FALSE)))

	SingleCellExperiment::reducedDim(sce.norm, "S_r") <- mnn.out$corrected

	V = mnn.out$rotation

	colnames(V) = sapply(1:dim(V)[2], function(i) sprintf('PC%d', i))

	X = rowData(sce.norm)
	PC.idx = grep('^PC', colnames(X))
	if(length(PC.idx) > 0)
		X = X[, -PC.idx]	
	
	rowData(sce.norm) = cbind(V, X)    

	return(sce.norm)
}


construct.signature.profile <- function(sce, ACTIONet.out) {
	require(SingleCellExperiment)

	# eigengene x archetypes
	reduced.archetype.profile = (t(sce@reducedDims[["S_r"]]) %*% ACTIONet.out$reconstruct.out$C_stacked)
	
	X = rowData(sce)
	cnames = colnames(X)
	idx = grep('^PC', cnames)
	V = as.matrix(X[, idx])
	
	perm = order(sapply(cnames[idx], function(str) as.numeric(stringr::str_replace(str, "PC", ""))))
	V = V[, perm]
	
	# gene x archetypes
	archetype.signature.profile = V %*% reduced.archetype.profile 
	rownames(archetype.signature.profile) = rownames(sce)
	
	return(archetype.signature.profile)
}

identify.core.archetypes <- function(ACTIONet.out, pruning.zscore.threshold = 3) {
	require(igraph)
	require(ACTIONet)
	
	mergeTree = mergeArchetypes(C_stacked = ACTIONet.out$reconstruct.out$C_stacked, H_stacked = ACTIONet.out$reconstruct.out$H_stacked)

	A = matrix(as.numeric(mergeTree < pruning.zscore.threshold & mergeTree != 0), ncol = ncol(mergeTree))
	A = A + t(A)
	merged.graph = graph_from_adjacency_matrix(A, weighted = TRUE, mode = "undirected")

	comps = components(merged.graph)
	core.archs = sapply(unique(comps$membership), function(c) min(which(comps$membership == c)))
	
	out.list = list(core.archs = core.archs, arch.membership = comps$membership, mergeTree=mergeTree)
	
	return(out.list)
}


run.ACTIONet <- function(sce, k_max = 20, compactness_level = 50, thread_no = 8, epsilon = 3.0, LC = 1.0, auto_adjust_LC = FALSE, arch.specificity.z = -1, core.z = 3.0, n_epochs = 500, sce.data.attr = "logcounts") {
	require(Matrix)
	require(igraph)
	require(ACTIONet)

	if( !(sce.data.attr %in% names(sce@assays)) ) {
		R.utils::printf("Attribute %s is not an assay of the input sce\n", sce.data.attr)
		return();
	}

	# Run ACTION
	ACTION.out = runACTION(t(sce@reducedDims[["S_r"]]), k_max = k_max, thread_no = thread_no)

	# Reconstruct archetypes in the original space  
	reconstruct.out = reconstructArchetypes(as(sce@assays[[sce.data.attr]], 'sparseMatrix'), ACTION.out$C, ACTION.out$H, z_threshold = arch.specificity.z)
	rownames(reconstruct.out$archetype_profile) = rownames(sce)

	# Build ACTIONet
	build.out = buildAdaptiveACTIONet(H_stacked = reconstruct.out$H_stacked, thread_no = thread_no, LC = LC, auto_adjust_LC = auto_adjust_LC, epsilon = epsilon)

	# Layout ACTIONet
	vis.out = layoutACTIONet(build.out$ACTIONet, S_r = t(sce@reducedDims[["S_r"]]), compactness_level = compactness_level, n_epochs = n_epochs)

	# Construct igraph object
	ACTIONet = graph_from_adjacency_matrix(build.out$ACTIONet, mode="undirected", weighted = TRUE)
	coor = vis.out$coordinates;
	coor3D = vis.out$coordinates_3D
	V(ACTIONet)$x = coor[, 1]
	V(ACTIONet)$y = coor[, 2]
	V(ACTIONet)$x3D = coor3D[, 1]
	V(ACTIONet)$y3D = coor3D[, 2]
	V(ACTIONet)$z3D = coor3D[, 3]
	V(ACTIONet)$color = rgb(vis.out$colors)

	arch.Lab = t(reconstruct.out$C_stacked) %*% grDevices::convertColor(color= vis.out$colors, from = 'sRGB', to = 'Lab')
	arch.colors = rgb(grDevices::convertColor(color=arch.Lab, from = 'Lab', to = 'sRGB'))
	arch.coordinates = t(reconstruct.out$C_stacked) %*% vis.out$coordinates
	arch.coordinates_3D = t(reconstruct.out$C_stacked) %*% vis.out$coordinates_3D
	arch.vis.out = list(colors = arch.colors, coordinates = arch.coordinates, coordinates_3D = arch.coordinates_3D)

	ACTIONet.out = list(ACTION.out = ACTION.out, reconstruct.out = reconstruct.out, build.out = build.out, vis.out = vis.out, ACTIONet = ACTIONet, arch.vis.out = arch.vis.out)

	# Add signature profile
	signature.profile = construct.signature.profile(sce = sce, ACTIONet.out = ACTIONet.out) 
	ACTIONet.out$signature.profile = signature.profile
		
	core.out = identify.core.archetypes(ACTIONet.out, core.z)
	H.core = runsimplexRegression(t(sce@reducedDims[["S_r"]]) %*% ACTIONet.out$reconstruct.out$C_stacked[, core.out$core.archs], t(sce@reducedDims[["S_r"]]))
	core.out$H = H.core	
	ACTIONet.out$core.out = core.out
	
	return(ACTIONet.out)
}

filter.ACTIONet.factors <- function(ACTIONet.out, sce, k_max = 20, compactness_level = 50, thread_no = 8, epsilon = 3.0, LC = 1.0, auto_adjust_LC = FALSE, arch.specificity.z = -1, core.z = 3.0, n_epochs = 500, sce.data.attr = "logcounts") {
	C = ACTIONet.out$ACTION.out$C
	H = ACTIONet.out$ACTION.out$H
	
	k_actual = length(C)
	if(k_actual <= k_max) {
		print("Requested k is less than or equal to the max k");
		return(ACTIONet.out)
	}
	
	ACTION.out = list(C = C[1:k_max], H = H[1:k_max])
	

	# Reconstruct archetypes in the original space  
	reconstruct.out = reconstructArchetypes(as(sce@assays[[sce.data.attr]], 'sparseMatrix'), ACTION.out$C, ACTION.out$H, z_threshold = arch.specificity.z)
	rownames(reconstruct.out$archetype_profile) = rownames(sce)

	# Build ACTIONet
	build.out = buildAdaptiveACTIONet(H_stacked = reconstruct.out$H_stacked, thread_no = thread_no, LC = LC, auto_adjust_LC = auto_adjust_LC, epsilon = epsilon)

	# Layout ACTIONet
	vis.out = layoutACTIONet(build.out$ACTIONet, S_r = t(sce@reducedDims[["S_r"]]), compactness_level = compactness_level, n_epochs = n_epochs)

	# Construct igraph object
	ACTIONet = graph_from_adjacency_matrix(build.out$ACTIONet, mode="undirected", weighted = TRUE)
	coor = vis.out$coordinates;
	coor3D = vis.out$coordinates_3D
	V(ACTIONet)$x = coor[, 1]
	V(ACTIONet)$y = coor[, 2]
	V(ACTIONet)$x3D = coor3D[, 1]
	V(ACTIONet)$y3D = coor3D[, 2]
	V(ACTIONet)$z3D = coor3D[, 3]
	V(ACTIONet)$color = rgb(vis.out$colors)

	arch.Lab = t(reconstruct.out$C_stacked) %*% grDevices::convertColor(color= vis.out$colors, from = 'sRGB', to = 'Lab')
	arch.colors = rgb(grDevices::convertColor(color=arch.Lab, from = 'Lab', to = 'sRGB'))
	arch.coordinates = t(reconstruct.out$C_stacked) %*% vis.out$coordinates
	arch.coordinates_3D = t(reconstruct.out$C_stacked) %*% vis.out$coordinates_3D
	arch.vis.out = list(colors = arch.colors, coordinates = arch.coordinates, coordinates_3D = arch.coordinates_3D)

	ACTIONet.out = list(ACTION.out = ACTION.out, reconstruct.out = reconstruct.out, build.out = build.out, vis.out = vis.out, ACTIONet = ACTIONet, arch.vis.out = arch.vis.out)

	# Add signature profile
	signature.profile = construct.signature.profile(sce = sce, ACTIONet.out = ACTIONet.out) 
	ACTIONet.out$signature.profile = signature.profile
		
	core.out = identify.core.archetypes(ACTIONet.out, core.z)
	H.core = runsimplexRegression(t(sce@reducedDims[["S_r"]]) %*% ACTIONet.out$reconstruct.out$C_stacked[, core.out$core.archs], t(sce@reducedDims[["S_r"]]))
	core.out$H = H.core	
	ACTIONet.out$core.out = core.out
	
	return(ACTIONet.out)
}

reconstruct.ACTIONet <- function(ACTIONet.out, sce, compactness_level = 50, thread_no = 8, epsilon = 3.0, LC = 1.0, auto_adjust_LC = FALSE, n_epochs = 500) {	  
	if( !("S_r" %in% names(sce@reducedDims)) ) {
		print("Please provide reduced sce with S_r slot")
		return(ACTIONet.out)
	}
	# Build ACTIONet
	build.out = buildAdaptiveACTIONet(H_stacked = ACTIONet.out$reconstruct.out$H_stacked, thread_no = thread_no, LC = LC, auto_adjust_LC = auto_adjust_LC, epsilon = epsilon)
	#build.out = ACTIONet.out$build.out
	
	# Layout ACTIONet
	vis.out = layoutACTIONet(build.out$ACTIONet, S_r = t(sce@reducedDims[["S_r"]]), compactness_level = compactness_level, n_epochs = n_epochs)

	# Construct igraph object
	ACTIONet = graph_from_adjacency_matrix(build.out$ACTIONet, mode="undirected", weighted = TRUE)
	coor = vis.out$coordinates;
	coor3D = vis.out$coordinates_3D
	V(ACTIONet)$x = coor[, 1]
	V(ACTIONet)$y = coor[, 2]
	V(ACTIONet)$x3D = coor3D[, 1]
	V(ACTIONet)$y3D = coor3D[, 2]
	V(ACTIONet)$z3D = coor3D[, 3]
	V(ACTIONet)$color = rgb(vis.out$colors)

	arch.Lab = t(ACTIONet.out$reconstruct.out$C_stacked) %*% grDevices::convertColor(color= vis.out$colors, from = 'sRGB', to = 'Lab')
	arch.colors = rgb(grDevices::convertColor(color=arch.Lab, from = 'Lab', to = 'sRGB'))
	arch.coordinates = t(ACTIONet.out$reconstruct.out$C_stacked) %*% vis.out$coordinates
	arch.coordinates_3D = t(ACTIONet.out$reconstruct.out$C_stacked) %*% vis.out$coordinates_3D
	arch.vis.out = list(colors = arch.colors, coordinates = arch.coordinates, coordinates_3D = arch.coordinates_3D)

	ACTIONet.out$build.out = build.out
	ACTIONet.out$vis.out = vis.out
	ACTIONet.out$ACTIONet = ACTIONet
	ACTIONet.out$arch.vis.out = arch.vis.out
	
	return(ACTIONet.out)
}

rerun.layout <- function(ACTIONet.out, sce, compactness = 50, init.slot = "S_r") {
	require(igraph)
	require(ACTIONet)

	# Layout ACTIONet
	vis.out = layoutACTIONet(ACTIONet.out$build.out$ACTIONet, S_r = t(sce@reducedDims[[init.slot]]), compactness_level = compactness, n_epochs = 500)
	ACTIONet.out$vis.out = vis.out
	
	# Construct igraph object
	ACTIONet = ACTIONet.out$ACTIONet
	coor = vis.out$coordinates;
	coor3D = vis.out$coordinates_3D
	V(ACTIONet)$x = coor[, 1]
	V(ACTIONet)$y = coor[, 2]
	V(ACTIONet)$x3D = coor3D[, 1]
	V(ACTIONet)$y3D = coor3D[, 2]
	V(ACTIONet)$z3D = coor3D[, 3]
	V(ACTIONet)$color = rgb(vis.out$colors)
	ACTIONet.out$ACTIONet = ACTIONet
	
	arch.Lab = t(ACTIONet.out$reconstruct.out$C_stacked) %*% grDevices::convertColor(color= vis.out$colors, from = 'sRGB', to = 'Lab')
	arch.colors = rgb(grDevices::convertColor(color=arch.Lab, from = 'Lab', to = 'sRGB'))
	arch.coordinates = t(ACTIONet.out$reconstruct.out$C_stacked) %*% vis.out$coordinates
	arch.coordinates_3D = t(ACTIONet.out$reconstruct.out$C_stacked) %*% vis.out$coordinates_3D
	arch.vis.out = list(colors = arch.colors, coordinates = arch.coordinates, coordinates_3D = arch.coordinates_3D)

	ACTIONet.out$arch.vis.out = arch.vis.out
	
	return(ACTIONet.out)
}

prune.ACTIONet <- function(ACTIONet.out, z.threshold = -1) {
	cn = coreness(ACTIONet.out$ACTIONet)
	z = scale(cn)
	filtered.cells = which(z < z.threshold)
	ACTIONet.out = remove.cells(ACTIONet.out, filtered.cells)
}

annotate.cells.using.markers <- function(ACTIONet.out, sce, marker.genes, alpha_val = 0.9, rand.sample.no = 100, thread_no = 8, imputation = "PageRank") {
	require(ACTIONet)
	require(igraph)
	require(Matrix)
	require(stringr)
	
	rownames(sce) = toupper(rownames(sce))
	
	GS.names = names(marker.genes)
	if(is.null(GS.names)) {
		GS.names = sapply(1:length(GS), function(i) sprintf('Celltype %s', i))
	}

	markers.table = do.call(rbind, lapply(names(marker.genes), function(celltype) {
		genes = marker.genes[[celltype]]
		if(length(genes) == 0)
			return(data.frame())
		
		is.signed = sum(sapply(genes, function(gene) {sgn_mark = stringr::str_sub(gene, start = -1); return(sgn_mark == "-" | sgn_mark == "+")}))
		if(! is.signed ) {
			df = data.frame(Gene = toupper(genes), Direction = +1, Celltype = celltype)
		} else {
			pos.genes = toupper(as.character(sapply(genes[grepl('+', genes, fixed = TRUE)], function(gene) stringr::str_replace(gene, stringr::fixed("+"), ""))))
			neg.genes = toupper(as.character(sapply(genes[grepl('-', genes, fixed = TRUE)], function(gene) stringr::str_replace(gene, stringr::fixed("-"), ""))))
			df = data.frame(Gene = c(pos.genes, neg.genes), Direction = c(rep(+1, length(pos.genes)), rep(-1, length(neg.genes))), Celltype= celltype)
		}
	}))
	markers.table = markers.table[markers.table$Gene %in% rownames(sce), ]
	if(dim(markers.table)[1] == 0) {
		print("No markers are left")
		return()
	}

	rows = match(markers.table$Gene, rownames(sce))
	if(imputation == "PageRank") { # PageRank-based imputation
		print("Using PageRank for imptation of marker genes")
		imputed.marker.expression = impute.genes.using.ACTIONet(ACTIONet.out, sce, markers.table$Gene, alpha_val, thread_no, prune = FALSE, rescale = FALSE)
	} else { # PCA-based imputation
		print("Using archImpute for imptation of marker genes")
		imputed.marker.expression = t(ACTIONet.out$signature.profile[rows, ACTIONet.out$core.out$core.archs] %*% ACTIONet.out$core.out$H)
	}
	colnames(imputed.marker.expression) = toupper(colnames(imputed.marker.expression))
	
	
	IDX = split(1:dim(markers.table)[1], markers.table$Celltype)		
	
	print("Computing significance scores")
	set.seed(0)
	Z = sapply(IDX, function(idx) {
		markers = toupper(as.character(markers.table$Gene[idx]))
		directions = markers.table$Direction[idx]		
		mask = markers %in% colnames(imputed.marker.expression)
		
		A = imputed.marker.expression[, markers[mask]]
		sgn = as.numeric(directions[mask])	
		stat = A %*% sgn		
		
		rand.stats = sapply(1:rand.sample.no, function(i) {
			rand.samples = sample.int(dim(imputed.marker.expression)[2], sum(mask))
			rand.A = imputed.marker.expression[, rand.samples]
			rand.stat = rand.A %*% sgn
		})
		
		cell.zscores = as.numeric((stat - apply(rand.stats, 1, mean) ) / apply(rand.stats, 1, sd))
		
		return(cell.zscores)
	})

	Z[is.na(Z)] = 0
	Labels = colnames(Z)[apply(Z, 1, which.max)]

	L = names(marker.genes)
	L = L[L %in% Labels]
	Labels = factor(Labels, levels = L)
	Labels.conf = apply(Z, 1, max)

	out.list = list(Labels = Labels, Labels.confidence = Labels.conf, cell2celltype.mat = Z)

	return(out.list)
}

remove.cells <- function(ACTIONet.out, filtered.cells, force = FALSE) {	
	core.cells = which(rowSums(ACTIONet.out$reconstruct.out$C_stacked) > 0)
	if(!force)
		selected.cells = sort(unique(union(core.cells, setdiff(1:dim(ACTIONet.out$reconstruct.out$C_stacked)[1], filtered.cells))))
	else
		selected.cells = sort(unique(core.cells))

	ACTIONet.out.pruned = ACTIONet.out

	C_trace = lapply(ACTIONet.out.pruned$ACTION.out$C, function(C) C[selected.cells, ])
	ACTIONet.out.pruned$ACTION.out$C = C_trace

	H_trace = lapply(ACTIONet.out.pruned$ACTION.out$H, function(H) H[, selected.cells])
	ACTIONet.out.pruned$ACTION.out$H = H_trace


	C_stacked.pruned = ACTIONet.out.pruned$reconstruct.out$C_stacked[selected.cells, ]
	ACTIONet.out.pruned$reconstruct.out$C_stacked = C_stacked.pruned

	H_stacked.pruned = ACTIONet.out.pruned$reconstruct.out$H_stacked[, selected.cells]
	ACTIONet.out.pruned$reconstruct.out$H_stacked = H_stacked.pruned


	ACTIONet.pruned = induced.subgraph(ACTIONet.out.pruned$ACTIONet, V(ACTIONet.out.pruned$ACTIONet)[selected.cells])
	ACTIONet.out.pruned$ACTIONet = ACTIONet.pruned


	ACTIONet.out.pruned$vis.out$coordinates = ACTIONet.out.pruned$vis.out$coordinates[selected.cells, ]
	ACTIONet.out.pruned$vis.out$coordinates_3D = ACTIONet.out.pruned$vis.out$coordinates_3D[selected.cells, ]
	ACTIONet.out.pruned$vis.out$colors = ACTIONet.out.pruned$vis.out$colors[selected.cells, ]


	ACTIONet.out.pruned$build.out$ACTIONet = ACTIONet.out.pruned$build.out$ACTIONet[selected.cells, selected.cells]
	ACTIONet.out.pruned$build.out$ACTIONet_asym = ACTIONet.out.pruned$build.out$ACTIONet_asym[selected.cells, selected.cells]

	ACTIONet.out.pruned$selected.cells = selected.cells
	
	return(ACTIONet.out.pruned)
}



cluster.ACTIONet <- function(ACTIONet.out, resolution_parameter = 1.0, arch.init =TRUE, thread_no = 8) {
	if(arch.init == TRUE) {
		U = ACTIONet.out$reconstruct.out$C_stacked[, ACTIONet.out$core.out$core.archs]
		U.smoothed = batchPR(ACTIONet.out$build.out$ACTIONet, U, thread_no = thread_no)
		
		initial.clusters = apply(U.smoothed, 1, which.max)
		
		initial.clusters.updated = cluster_label_prop(ACTIONet.out$ACTIONet, initial = initial.clusters)$membership
		initial.clusters.updated = match(initial.clusters.updated, sort(unique(initial.clusters.updated))) - 1
		
		clusters = unsigned_cluster(ACTIONet.out$build.out$ACTIONet, resolution_parameter, 0, initial.clusters.updated)
	} else {
		clusters = unsigned_cluster(ACTIONet.out$build.out$ACTIONet, resolution_parameter, 0)
	}
	
	counts = table(clusters)
	clusters[clusters %in% as.numeric(names(counts)[counts < 10])] = NA
	clusters = as.numeric(infer.missing.Labels(ACTIONet.out, clusters))
	
	clusters = factor(clusters, as.character(sort(unique(clusters))))

	return(clusters)
}

construct.sparse.backbone <- function(ACTIONet.out, reduction.slot = "S_r", stretch.factor = 10) {
	core.index = ACTIONet.out$core.out$core.archs 

	# Construct core-backbone
	#X = ACTIONet.out$reconstruct.out$C_stacked[, core.index]
	X = t(ACTIONet.out$reconstruct.out$H_stacked[core.index, ])
	
	mat = t(reducedDims(sce)[[reduction.slot]]) %*% X
	
	#mat = t(ACTIONet.out_ACTION$reconstruct.out$H_stacked[ACTIONet.out$core.out$core.archs, ])
	#mat = orthoProject(mat, Matrix::rowMeans(mat))
	
	backbone = cor(mat)
		
	backbone[backbone < 0] = 0
	diag(backbone) = 0

	backbone.graph = graph_from_adjacency_matrix(backbone, mode = "undirected", weighted = TRUE)

	# Construct t-spanner
	t = (2*stretch.factor-1)

	d = 1 - E(backbone.graph)$weight
	EL = get.edgelist(backbone.graph, names = FALSE)
	perm = order(d, decreasing = FALSE)

	backbone.graph.sparse = delete.edges(backbone.graph, E(backbone.graph))
	for(i in 1:length(d)) {
		u = EL[perm[i], 1]
		v = EL[perm[i], 2]
		sp = distances(backbone.graph.sparse, v = u, to = v)[1, 1]

		if(sp > t*d[perm[i]]) {
			backbone.graph.sparse = add.edges(backbone.graph.sparse, EL[perm[i], ], attr = list(weight = 1 - d[perm[i]]))
		}
	}

	backbone.graph = backbone.graph.sparse
	
	return(backbone.graph)
}
