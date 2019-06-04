reduce.sce <- function(sce, reduced_dim = 50, max.iter = 5, normalize = TRUE) {
	require(scran)
	require(scater)
	require(ACTIONet)

	if(normalize) {
		sce.norm = sce    
		A = as(sce@assays[["counts"]], 'dgTMatrix')
		cs = Matrix::colSums(A)    
		B = Matrix::sparseMatrix(i = A@i+1, j = A@j+1, x = log(1 + median(cs)*(A@x / cs[A@j + 1])), dims = dim(A))
		sce.norm@assays[["logcounts"]] = B   
	} else {
		sce.norm = sce
	}
	
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


reduce.and.batch.correct.sce <- function(sce, batch.attr=NA, reduced_dim = 50, k=20) {
	require(scran)
	require(scater)
	require(ACTIONet)

	if(is.na(batch.attr)) {
		print("You need to provide the batch vector/attr");
		return()
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

identify.core.archetypes <- function(ACTIONet.out, pruning.zscore.threshold = 0) {
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


run.ACTIONet <- function(sce, k_max = 20, compactness_level = 50, thread_no = 8, LC = 1.0, arch.specificity.z = -1, core.z = 1, sce.data.attr = "logcounts") {
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
	build.out = buildAdaptiveACTIONet(H_stacked = reconstruct.out$H_stacked, thread_no = 8, LC = 1.0)

	# Layout ACTIONet
	vis.out = layoutACTIONet(build.out$ACTIONet, S_r = t(sce@reducedDims[["S_r"]]), compactness_level = 50, n_epochs = 500)

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
		
	core.out = identify.core.archetypes(ACTIONet.out, 1)
	H.core = runsimplexRegression(t(sce@reducedDims[["S_r"]]) %*% ACTIONet.out$reconstruct.out$C_stacked[, core.out$core.archs], t(sce@reducedDims[["S_r"]]))
	core.out$H = H.core	
	ACTIONet.out$core.out = core.out
	
	return(ACTIONet.out)
}

annotate.cells.using.markers <- function(ACTIONet.out, sce, marker.genes, alpha_val = 0.9, rand.sample.no = 100, thread_no = 8) {
	require(ACTIONet)
	require(igraph)
	require(Matrix)

	
	GS.names = names(marker.genes)
	if(is.null(GS.names)) {
		GS.names = sapply(1:length(GS), function(i) sprintf('Celltype %s', i))
	}

	markers.table = do.call(rbind, lapply(names(marker.genes), function(celltype) {
		genes = marker.genes[[celltype]]
		if(length(genes) == 0)
			return(data.frame())
			
		if(sum(grepl('-|+', genes, fixed = TRUE)) == 0) {
			df = data.frame(Gene = genes, Direction = +1, Celltype = celltype)
		} else {
			pos.genes = as.character(sapply(genes[grepl('+', genes, fixed = TRUE)], function(gene) stringr::str_replace(gene, stringr::fixed("+"), "")))
			neg.genes = as.character(sapply(genes[grepl('-', genes, fixed = TRUE)], function(gene) stringr::str_replace(gene, stringr::fixed("-"), "")))
			df = data.frame(Gene = c(pos.genes, neg.genes), Direction = c(rep(+1, length(pos.genes)), rep(-1, length(neg.genes))), Celltype= celltype)
		}
	}))
	markers.table = markers.table[markers.table$Gene %in% rownames(sce), ]
	if(dim(markers.table)[1] == 0) {
		print("No markers are left")
		return()
	}

	imputed.marker.expression = impute.genes.using.ACTIONet(ACTIONet.out, sce, markers.table$Gene, alpha_val, thread_no, prune = TRUE)
	
	
	
	IDX = split(1:dim(markers.table)[1], markers.table$Celltype)		
	
	set.seed(0)
	Z = sapply(IDX, function(idx) {
		markers = as.character(markers.table$Gene[idx])
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

	
	Labels = colnames(Z)[apply(Z, 1, which.max)]

	L = names(marker.genes)
	L = L[L %in% Labels]
	Labels = factor(Labels, levels = L)
	Labels.conf = apply(Z, 1, max)

	out.list = list(Labels = Labels, Labels.confidence = Labels.conf, cell2celltype.mat = Z)

	return(out.list)
}

remove.cells <- function(ACTIONet.out, filtered.cells) {	
	core.cells = which(rowSums(ACTIONet.out$reconstruct.out$C_stacked) > 0)
	selected.cells = sort(unique(union(core.cells, setdiff(1:dim(ACTIONet.out$reconstruct.out$C_stacked)[1], filtered.cells))))

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

