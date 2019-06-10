
map.clusters <- function(Labels, clusters) {
  N = length(Labels)
  cluster.ids = sort(unique(clusters))
  if(is.factor(Labels))
	Label.ids = levels(Labels)
  else
	Label.ids = sort(unique(Labels))
  
  W = sapply(Label.ids, function(label) {
    idx1 = which(Labels == label)
    n1 = length(idx1)
    if(n1 == 0)
		return(array(0, length(cluster.ids)))
		
    log.pvals = sapply(cluster.ids, function(cluster.id) {
		idx2 = which(clusters == cluster.id)
		n2 = length(idx2)
		if(n2  == 0)
			return(0)
			
		success = intersect(idx1, idx2)
		
		pval = phyper(length(success)-1, n1, N - n1, n2, lower.tail = FALSE)
		return(-log10(pval))
    })
    return(log.pvals)
  })

  W[is.na(W)] = 0
  W[W > 300] = 300
  
  W.matched = MWM(W)
  W.matched = as(W.matched, 'dgTMatrix')
  
  updated.Labels = rep(NA, length(Labels))
  matched.clusters = W.matched@i+1
  matched.celltype = Label.ids[W.matched@j+1]
  for(k in 1:length(matched.clusters)) {
    updated.Labels[clusters == matched.clusters[k]] = matched.celltype[[k]]    
  }
  
  return(updated.Labels)
}

impute.genes.using.ACTIONet <- function(ACTIONet.out, sce, genes, alpha_val = 0.9, thread_no = 8, prune = FALSE, rescale = TRUE, expr.slot = "logcounts") {
	
	if(is.igraph(ACTIONet.out))
		ACTIONet = ACTIONet.out
	else
		ACTIONet = ACTIONet.out$ACTIONet
		
	if(length(V(ACTIONet)) != dim(sce)[2]) {
		R.utils::printf("Number of cells in the input sce (%d) doesn't match the number of vertices in the ACTIONet (%d)\n", dim(sce)[2], length(V(ACTIONet)))
	}
	G = as(get.adjacency(ACTIONet, attr = "weight"), 'dgTMatrix')


	matched.genes = intersect(genes, rownames(sce))	
	matched.idx = match(matched.genes, rownames(sce))

	# Smooth/impute gene expressions
	if(! (expr.slot %in% names(sce@assays)) ) {
		R.utils::printf('%s is not in assays of sce\n', expr.slot)
	}
	raw.gene.expression = Matrix::t(as(sce@assays[[expr.slot]][matched.idx, ], 'dgTMatrix'))
	U = raw.gene.expression
	U[U < 0] = 0
	cs = Matrix::colSums(U)
	cs[cs == 0] = 1
	U = Matrix::sparseMatrix(i = U@i+1, j = U@j+1, x = U@x / cs[U@j+1], dims = dim(U))

	imputed.gene.expression = batchPR(G, as.matrix(U), alpha_val, thread_no)

	# Prune values
	if(prune == TRUE) {
		imputed.gene.expression = apply(imputed.gene.expression, 2, function(x) {
			cond = sweepcut(ACTIONet.out$build.out$ACTIONet, x)
			idx = which.min(cond)

			perm = order(x, decreasing = TRUE)
			x[perm[(idx+1):length(x)]] = 0
			
			return(x)
		})
	}
	
	
	# rescale 
	if(rescale) {
		imputed.gene.expression = sapply(1:dim(raw.gene.expression)[2], function(col) {
			x = raw.gene.expression[, col]
			y = imputed.gene.expression[, col]
			
			x.Q = quantile(x, 0.99)
			y.Q = quantile(y, 0.99)

			if(y.Q == 0) {
				return(array(0, length(x)))
			}
			
			y = y * x.Q / y.Q;
			
			y[y > max(x)] = max(x);
			
			return(y)
		})
	}
	
	colnames(imputed.gene.expression) = matched.genes

	# Prune unobserved genes
	cs = Matrix::colSums(imputed.gene.expression)
	filtered.genes = which(cs == 0)
	if(length(filtered.genes) > 0)
		imputed.gene.expression = imputed.gene.expression[, -filtered.genes]
		
	return(imputed.gene.expression)
}


smooth.archetype.footprint <- function(ACTIONet.out, alpha_val = 0.9, thread_no = 8) {
  H.norm = as(t(ACTIONet.out$reconstruct.out$H_stacked), 'dgTMatrix')
  cs = Matrix::colSums(H.norm)
  cs[cs == 0] = 1
  H.norm = as.matrix(sparseMatrix(i = H.norm@i+1, j = H.norm@j+1, x = H.norm@x / cs[H.norm@j+1], dims = dim(H.norm)))
  
  pr.scores = batchPR(G = ACTIONet.out$build.out$ACTIONet, H.norm, alpha = alpha_val, thread_no = thread_no)
  
  arch.projections = apply(pr.scores, 2, function(x) {
  	cond = sweepcut(ACTIONet.out$build.out$ACTIONet, x)
  	idx = which.min(cond)
  
  	perm = order(x, decreasing = TRUE)
  	x[perm[(idx+1):length(x)]] = 0

  	return(x)
  })

  return(arch.projections)
}

impute.geneset.activity <- function(ACTIONet.out, sce, genes, alpha_val = 0.9, thread_no = 8) {
	imputed.gene.expression = impute.genes.using.ACTIONet(ACTIONet.out, sce, genes, alpha_val, thread_no, prune = FALSE, rescale = FALSE)
	

	gs.score = Matrix::rowMeans(imputed.gene.expression)
	
	return(gs.score)		
}



cluster.ACTIONet <- function(ACTIONet.out, resolution = 1.0, method = "ModularityVertexPartition") {
	library(leiden)
	
	clusters <- leiden(ACTIONet.out$build.out$ACTIONet, resolution_parameter = resolution, partition_type = method)
	
	return(clusters)	
}



infer.missing.Labels <- function(ACTIONet.out, Labels) {	
	require(igraph)
	
	if(is.igraph(ACTIONet.out))
		ACTIONet = ACTIONet.out
	else
		ACTIONet = ACTIONet.out$ACTIONet
		
	G = as(get.adjacency(ACTIONet), 'dgTMatrix')
	cs = Matrix::colSums(G)
	cs[cs == 0] = 1
	P = Matrix::sparseMatrix(i = G@i+1, j = G@j+1, x = G@x / cs[G@i+1], dims = dim(G))

	
	na.mask = is.na(Labels)
	if(sum(na.mask) > 0) {
		ct = sort(unique(Labels[!na.mask]))
		
		celltype.mask = sapply(ct, function(l) as.numeric(Labels == l))
		colnames(celltype.mask) = ct
		celltype.mask[is.na(celltype.mask)] = 0
		
		Neighbors.vote = as.matrix(P %*% celltype.mask)		
		Neighbors.Labels = ct[apply(Neighbors.vote, 1, which.max)]
		
		Labels[na.mask] = Neighbors.Labels[na.mask]
	}	
	
	return(Labels)
}


prop.Labels <- function(ACTIONet.out, Labels, max_it = 3) {
	require(Matrix)
	require(igraph)
	
	G = as(get.adjacency(ACTIONet.out$ACTIONet, attr = "weight"), 'dgTMatrix')
	cs = Matrix::colSums(G)
	cs[cs == 0] = 1
	P = Matrix::sparseMatrix(i = G@i+1, j = G@j+1, x = G@x / cs[G@i+1], dims = dim(G))

	updated_labels = Labels
	if(is.factor(updated_labels)) {
		Ll = levels(updated_labels)
	} else {
		Ll = sort(unique(updated_labels))
	}
	
	for(i in 1:max_it) {
		celltype.mask = sapply(Ll, function(l) as.numeric(updated_labels == l))
		X = celltype.mask
#		cs = Matrix::colSums(celltype.mask)
#		cs[cs == 0] = 1
#		X = scale(celltype.mask, center = FALSE, scale = cs)


		Neighbors.vote = as.matrix(P %*% X)		
		Neighbors_labels = Ll[apply(Neighbors.vote, 1, which.max)]
		idx = which(Neighbors_labels != updated_labels)
		updated_labels[idx] = Neighbors_labels[idx]
	}
	
	return(updated_labels)
}

update.Labels <- function(ACTIONet.out, Labels) {
	
	ACTIONet = ACTIONet.out$ACTIONet
	if(is.character(Labels))
		Labels = factor(Labels, levels = sort(unique(Labels)))
	
	require(igraph)
	require(Matrix)
	
	cl.out = cluster_label_prop(ACTIONet, initial = as.numeric(Labels))
	clusters = cl.out$membership

	updated.Labels = map.clusters(Labels, clusters)

	while(sum(is.na(updated.Labels)) > 0)
		updated.Labels = infer.missing.Labels(ACTIONet, updated.Labels)
	
	
	updated.Labels = factor(updated.Labels, levels = levels(Labels))

	return(updated.Labels) 
}

annotate.archetype <- function(ACTIONet.out, Labels, rand_perm_no = 1000) {
	require(ACTIONet)
	require(SingleCellExperiment)
	
	if(!is.factor(Labels)) {
		Labels = factor(Labels, levels = sort(unique(Labels)))
	}
	
	Labels.mask = sapply(levels(Labels), function(l) as.numeric(Labels == l))
	Enrichment.Z = phenotypeEnrichment(ACTIONet.out$reconstruct.out$H_stacked, phenotype_associations = Labels.mask, rand_perm_no = 1000)
	colnames(Enrichment.Z) = levels(Labels)

	archetypeLabels = levels(Labels)[apply(Enrichment.Z, 1, which.max)]
	archetypeLabels = factor(archetypeLabels, levels = levels(Labels))
	
	out.list = list(archetypeLabels = archetypeLabels, Enrichment = Enrichment.Z)
	
	return(out.list)
}


cluster.ACTIONet <- function(ACTIONet.out, resolution_parameter = 1.0) {
	clusters = unsigned_cluster(ACTIONet.out$build.out$ACTIONet, resolution_parameter)
	
	return(clusters)
}
