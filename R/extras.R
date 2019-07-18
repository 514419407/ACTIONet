
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


impute.genes.using.ACTIONet <- function(ACTIONet.out, sce, genes, alpha_val = 0.9, thread_no = 8, prune = FALSE, rescale = FALSE, expr.slot = "logcounts") {
	require(igraph)
	
	genes = unique(genes)
	
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
	
	if(length(matched.idx) > 1) {
		raw.gene.expression = Matrix::t(as(sce@assays[[expr.slot]][matched.idx, ], 'dgTMatrix'))
		U = raw.gene.expression
		U[U < 0] = 0
		cs = Matrix::colSums(U)
		U = as.matrix(Matrix::sparseMatrix(i = U@i+1, j = U@j+1, x = U@x / cs[U@j+1], dims = dim(U)))
		U = U[, cs > 0]
		gg = matched.genes[cs > 0]
	}
	else {		
		raw.gene.expression = matrix(sce@assays[[expr.slot]][matched.idx, ])
		U = raw.gene.expression / sum(raw.gene.expression)
		gg = matched.genes
	}	
	
	#clusters = cluster.ACTIONet(ACTIONet.out)
	#imputed.gene.expression = zoned_diffusion(ACTIONet.out$build.out$ACTIONet, clusters, U, alpha = 0.95)
	
	imputed.gene.expression = batchPR(G, U, alpha_val, thread_no)
	
	imputed.gene.expression[is.na(imputed.gene.expression)] = 0
	
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
			
			x.Q = max(x, 1)
			y.Q = quantile(y, 1)

			if(y.Q == 0) {
				return(array(0, length(x)))
			}
			
			y = y * x.Q / y.Q;
			
			y[y > max(x)] = max(x);
			
			return(y)
		})
	}
	
	colnames(imputed.gene.expression) = gg
		
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
	set.seed(0)	
	
	if(is.igraph(ACTIONet.out))
		ACTIONet = ACTIONet.out
	else
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
	
	
	# if(is.igraph(ACTIONet.out))
	# 	ACTIONet = ACTIONet.out
	# else
	# 	ACTIONet = ACTIONet.out$ACTIONet
	# 
	# if(is.character(Labels)) {
	# 	initial.clusters = match(Labels, sort(unique(Labels))) - 1
	# } else if(is.factor(Labels)) {
	# 	initial.clusters = as.numeric(Labels)-1
	# } else if(is.numeric(Labels)) {
	# 	initial.clusters = Labels - min(Labels)
	# }
	# 
	# G = get.adjacency(ACTIONet, attr = "weight")
	# clusters = unsigned_cluster(G, resolution_parameter, seed, initial.clusters)
	# 
	# updated.Labels = map.clusters(Labels, clusters)
	# 
	# while(sum(is.na(updated.Labels)) > 0)
	# 	updated.Labels = infer.missing.Labels(ACTIONet, updated.Labels)
	# 
	# 
	# updated.Labels = factor(updated.Labels, levels = levels(Labels))
	# 
	# return(updated.Labels) 
	
	
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

orthoProject <- function(A, S) {
  A = scale(A)
  S = scale(S)
  A_r = A - S %*% MASS::ginv(t(S) %*% S) %*% (t(S)%*%A)
  A_r = scale(A_r)
  return(A_r)
}
