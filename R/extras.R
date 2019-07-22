
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
	
	if('clusters' %in% names(ACTIONet.out)) {
		imputed.gene.expression = zoned_diffusion(ACTIONet.out$build.out$ACTIONet, clusters, U, alpha = alpha_val)
		
	} else {
		imputed.gene.expression = batchPR(G, U, alpha_val, thread_no)
	}
	
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
			
			x.Q = quantile(x, 1)
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


infer.missing.Labels <- function(ACTIONet.out, Labels, double.stochastic = FALSE) {	
	if(is.igraph(ACTIONet.out))
		ACTIONet = ACTIONet.out
	else
		ACTIONet = ACTIONet.out$ACTIONet
	
	
  if(!is.factor(Labels)) {
    Labels = factor(Labels, levels = sort(unique(Labels)))
  }
  
  A = as(get.adjacency(ACTIONet, attr = "weight"), 'dgTMatrix')
  eps = 1e-16
  rs = Matrix::rowSums(A)
  P = sparseMatrix(i = A@i+1, j = A@j+1, x = A@x/rs[A@i+1], dims = dim(A))  
  if(double.stochastic == TRUE) {
    w = sqrt(Matrix::colSums(P)+eps)
    W = P %*% Matrix::Diagonal(x = w, n = length(w))
    P = W %*% Matrix::t(W)
  } 
  
  
  na.mask = is.na(Labels)
  
  if (sum(na.mask) > 0) {
	counts = table(Labels)
	p = counts/sum(counts)
	X = sapply(names(p), function(celltype) {
	x = as.numeric(Matrix::sparseVector(x = 1, i = which(Labels == celltype), length = length(Labels)))
	})

	Exp = array(1, nrow(A)) %*% t(p)
	Obs = P %*% X 

	Lambda = Obs - Exp


	w2 = Matrix::rowSums(P^2)
	Nu = w2 %*% t(p)

	a = as.numeric(qlcMatrix::rowMax(P)) %*% t(array(1, length(p)))


	logPval = (Lambda^2) / (2 * (Nu + (a*Lambda)/3))
	logPval[Lambda < 0] = 0
	logPval[is.na(logPval)] = 0
	updated.Labels = factor(levels(Labels)[apply(logPval, 1, which.max)], levels = levels(Labels))

	Labels[na.mask] = updated.Labels[na.mask]
  }
  
  return(Labels)  
}


update.Labels <- function(ACTIONet.out, Labels, max.iter = 3, double.stochastic = FALSE) {	
	if(is.igraph(ACTIONet.out))
		ACTIONet = ACTIONet.out
	else
		ACTIONet = ACTIONet.out$ACTIONet
	
	
  if(!is.factor(Labels)) {
    Labels = factor(Labels, levels = sort(unique(Labels)))
  }
  
  A = as(get.adjacency(ACTIONet, attr = "weight"), 'dgTMatrix')
  eps = 1e-16
  rs = Matrix::rowSums(A)
  P = sparseMatrix(i = A@i+1, j = A@j+1, x = A@x/rs[A@i+1], dims = dim(A))  
  if(double.stochastic == TRUE) {
    w = sqrt(Matrix::colSums(P)+eps)
    W = P %*% Matrix::Diagonal(x = w, n = length(w))
    P = W %*% Matrix::t(W)
  } 
  

	for(it in 1:max.iter) {
		R.utils::printf("iter %d\n", it)
		counts = table(Labels)
		p = counts/sum(counts)
		X = sapply(names(p), function(celltype) {
			x = as.numeric(Matrix::sparseVector(x = 1, i = which(Labels == celltype), length = length(Labels)))
		})

		Exp = array(1, nrow(A)) %*% t(p)
		Obs = P %*% X 

		Lambda = Obs - Exp


		w2 = Matrix::rowSums(P^2)
		Nu = w2 %*% t(p)

		a = as.numeric(qlcMatrix::rowMax(P)) %*% t(array(1, length(p)))


		logPval = (Lambda^2) / (2 * (Nu + (a*Lambda)/3))
		logPval[Lambda < 0] = 0
		logPval[is.na(logPval)] = 0
		Labels = factor(levels(Labels)[apply(logPval, 1, which.max)], levels = levels(Labels))
	}
  
  return(Labels)  
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
