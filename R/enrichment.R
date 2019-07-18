ACTIONet.autocorrelation <- function(ACTIONet.out, variable) {
	
}

ACTIONet.categorical.autocorrelation <- function(ACTIONet.out, labels) {
  A = as(ACTIONet.out$build.out$ACTIONet, 'dgTMatrix')
  
  n = length(labels)
  counts = table(labels)
  categoties = as.numeric(names(counts))
  
  pvec = as.vector(counts) / n
  
  k = length(pvec)
  
  w = A@x
  s0 = sum(w)
  s1 = sum(4*w^2) / 2
  s2 = sum( (colSums(A) + rowSums(A))^2 )
  
  m1.rawphi =  (s0/(n*(n-1)))*(n^2*k*(2-k) - n*sum(1/pvec))
  
  Q1 = sum(1/pvec)
  Q2 = sum(1/pvec^2)
  Q3 = sum(1/pvec^3)
  Q22 = sum((1/pvec)%*%t(1/pvec))
  E1 = (n^2*Q22 - n*Q3)/(n*(n-1))
  E2 = 4*n^3*Q1 - 4*n^3*k*Q1 + n^3*k^2*Q1 - 2*( 2*n^2*Q2 - n^2*k*Q2) + 2*n*Q3 - n^2*Q22
  E2 = E2/(n*(n-1)*(n-2))
  
  A1 = 4*n^4*k^2 - 4*n^4*k^3 + n^4*k^4 - (2*n^3*k*Q1 - n^3*k^2*Q1)
  A2 = 4*n^3*Q1 - 4*n^3*k*Q1 + n^3*k^2*Q1 - (2*n^2*Q2 - n^2*k*Q2)
  Apart = A1 - 2*A2
  
  B1 = 4*n^3*Q1 - 4*n^3*k*Q1 + n^3*k^2*Q1 - (2*n^2*Q2 - n^2*k*Q2)
  B2 = 2*n^2*Q2 - n^2*k*Q2 - n*Q3
  B3 = n^2*Q22 - n*Q3
  Bpart = B1 - B2 - B3
  
  C1 = 2*n^3*k*Q1 - n^3*k^2*Q1 - n^2*Q22
  C2 = 2*n^2*Q2 - n^2*k*Q2 - n*Q3
  Cpart = C1 - 2*C2
  
  E3 = (Apart - 2*Bpart - Cpart) / (n*(n-1)*(n-2)*(n-3))
  
  m2.rawphi = s1*E1 + (s2 - 2*s1)*E2 + (s0^2 - s2 + s1)*E3
  
  v_i = v[A@i+1];
  v_j = v[A@j+1];
  
  p_i = pvec[match(v_i, categoties)]
  p_j = pvec[match(v_j, categoties)]

  rawphi = sum(w*(2*(v_i == v_j)-1) / (p_i*p_j))
  
  mean.rawphi = m1.rawphi
  var.rawphi = m2.rawphi - mean.rawphi^2
  phi.z = (rawphi - mean.rawphi) / sqrt(var.rawphi)
  phi.logPval = -log10(pnorm(phi.z, lower.tail = FALSE))
    
  return(list(z = phi.z, logPval = phi.logPval, phi = rawphi))
}

# HGT tail bound
Kappa <- function(p, q) {
  kl = array(1, length(p))

  suppressWarnings( {a = p*log(p / q)} )
  a[p == 0] = 0

  suppressWarnings( {b = (1-p)*log((1-p)/(1-q))} )
  b[p == 1] = 0    

  k = a + b
  return(k)
}

HGT_tail <- function(population.size, success.count, sample.size, observed.success) {
  if(success.count == 0)
    return(1)

  success.rate = success.count / population.size
  expected.success = sample.size * success.rate
  delta = (observed.success/expected.success)-1

  log.tail_bound = sample.size*Kappa((1+delta)*success.rate, success.rate)
  log.tail_bound[delta < 0] = 0
  log.tail_bound[is.na(log.tail_bound)] = 0
  
  return(log.tail_bound)
}

geneset.enrichment.gProfiler <- function(genes, top.terms = 10, col = "tomato", organism = "hsapiens") {
	require(gProfileR)
	require(ggpubr)

	terms = gprofiler(genes, ordered_query = FALSE, hier_filtering = 'none', exclude_iea=FALSE, correction_method='fdr', src_filter = c('GO:REAC'), organism = organism)

	terms$logPval = -log10(terms$p.value)


	too.long = which(sapply(terms$term.name, function(x) stringr::str_length(x)) > 50)
	terms = terms[-too.long,]

	terms = terms[order(terms$logPval, decreasing = TRUE),]


	p = ggbarplot(terms[1:min(top.terms, sum(terms$logPval > 2)), ], x="term.name", y="logPval", sort.val = "asc", orientation = "horiz", fill=col, xlab = "", ylab="") + geom_hline(yintercept = -log10(0.05), col="gray", lty=2)

	plot(p)    

}

archetype.geneset.enrichment <- function(ACTIONet.out, genesets) {
	signature.profile = ACTIONet.out$signature.profile
	
	if(is.null(rownames(signature.profile))) {
		print("Rows of the signature profile have to be named with genes.")
	}

	require(Matrix)
	
	if(is.list(genesets)) {
		ind.mat = as(sapply(genesets, function(gs) as.numeric(rownames(signature.profile) %in% gs)), 'sparseMatrix')
		rownames(ind.mat) = rownames(signature.profile)
	}
	else {
		ind.mat = genesets
	}
	common.genes = intersect(rownames(signature.profile), rownames(ind.mat))
	
	idx = match(common.genes, rownames(signature.profile))
	signature.profile = signature.profile[idx, ]
	
	idx = match(common.genes, rownames(ind.mat))
	X = ind.mat[idx, ]
		

	# Normalize scores to avoid heavy-tail side-effect(s)
	pos.scores = signature.profile
	pos.scores[pos.scores < 0] = 0	
	A = pos.scores  / max(pos.scores)

	
	p_c = Matrix::colMeans(X)

	Obs = as.matrix(Matrix::t(X) %*% A)
	Exp = as.matrix(p_c %*% Matrix::t(Matrix::colSums(A)))
	Nu = as.matrix(p_c %*% Matrix::t(Matrix::colSums(A^2)))

	Lambda = Obs - Exp
	
	a = apply(A, 2, max)
	ones = array(1, dim = dim(Lambda)[1])
		
	logPvals = Lambda^2 / (2*(Nu + (ones%*% Matrix::t(a))*Lambda/3))
	
	rownames(logPvals) = colnames(ind.mat)
	colnames(logPvals) = colnames(signature.profile)

	logPvals[Lambda  < 0] = 0
	logPvals[is.na(logPvals)] = 0
	
	logPvals = logPvals / log(10)

	
	return(logPvals)
}

