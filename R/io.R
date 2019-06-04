import.sce.from.table <- function(fname, sep = '\t', header = TRUE, prefilter = FALSE, min.cells.per.gene = 10, min.genes.per.cell = 300) {
	require(Matrix)
	require(SingleCellExperiment)

	counts = read.table(fname, header = TRUE, sep = sep, as.is = TRUE)

	if(!is.numeric(counts[1, 1])) {
	  row.names = counts[, 1]
	  counts = counts[, -1]
	  rownames(counts) = row.names
	}

	counts = as(as.matrix(counts), 'sparseMatrix')
	
	sce <- SingleCellExperiment(
	  assays = list(counts = counts)
	) 		

	if(prefilter) {
		cell.counts = Matrix::rowSums(sce@assays[["counts"]] > 0)
		sce = sce[cell.counts > min.cells.per.gene, ]
		
		feature.counts = Matrix::colSums(sce@assays[["counts"]] > 0)
		sce = sce[, feature.counts > min.genes.per.cell]
	}	
	
	return(sce)
}

import.sce.from.10X <- function(input_path, mtx_file = 'matrix.mtx.gz', gene_annotations = 'features.tsv.gz', sample_annotations = 'barcodes.tsv.gz', sep = '\t', header = FALSE, prefilter = FALSE, min.cells.per.gene = 10, min.genes.per.cell = 300) {
	require(Matrix)
	require(scran)
	
	counts = readMM(paste(input_path, mtx_file, sep='/'));
	gene.table = read.table(paste(input_path, gene_annotations, sep='/'), header = header, sep=sep, as.is = TRUE)
	sample_annotations = read.table(paste(input_path, sample_annotations, sep='/'), header = header, sep=sep, as.is = TRUE)
	
	rownames(counts) = gene.table[, 1]
	colnames(counts) = sample_annotations[, 1]

	sce <- SingleCellExperiment(
	  assays = list(counts = counts), 
	  colData = sample_annotations,
	  rowData = gene.table
	) 		
	
	if(prefilter) {
		cell.counts = Matrix::rowSums(sce@assays[["counts"]] > 0)
		sce = sce[cell.counts > min.cells.per.gene, ]
		
		feature.counts = Matrix::colSums(sce@assays[["counts"]] > 0)
		sce = sce[, feature.counts > min.genes.per.cell]
	}
	
	return(sce)
}

# Uses loomR: https://satijalab.org/loomR/loomR_tutorial.html
import.sce.from.loom <- function(input_file, prefilter = FALSE, min.cells.per.gene = 10, min.genes.per.cell = 300) {
	require(loomR)
	require(scran)
	
	lfile <- loomR::connect(filename = input_file, mode = "r+")

	counts <- lfile$get.sparse(dataset = 'matrix')

	gene_annotations <- lfile$get.attribute.df(MARGIN = 1)

	meta.data <- setdiff(names(x = lfile[['col_attrs']]), c("ClusterID", "ClusterName", "Clusterings"))
	sample_annotations <- lfile$get.attribute.df(MARGIN = 2, attributes = meta.data)

	rownames(counts) = rownames(gene_annotations)
	colnames(counts) = colnames(sample_annotations)
		
	meta.df <- data.frame('NA' = rep.int(x = NA, times = lfile$shape[2]))
	rownames(x = meta.df) <- colnames(x = counts)
	for (i in meta.data) {
		X = lfile[[sprintf("col_attrs/%s", i)]][]
		if(is.data.frame(X))
			next
		meta.df[, i] <- X
	}
	meta.df <- meta.df[, -1]

	sce <- SingleCellExperiment(
	  assays = list(counts = counts), 
	  colData = meta.df
	) 		
	
	if(prefilter) {
		cell.counts = Matrix::rowSums(sce@assays[["counts"]] > 0)
		sce = sce[cell.counts > min.cells.per.gene, ]
		
		feature.counts = Matrix::colSums(sce@assays[["counts"]] > 0)
		sce = sce[, feature.counts > min.genes.per.cell]
	}	
	
	lfile$close_all()

	return(sce)
}


# https://satijalab.org/seurat/v3.0/conversion_vignette.html
import.sce.from.Seurat <- function(Seurat.obj) {
	require(Seurat)
	
	sce <- as.SingleCellExperiment(Seurat.obj)
	
	return(sce)
}

add.HTO.to.sce <- function(sce, fname, sep = '\t', summary.lines = 1) {
	HTO = read.csv(fname, sep = sep)
	
	if(!is.numeric(HTO[1, 1])) {
		rnames = HTO[, 1]
		HTO = HTO[, -1]
		rownames(HTO) = rnames
	}
	HTO = HTO[1:(dim(HTO)[1]-summary.lines), ]
	
	idx = match(colnames(sce), colnames(HTO))
	
	sce = sce[, !is.na(idx)]
	sub.HTO = HTO[, idx[!is.na(idx)]]

	sce@reducedDims[["HTO"]] = t(sub.HTO)
	
	return(sce)
}

add.ADT.to.sce <- function(sce, fname, sep = '\t', summary.lines = 1) {
	ADT = read.csv(fname, sep = sep)
	
	if(!is.numeric(ADT[1, 1])) {
		rnames = ADT[, 1]
		ADT = ADT[, -1]
		rownames(ADT) = rnames
	}
	ADT = ADT[1:(dim(ADT)[1]-summary.lines), ]
	
	idx = match(colnames(sce), colnames(ADT))
	
	sce = sce[, !is.na(idx)]
	sub.ADT = ADT[, idx[!is.na(idx)]]

	sce@reducedDims[["ADT"]] = t(sub.ADT)
	
	return(sce)
}

split.CITESeq.sce <- function(sce) {
  annot = rowData(sce)
  expression.rows = (annot$V3 == 'Gene Expression')
  sce.expression = sce[expression.rows, ]
  
  Ab.mat = as.matrix(sce@assays[["counts"]][!expression.rows, ])
  sce.expression@reducedDims[["ADT"]] = t(Ab.mat)
  rownames(sce.expression) = annot$V2[expression.rows]	
  
  return(sce.expression)
}

rename.sce.rows <- function(sce, from = "ENSEMBL", to = "SYMBOL") {
require(org.Hs.eg.db)

suppressWarnings(ids <- mapIds(org.Hs.eg.db,
	 keys=row.names(sce),
	 keytype=from,
	 column=to,
	 multiVals="first"))
	 
	 ids[is.na(ids)] = ''
	 
	 rownames(sce) = ids
	 
	 return(sce)
	 
	 
}
