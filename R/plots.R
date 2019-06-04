plot.ACTIONet <- function(ACTIONet.out, color.attr = NA, transparency.attr = NA, size.attr = NA, cex = 1, CPal = "Dark2", legend.pos = "bottomright") {
	require(ggplot2)
	require(ggpubr)
	
	if(is.igraph(ACTIONet.out))
		ACTIONet = ACTIONet.out
	else
		ACTIONet = ACTIONet.out$ACTIONet
		
	nV = length(V(ACTIONet))
	coor = cbind(V(ACTIONet)$x, V(ACTIONet)$y)


	if(is.numeric(color.attr)) {
		v = sort(unique(color.attr))
		Annot = as.character(v)
		Pal = ggpubr::get_palette(CPal, length(Annot))
		vCol = Pal[color.attr]
	}
	else if(is.factor(color.attr)) {
		Annot = levels(color.attr)
		Pal = ggpubr::get_palette(CPal, length(Annot))
		vCol = Pal[color.attr]
	} 
	else if(is.character(color.attr)) {
		color.attr = factor(color.attr, levels = sort(unique(color.attr)))
		Annot = levels(color.attr)
		Pal = ggpubr::get_palette(CPal, length(Annot))
		vCol = Pal[color.attr]
	}
	else {
		vCol = V(ACTIONet)$color
		Annot = NA
	}		

	HSV = rgb2hsv(col2rgb(vCol))
	HSV[3, ] = HSV[3, ]*0.7
	vCol.border = apply(HSV, 2, function(v) do.call(hsv, as.list(v)))

	if(is.numeric(transparency.attr)) {
		vCol = ggplot2::alpha(vCol, 1 / (1 + exp(-scale(transparency.attr))))
		vCol.border = ggplot2::alpha(vCol.border, 1 / (1 + exp(-scale(transparency.attr))))
	}


	if(is.numeric(size.attr)) {
		cex = cex * (1 / (1 + exp(-scale(size.attr))))
	}


	plot(coor, pch=21, bg=vCol, col=vCol.border, cex=cex, axes=FALSE, xlab="", ylab="");
	if(length(Annot) > 1) 
		legend(legend.pos, legend=Annot, col = Pal, bty = "n", pch=20 , pt.cex = 1.5, cex = 0.5, horiz = FALSE, inset = c(0.1, 0.1))

}


plot.ACTIONet.igraph <- function(ACTIONet.out, color.attr = NA, transparency.attr = NA, size.attr = NA, cex = 2, CPal = "d3", legend.pos = "bottomright") {
	require(ggplot2)
	require(ggpubr)
	
	if(is.igraph(ACTIONet.out))
		ACTIONet = ACTIONet.out
	else
		ACTIONet = ACTIONet.out$ACTIONet
		
	nV = length(V(ACTIONet))
	coor = cbind(V(ACTIONet)$x, V(ACTIONet)$y)


	if(is.numeric(color.attr)) {
		v = sort(unique(color.attr))
		Annot = as.character(v)
		Pal = ggpubr::get_palette(CPal, length(Annot))
		vCol = Pal[color.attr]
	}
	else if(is.factor(color.attr)) {
		Annot = levels(color.attr)
		Pal = ggpubr::get_palette(CPal, length(Annot))
		vCol = Pal[color.attr]
	} 
	else if(is.character(color.attr)) {
		color.attr = factor(color.attr, levels = sort(unique(color.attr)))
		Annot = levels(color.attr)
		Pal = ggpubr::get_palette(CPal, length(Annot))
		vCol = Pal[color.attr]
	}
	else {
		vCol = V(ACTIONet)$color
		Annot = NA
	}		


	HSV = rgb2hsv(col2rgb(vCol))
	HSV[3, ] = HSV[3, ]*0.7
	vCol.border = apply(HSV, 2, function(v) do.call(hsv, as.list(v)))

	if(is.numeric(transparency.attr)) {
		trans.factor = 1 / ( 1+ exp(-scale(transparency.attr)) )
		trans.factor = trans.factor ^ 2
		
		vCol = ggplot2::alpha(vCol, trans.factor)
		vCol.border = ggplot2::alpha(vCol.border, trans.factor)
	}

	if(is.numeric(size.attr)) {
		size.factor = 1 / ( 1+ exp(-scale(transparency.attr)) )
		size.factor = trans.factor ^ 2

		cex = cex * size.attr
		 
	}

	sketch.graph = ACTIONet
	sketch.graph = delete.edges(sketch.graph, E(sketch.graph))
	V(sketch.graph)$size = cex
	V(sketch.graph)$color = vCol
	V(sketch.graph)$frame.color = vCol.border
	

	plot(sketch.graph, vertex.label=NA, layout=coor)
	if(length(Annot) > 1) 
		legend(legend.pos, legend = Annot, fill=Pal, cex = 0.5)

}

plot.ACTIONet.3D <- function(ACTIONet.out, color.attr = NA, transparency.attr = NA, size.attr = NA, node.size = 0.2, CPal = "d3") {
	require(ggplot2)
	require(ggpubr)
	require(threejs)
	
	if(is.igraph(ACTIONet.out))
		ACTIONet = ACTIONet.out
	else
		ACTIONet = ACTIONet.out$ACTIONet
		
	nV = length(V(ACTIONet))
	coor = cbind(V(ACTIONet)$x3D, V(ACTIONet)$y3D, V(ACTIONet)$z3D)

	if(is.numeric(color.attr)) {
		v = sort(unique(color.attr))
		Annot = as.character(v)
		Pal = ggpubr::get_palette(CPal, length(Annot))
		vCol = Pal[color.attr]
	}
	else if(is.factor(color.attr)) {
		Annot = levels(color.attr)
		Pal = ggpubr::get_palette(CPal, length(Annot))
		vCol = Pal[color.attr]
	} 
	else if(is.character(color.attr)) {
		color.attr = factor(color.attr, levels = sort(unique(color.attr)))
		Annot = levels(color.attr)
		Pal = ggpubr::get_palette(CPal, length(Annot))
		vCol = Pal[color.attr]
	}
	else {
		vCol = V(ACTIONet)$color
		Annot = NA
	}

	HSV = rgb2hsv(col2rgb(vCol))
	HSV[3, ] = HSV[3, ]*0.7
	vCol.border = apply(HSV, 2, function(v) do.call(hsv, as.list(v)))

	if(is.numeric(transparency.attr)) {
		vCol = ggplot2::alpha(vCol, 1 / (1 + exp(-scale(transparency.attr))))
		vCol.border = ggplot2::alpha(vCol.border, 1 / (1 + exp(-scale(transparency.attr))))
	}

	if(is.numeric(size.attr)) {
		node.size = node.size * (1 / (1 + exp(-scale(size.attr))))
	}


	scatterplot3js(x = coor[, 1], y = coor[, 2], z = coor[, 3], axis.scales = FALSE, size = node.size, axis = F, grid = F, color = vCol, stroke = vCol.border, bg='black')
}


plot.ACTIONet.gradient <- function(ACTIONet.out, cell.scores, alpha_val = 0, transform = FALSE, subset.scores = NULL, prefilter = FALSE, node.size = 2, CPal = "inferno", title = "scores") {
	require(viridis)
	require(scales)

	sketch.graph = ACTIONet.out$ACTIONet
	sketch.graph = delete.edges(sketch.graph, E(sketch.graph))

	if( !is.null(subset.scores) ) {
		print("Subsetting cells");	
		x = cell.scores[subset.scores]		
		G = induced_subgraph(ACTIONet.out$ACTIONet, V(ACTIONet.out$ACTIONet)[subset.scores])
	}
	else {
		x = cell.scores
		G = ACTIONet.out$ACTIONet
	}

	#z = (x - median(x)) / mad(x)
	z = scale(x)
	# Adjust outliers
	z[z > 3] = 3
	z[z < -3] = -3
	
	if(transform) {
		print("Sigmoid transform scores")
		x = 1 / (1 + exp(-z))
	}
	else
		x = (z - min(z)) / (max(z) - min(z))
		
	
	if(alpha_val != 0)
		cell.scores = page_rank(G, personalized =  x, damping = alpha_val)$vector		
	else 
		cell.scores = x
	
	
	cell.scores.scaled = 0.01 + 0.99*cell.scores / max(cell.scores)

	Pal_grad = switch(CPal, 
	"inferno" = inferno(500, alpha = 0.5),
	"magma" = magma(500, alpha = 0.5),
	"viridis" = viridis(500, alpha = 0.5),
	"Reds" = ggpubr::get_palette("Reds", 500),
	"BlGrRd" = colorRampPalette(c('blue', 'lightgrey', 'red'))(500))
	#Pal_grad = c(rep(rgb(1, 1, 1, 0.5), 50), Pal_grad)
	


	if( !is.null(subset.scores) ) {
		vCol = rep('lightgrey', length(cell.scores))
		vCol[subset.scores] = scales::col_bin(Pal_grad, domain = NULL, bins = 100)(cell.scores.scaled)
	} else {
		vCol = scales::col_bin(Pal_grad, domain = NULL, bins = 100)(cell.scores.scaled)
	}
	
	
	HSV = rgb2hsv(col2rgb(vCol))
	HSV[3, ] = HSV[3, ]*0.8
	vCol.border = apply(HSV, 2, function(v) do.call(hsv, as.list(v)))	
	

	V(sketch.graph)$size = node.size
	coor = cbind(V(sketch.graph)$x, V(sketch.graph)$y)



	V(sketch.graph)$color = vCol
	V(sketch.graph)$frame.color =  vCol.border

	plot(sketch.graph, vertex.label=NA, layout=coor, main = title)
	
}





plot.ACTIONet.cell.states <- function(ACTIONet.out, cex = 2, min.cor = 0.3) {
	require(ggplot2)
	require(ggpubr)
	
	sketch.graph = ACTIONet.out$ACTIONet
	sketch.graph = delete.edges(sketch.graph, E(sketch.graph))
	
		
	V(sketch.graph)$color = 'lightgrey'
	V(sketch.graph)$size = cex

	attr = list(color = ACTIONet.out$arch.vis.out$colors[ACTIONet.out$core.out$core.archs], x = ACTIONet.out$arch.vis.out$coordinates[ACTIONet.out$core.out$core.archs, 1], y = ACTIONet.out$arch.vis.out$coordinates[ACTIONet.out$core.out$core.archs, 2], size = rep(cex*1.5, length(ACTIONet.out$core.out$core.archs)))			  
	sketch.graph = add.vertices(sketch.graph, nv = length(ACTIONet.out$core.out$core.archs), attr = attr)

	HSV = rgb2hsv(col2rgb(V(sketch.graph)$color))
	HSV[3, ] = HSV[3, ]*0.8
	V(sketch.graph)$frame.color = apply(HSV, 2, function(v) do.call(hsv, as.list(v)))

	coor = cbind(V(sketch.graph)$x, V(sketch.graph)$y)

	
	plot(sketch.graph, vertex.label=NA, layout=coor)
	
	
	## Add cell state graph on top
# 	pcors = suppressWarnings(ppcor::pcor(t(ACTIONet.out$reconstruct.out$H_stacked[ACTIONet.out$core.out$core.archs, ])))	
# 	CC = pcors$estimate
# 	
# 	CC[CC < min.cor] = 0
# 	CC = CC - diag(diag(CC))
# 	D = CC
# 	D[CC > 0] = 1 - D[CC > 0]
# 	core.graph.dist = graph_from_adjacency_matrix(D, mode = "undirected", weighted = TRUE)
# 	
# 	core.graph = mst(core.graph.dist, algorithm = "prim", weights = E(core.graph.dist)$weight)	
# 	#core.graph = graph_from_adjacency_matrix(CC, mode = "undirected", weighted = TRUE)
# 
# 
# 	V(core.graph)$color = ACTIONet.out$arch.vis.out$colors[ACTIONet.out$core.out$core.archs]
# 
# 	HSV = rgb2hsv(col2rgb(V(core.graph)$color))
# 	HSV[3, ] = HSV[3, ]*0.7
# 	V(core.graph)$frame.color = apply(HSV, 2, function(v) do.call(hsv, as.list(v)))
# 	
# 	V(core.graph)$size = 1.5*cex
# 
# 	E(core.graph)$color = ggplot2::alpha('black', 0.5)
# 	E(core.graph)$width = 2
# 
#   C = ACTIONet.out$reconstruct.out$C_stacked
#   arch.coors = t(sapply(1:dim(C)[2], function(col) {
#     cm.coor = ACTIONet.out$arch.vis.out$coordinates[col, ]
#     inf.cells = which(C[, col] > 0)
#     if(length(inf.cells) == 0)
#       return(c(0, 0))
#     
#     inf.cells.coor = ACTIONet.out$vis.out$coordinates[inf.cells, ]
#     if(length(inf.cells) == 1)
#       return(inf.cells.coor)
#     
#     delta = apply(inf.cells.coor, 1, function(x) sum(abs(x - cm.coor)))
#     rep.cell.coor = inf.cells.coor[which.min(delta), ]
#     
#     return(rep.cell.coor)
#   }))
#   core.arch.coors = arch.coors[ACTIONet.out$core.out$core.archs, ]  
# 
# 	plot(core.graph, vertex.label=NA, layout=core.arch.coors, add=T)

}


visualize.markers <- function(ACTIONet.out, sce, marker.genes, alpha_val = 0.9, node.size = 2, CPal = "d3", export_path = NA, thread_no = 8) {
	require(igraph)
	require(ACTIONet)
	require(viridis)
	require(ggpubr)
	
	
	if(!sum(sapply(marker.genes, length) != 1) & is.null(names(marker.genes))) {
		names(marker.genes) = marker.genes
	}
	 
	gg = unique(unlist(marker.genes))
	all.marker.genes = sort(intersect(gg, rownames(sce)))

	imputed.marker.expression = impute.genes.using.ACTIONet(ACTIONet.out, sce, all.marker.genes, alpha_val, thread_no, prune = TRUE)
	
	sketch.graph = ACTIONet.out$ACTIONet
	sketch.graph = delete.edges(sketch.graph, E(sketch.graph))
	V(sketch.graph)$size = 2.5
	coor = cbind(V(sketch.graph)$x, V(sketch.graph)$y)


	if(!(CPal %in% c("inferno", "magma", "viridis", "BlGrRd"))) {
		lg = hsv(h = 0, s = 0, v = 0.95)
		Pal = ggpubr::get_palette(CPal, length(names(marker.genes)))
		names(Pal) = names(marker.genes)
	}
	
	V(sketch.graph)$size = node.size
	lapply(all.marker.genes, function(gene) {
		if(! (gene %in% colnames(imputed.marker.expression)) )
			return()
			
		idx = which(sapply(marker.genes, function(gs) gene %in% gs))[1]
		celltype.name = names(marker.genes)[idx]

		if(CPal %in% c("inferno", "magma", "viridis", "BlGrRd")) {
			Pal_grad = switch(CPal, 
			"inferno" = inferno(500, alpha = 0.5),
			"magma" = magma(500, alpha = 0.5),
			"viridis" = viridis(500, alpha = 0.5),
			"BlGrRd" = colorRampPalette(c('blue', 'grey', 'red'))(500))
		} else {
			Pal_grad = colorRampPalette(c(lg, Pal[celltype.name]), bias = 100)(500)
		}
		x = imputed.marker.expression[, gene]
		if(sum(x) == 0)
			return()
		
		vCol = rep(Pal_grad[[1]], length(x))
		#vCol[x > 0] = scales::col_bin(Pal_grad, domain = NULL, bins = 100)(x[x > 0])
		v = x[x > 0]
		v = 1 / (1 + exp(-scale(v)))
		vCol[x > 0] = scales::col_numeric(Pal_grad, domain = NULL)(v)
		V(sketch.graph)$color = vCol


		HSV = rgb2hsv(col2rgb(vCol))
		HSV[3, ] = HSV[3, ]*0.7
		vCol.border = apply(HSV, 2, function(v) do.call(hsv, as.list(v)))
		V(sketch.graph)$frame.color = vCol.border
		
		plot(sketch.graph, vertex.label=NA, layout=coor, main = ifelse(celltype.name == gene, gene, sprintf('%s (%s)', celltype.name, gene)))

		if(!is.na(export_path)) {
			fname = sprintf('%s/%s.pdf', export_path, ifelse(celltype.name == gene, gene, sprintf('%s_%s', celltype.name, gene)));
			pdf(fname)
			plot(sketch.graph, vertex.label=NA, layout=coor, main = gene)
			dev.off()
		}
	});
}


visualize.ADTs <- function(ACTIONet.out, sce, alpha_val = 0.95, node.size = 2, CPal = "inferno", export_path = NA) {
	require(igraph)
	require(ACTIONet)
	require(viridis)
	require(ggpubr)
	
	if(! ('ADT' %in% names(sce@reducedDims)) ) {
		print("ADT is not a stored in the input sce");
		return()
	}
	
	G = get.adjacency(ACTIONet.out$ACTIONet, attr = "weight")
	 
	U = log(1 + sce@reducedDims[["ADT"]])
	U = scale(U, center = FALSE, scale = Matrix::colSums(U))

	smoothed.scores = batchPR(G, U, alpha_val, 8)	
	colnames(smoothed.scores) = colnames(sce@reducedDims[["ADT"]])

	sketch.graph = ACTIONet.out$ACTIONet
	sketch.graph = delete.edges(sketch.graph, E(sketch.graph))
	V(sketch.graph)$size = 2.5
	coor = cbind(V(sketch.graph)$x, V(sketch.graph)$y)

	
	V(sketch.graph)$size = node.size
	lapply(colnames(smoothed.scores), function(Ab) {

		Pal_grad = switch(CPal, 
		"inferno" = inferno(500, alpha = 0.5),
		"magma" = magma(500, alpha = 0.5),
		"viridis" = viridis(500, alpha = 0.5),
		"BlGrRd" = colorRampPalette(c('blue', 'grey', 'red'))(500))
		
		x = smoothed.scores[, Ab]

		nnz = sum(x^2)^2/sum(x^4)
		perm = order(x, decreasing = TRUE)
		mask = (1:length(perm)) %in% perm[1:nnz]
		
		#mask = LFR[, gene] > LFR.threshold
		
		u = array(0, length(x))
		u[mask] = scale(x[mask]) #1 / (1 + exp(-0.5*scale(x[mask])))


		V(sketch.graph)$color[mask] = scales::col_bin(Pal_grad, domain = NULL, bins = 100)(u[mask])
		V(sketch.graph)$color[!mask] =  Pal_grad[[1]]


		V(sketch.graph)$frame.color = V(sketch.graph)$color#apply(HSV, 2, function(v) 
		plot(sketch.graph, vertex.label=NA, layout=coor, main = sprintf('%s', Ab))

		if(!is.na(export_path)) {
			fname = sprintf('%s/ADT_%s.pdf', export_path, Ab);
			pdf(fname)
			plot(sketch.graph, vertex.label=NA, layout=coor, main = gene)
			dev.off()
		}
	})
}


plot.ACTIONet.gene.view <- function(ACTIONet.out, top.gene.count = 10, blacklist.pattern = '\\.|^RPL|^RPS|^MRP|MT-') {
	signature.profile = ACTIONet.out$signature.profile[, ACTIONet.out$core.out$core.archs]

	filtered.row.mask = grepl(blacklist.pattern, toupper(rownames(sce)))
	signature.profile = signature.profile[!filtered.row.mask, ]

	sorted.top.genes = apply(signature.profile, 2, function(x) rownames(signature.profile)[order(x, decreasing = TRUE)[1:top.gene.count]])


	selected.genes = sort(unique(as.character(sorted.top.genes)))

	arch.RGB = col2rgb(ACTIONet.out$arch.vis.out$colors[ACTIONet.out$core.out$core.archs]) / 256;
	arch.Lab = grDevices::convertColor(color= t(arch.RGB), from = 'sRGB', to = 'Lab')

	gene.color.Lab = t(sapply(selected.genes, function(gene) {
	  v = as.numeric(apply(sorted.top.genes, 2, function(gg) gene %in% gg))
	  v = v / sum(v)
	  gene.Lab = t(v) %*% arch.Lab
	  return(gene.Lab)
	}))
	gene.colors = rgb(grDevices::convertColor(color=gene.color.Lab, from = 'Lab', to = 'sRGB'))


	arch.coor = ACTIONet.out$arch.vis.out$coordinates[ACTIONet.out$core.out$core.archs, ]
	gene.coors = t(sapply(selected.genes, function(gene) {
	  v = as.numeric(apply(sorted.top.genes, 2, function(gg) gene %in% gg))
	  v = v / sum(v)
	  gene.coor = t(v) %*% arch.coor
	  return(gene.coor)
	}))


	genes.df = data.frame(gene = selected.genes, x = gene.coors[, 1], y = gene.coors[, 2])


	require(ggrepel)
	require(ggplot2)
	p <- ggplot(genes.df, aes(x, y, label = gene, color=gene)) + scale_colour_manual(values=gene.colors) + geom_point(show.legend = FALSE) + geom_label_repel(show.legend = FALSE) + theme_void()

	plot(p)	
}

plot.ACTIONet.interactive <- function(ACTIONet.out, sce, labels = NULL, top.arch.genes = 10, blacklist.pattern = '\\.|^RPL|^RPS|^MRP|MT-', marker.per.cell = 5, node.size = 2, CPal = "d3") {
	
	require(plotly)
	require(ACTIONet)

	signature.profile = ACTIONet.out$signature.profile[, ACTIONet.out$core.out$core.archs]

	filtered.row.mask = grepl(blacklist.pattern, toupper(rownames(sce)))
	signature.profile = signature.profile[!filtered.row.mask, ]

	sorted.top.genes = apply(signature.profile, 2, function(x) rownames(signature.profile)[order(x, decreasing = TRUE)[1:top.arch.genes]])

	selected.genes = sort(unique(as.character(sorted.top.genes)))

	#imputed.markers = impute.genes.using.ACTIONet(ACTIONet.out, sce, selected.genes, prune = TRUE, rescale = TRUE)
	imputed.markers = t(ACTIONet.out$signature.profile[selected.genes, ACTIONet.out$core.out$core.archs] %*% ACTIONet.out$core.out$H)

	node.annotations = apply(imputed.markers, 1, function(x) {
		top.genes = colnames(imputed.markers)[order(x, decreasing = TRUE)[1:marker.per.cell]]
		label = paste(top.genes, collapse = ';')
		return(label)
	})
	
	if(is.null(labels)) {
		labels = as.factor(array(1, dim(sce)[2]))
	} else if(is.character(labels) | is.numeric(labels)) {
		labels = factor(labels, levels = sort(unique(labels)))
	}
	
	# Setup visualization parameters
	sketch.graph = ACTIONet.out$ACTIONet
	sketch.graph = delete.edges(sketch.graph, E(sketch.graph))

	node.data <- get.data.frame(sketch.graph, what="vertices") 
	edge.data <- get.data.frame(sketch.graph, what="edges") 

	# Adjust parameters
	node.data$size = node.size
	Pal = ggpubr::get_palette(CPal, length(levels(labels)))
	names(Pal) = levels(labels)
	node.data$color = Pal[labels]

	Nv <- dim(node.data)[1]
	Ne <- dim(edge.data)[1]

	edge_shapes <- list()

	network <- plot_ly(node.data, x = ~x, y = ~y, marker = list(size = ~size, color = ~color, line = list(width = 0)), text = node.annotations, mode = "markers", type = 'scatter', hoverinfo = "text")

	axis <- list(title = "", showgrid = FALSE, showticklabels = FALSE, zeroline = FALSE)

	p <- plotly::layout(
	  network,
	  title = 'ACTIONet',
	  shapes = edge_shapes,
	  xaxis = axis,
	  yaxis = axis,
	  showlegend=FALSE
	)
}

## Function to draw Matrix
plotHeatMat <- function(x, ...){
  min <- min(x)
  max <- max(x)
  col.labels <- rownames(x)
  row.labels <- colnames(x)
  title <-c()
  
  # check for additional function arguments
  if( length(list(...)) ){
    Lst <- list(...)
    if( !is.null(Lst$zlim) ){
      min <- Lst$zlim[1]
      max <- Lst$zlim[2]
    }
    if( !is.null(Lst$col.labels) ){
      col.labels <- c(Lst$col.labels)
    }
    if( !is.null(Lst$row.labels) ){
      row.labels <- c(Lst$row.labels)
    }
    if( !is.null(Lst$title) ){
      title <- Lst$title
    }
  }
  # # check for null values
  # if( is.null(row.labels) ){
  #    row.labels <- c(1:ncol(x))
  # }
  # if( is.null(col.labels) ){
  #    col.labels <- c(1:nrow(x))
  # }
  
  layout(matrix(data=c(1,2), nrow=1, ncol=2), widths=c(4,1), heights=c(1,1))
  
  # Red and green range from 0 to 1 while Blue ranges from 1 to 0
  ColorRamp <- rgb( seq(0,1,length=256),  # Red
                    seq(0,1,length=256),  # Green
                    seq(1,0,length=256))  # Blue
  ColorLevels <- seq(min, max, length=length(ColorRamp))
  
  # Reverse Y axis
  reverse <- nrow(x) : 1
  col.labels <- col.labels[reverse]
  x <- x[reverse,]
  
  # Data Map
  par(mar = c(3,5,2.5,2))
  image(1:dim(x)[1], 1:dim(x)[2], x, col=ColorRamp, xlab="", ylab="", axes=FALSE, zlim=c(min,max))
  if( !is.null(title) ){
    title(main=title)
  }
  if(!is.na(row.labels)) {
    axis(BELOW<-1, at=1:length(row.labels), labels=row.labels, cex.axis=0.7)
  }
  if(!is.na(row.labels)) {
    axis(LEFT <-2, at=1:length(col.labels), labels=col.labels, las= HORIZONTAL<-1, cex.axis=0.7)
  }
  
  
  # Color Scale
  par(mar = c(3,2.5,2.5,2))
  image(1, ColorLevels,
        matrix(data=ColorLevels, ncol=length(ColorLevels),nrow=1),
        col=ColorRamp,
        xlab="",ylab="",
        xaxt="n")
  
  layout(1)
}
