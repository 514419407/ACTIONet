mycircle <- function(coords, v=NULL, params) {
  library(igraph)
  
  vertex.color <- params("vertex", "color")
  if (length(vertex.color) != 1 && !is.null(v)) {
    vertex.color <- vertex.color[v]
  }
  vertex.size  <- 1/200 * params("vertex", "size")
  if (length(vertex.size) != 1 && !is.null(v)) {
    vertex.size <- vertex.size[v]
  }
  vertex.frame.color <- params("vertex", "frame.color")
  if (length(vertex.frame.color) != 1 && !is.null(v)) {
    vertex.frame.color <- vertex.frame.color[v]
  }
  vertex.frame.width <- params("vertex", "frame.width")
  if (length(vertex.frame.width) != 1 && !is.null(v)) {
    vertex.frame.width <- vertex.frame.width[v]
  }
  
  mapply(coords[,1], coords[,2], vertex.color, vertex.frame.color,
         vertex.size, vertex.frame.width,
         FUN=function(x, y, bg, fg, size, lwd) {
           symbols(x=x, y=y, bg=bg, fg=fg, lwd=lwd,
                   circles=size, add=TRUE, inches=FALSE)
         })
}
add.vertex.shape("fcircle", clip=igraph.shape.noclip,
		 plot=mycircle, parameters=list(vertex.frame.color=1,
                                  vertex.frame.width=1))

plot.ACTIONet <- function(ACTIONet.out, color.attr = NA, transparency.attr = NA, size.attr = NA, cex = 1, CPal = "Dark2", legend.pos = "bottomright") {
	require(ggplot2)
	require(ggpubr)
	
	if(is.igraph(ACTIONet.out))
		ACTIONet = ACTIONet.out
	else
		ACTIONet = ACTIONet.out$ACTIONet
		
	nV = length(V(ACTIONet))
	coor = cbind(V(ACTIONet)$x, V(ACTIONet)$y)


	Annot = NA
	if(is.numeric(color.attr)) {
		color.attr = as.character(color.attr)
		color.attr = factor(color.attr, levels = sort(unique(color.attr)))
		Annot = levels(color.attr)
		
		if(is.list(CPal)) {
			Pal = CPal[1:length(Annot)]
		} else {
			Pal = ggpubr::get_palette(CPal, length(Annot))
		}
		names(Pal) = Annot
	}
	else if(is.factor(color.attr)) {
		Annot = levels(color.attr)
		if(is.list(CPal)) {
			Pal = CPal[1:length(Annot)]
		} else {
			Pal = ggpubr::get_palette(CPal, length(Annot))
		}
		names(Pal) = Annot
	} 
	else if(is.character(color.attr)) {
		color.attr = factor(color.attr, levels = sort(unique(color.attr)))
		Annot = levels(color.attr)
		if(is.list(CPal)) {
			Pal = CPal[1:length(Annot)]
		} else {
			Pal = ggpubr::get_palette(CPal, length(Annot))
		}
		names(Pal) = Annot
	}		
	if(length(Annot) > 1) {
		vCol = Pal[color.attr]
	} else {
		vCol = V(ACTIONet)$color
	}		


	if(is.numeric(color.attr)) {
		v = sort(unique(color.attr))
		if(length(v) < 30) {
			Annot = as.character(v)
			if(is.list(CPal)) {
				Pal = CPal[1:length(Annot)]
			} else {
				Pal = ggpubr::get_palette(CPal, length(Annot))
			}
			names(Pal) = Annot
			vCol = Pal[color.attr]
		}
	}
	else if(is.factor(color.attr)) {
		Annot = levels(color.attr)
		if(is.list(CPal)) {
			Pal = CPal[1:length(Annot)]
		} else {
			Pal = ggpubr::get_palette(CPal, length(Annot))
		}
		names(Pal) = Annot
		vCol = Pal[color.attr]
	} 
	else if(is.character(color.attr)) {
		color.attr = factor(color.attr, levels = sort(unique(color.attr)))
		Annot = levels(color.attr)
		if(is.list(CPal)) {
			Pal = CPal[1:length(Annot)]
		} else {
			Pal = ggpubr::get_palette(CPal, length(Annot))
		}
		names(Pal) = Annot
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


plot.ACTIONet.igraph <- function(ACTIONet.out, color.attr = NA, transparency.attr = NA, size.attr = NA, cex = 1.5, trans.fact = 1.0, size.fact = 1.0, CPal = "d3", legend.pos = "bottomright") {
	border.width = 0.5*cex
	
	require(ggplot2)
	require(ggpubr)
	require(igraph)
	
	if(is.igraph(ACTIONet.out))
		ACTIONet = ACTIONet.out
	else
		ACTIONet = ACTIONet.out$ACTIONet
		
	nV = length(V(ACTIONet))
	#coor = cbind(V(ACTIONet)$x, V(ACTIONet)$y)


	sketch.graph = ACTIONet
	sketch.graph = delete.edges(sketch.graph, E(sketch.graph))

	V(sketch.graph)$name = ''
	V(sketch.graph)$shape = 'fcircle'
	
	

	Annot = NA
	if(is.numeric(color.attr)) {
		color.attr = as.character(color.attr)
		color.attr = factor(color.attr, levels = sort(unique(color.attr)))
		Annot = levels(color.attr)
		
		if(is.list(CPal)) {
			Pal = CPal[1:length(Annot)]
		} else {
			Pal = ggpubr::get_palette(CPal, length(Annot))
		}
		names(Pal) = Annot
	}
	else if(is.factor(color.attr)) {
		Annot = levels(color.attr)
		if(is.list(CPal)) {
			Pal = CPal[1:length(Annot)]
		} else {
			Pal = ggpubr::get_palette(CPal, length(Annot))
		}
		names(Pal) = Annot
	} 
	else if(is.character(color.attr)) {
		color.attr = factor(color.attr, levels = sort(unique(color.attr)))
		Annot = levels(color.attr)
		if(is.list(CPal)) {
			Pal = CPal[1:length(Annot)]
		} else {
			Pal = ggpubr::get_palette(CPal, length(Annot))
		}
		names(Pal) = Annot
	}		
	if(length(Annot) > 0) {
		vCol = Pal[color.attr]
	} else {
		vCol = V(ACTIONet)$color
	}		


	#HSV = rgb2hsv(col2rgb(vCol))
	#HSV[3, ] = HSV[3, ]*0.85
	#vCol.border = apply(HSV, 2, function(v) do.call(hsv, as.list(v)))
	
	if(is.numeric(transparency.attr)) {
		#beta = 1 / (1 + exp(-trans.fact*(scale(transparency.attr)+1)))
		z = (transparency.attr - median(transparency.attr)) / mad(transparency.attr)		
		beta = 1 / ( 1+ exp(-trans.fact*(z)) )
		#beta = (transparency.attr - min(transparency.attr)) / (max(transparency.attr) - min(transparency.attr))
		beta = beta ^ 2
		
		#plot(density(1-beta))
		#vCol = colorspace::lighten(vCol, 1-beta, method = "relative", space = "HCL")
		#vCol.border = colorspace::darken(vCol, 0.2*beta)#colorspace::lighten(vCol.border, 1-beta, method = "relative", space = "HLS")

		vCol = scales::alpha(vCol, beta)		
		#vCol.border = scales::alpha(vCol.border, beta)
		
		vCol.border = colorspace::darken(vCol, 0.5*beta)#colorspace::lighten(vCol.border, 1-beta, method = "relative", space = "HLS")
	
	} else {
		vCol.border = colorspace::darken(vCol, 0.3)
	}		



	if(is.numeric(size.attr)) {
		z = (size.attr - median(size.attr)) / mad(size.attr)		
		beta = 1 / ( 1+ exp(-size.fact*(z)) )
		beta = beta ^ 2
				
		#z = (size.attr - median(size.attr)) / mad(size.attr)		
		#size.factor = 1 / ( 1+ exp(-size.fact*z) )
		#size.factor = size.factor ^ 2

		cex = cex * beta
		 
	}


	
	V(sketch.graph)$size = cex
	V(sketch.graph)$color = vCol
	V(sketch.graph)$frame.color = vCol.border
	V(sketch.graph)$frame.width = border.width
	


	plot(sketch.graph)
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

	Annot = NA
	if(is.numeric(color.attr)) {
		color.attr = as.character(color.attr)
		color.attr = factor(color.attr, levels = sort(unique(color.attr)))
		Annot = levels(color.attr)
		
		if(is.list(CPal)) {
			Pal = CPal[1:length(Annot)]
		} else {
			Pal = ggpubr::get_palette(CPal, length(Annot))
		}
		names(Pal) = Annot
	}
	else if(is.factor(color.attr)) {
		Annot = levels(color.attr)
		if(is.list(CPal)) {
			Pal = CPal[1:length(Annot)]
		} else {
			Pal = ggpubr::get_palette(CPal, length(Annot))
		}
		names(Pal) = Annot		
	} 
	else if(is.character(color.attr)) {
		color.attr = factor(color.attr, levels = sort(unique(color.attr)))
		Annot = levels(color.attr)
		if(is.list(CPal)) {
			Pal = CPal[1:length(Annot)]
		} else {
			Pal = ggpubr::get_palette(CPal, length(Annot))
		}
		names(Pal) = Annot
	}		
	if(length(Annot) > 1) {
		vCol = as.character(Pal[color.attr])
	} else {
		vCol = V(ACTIONet)$color
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



plot.ACTIONet.cell.state.map <- function(ACTIONet.out, sce, archetype.labels = NA, transparency.attr = NA, trans.fact = 1, stretch.factor = 10, CPal = "d3", cex = 2, node.scale.factor = 3, reduction.slot = "S_r") {
	# Plot main ACTIONet first
	sketch.graph = ACTIONet.out$ACTIONet
	sketch.graph = delete.edges(sketch.graph, E(sketch.graph))
	nV = length(V(sketch.graph))


	vCol = colorspace::lighten('black', 0.97)
	vCol.border = colorspace::lighten('black', 0.9)
	
	if(is.numeric(transparency.attr)) {
		z = (transparency.attr - median(transparency.attr)) / mad(transparency.attr)		
		beta = 1 / ( 1+ exp(-trans.fact*(z)) )
		beta = beta ^ 2

		vCol = scales::alpha(vCol, beta)				
		vCol.border = scales::alpha(vCol.border, beta)				
	
	}	
	


	V(sketch.graph)$color = vCol
	V(sketch.graph)$frame.color = vCol.border
	
	
	
	
	
	V(sketch.graph)$size = cex

	coor = cbind(V(sketch.graph)$x, V(sketch.graph)$y)

	plot(sketch.graph, vertex.label=NA, layout=coor)
	
	
	# Now overlay the core backbone connectome on top
	core.index = ACTIONet.out$core.out$core.archs 
	core.coor = ACTIONet.out$arch.vis.out$coordinates[core.index, ]	
		
	Annot = NA
	if(is.numeric(archetype.labels)) {
		archetype.labels = as.character(archetype.labels)
		archetype.labels = factor(archetype.labels, levels = sort(unique(archetype.labels)))
		Annot = levels(archetype.labels)
		
		if(is.list(CPal)) {
			Pal = CPal[1:length(Annot)]
		} else {
			Pal = ggpubr::get_palette(CPal, length(Annot))
		}
		names(Pal) = Annot
	}
	else if(is.factor(archetype.labels)) {
		Annot = levels(archetype.labels)
		if(is.list(CPal)) {
			Pal = CPal[1:length(Annot)]
		} else {
			Pal = ggpubr::get_palette(CPal, length(Annot))
		}
		names(Pal) = Annot
	} 
	else if(is.character(archetype.labels)) {
		archetype.labels = factor(archetype.labels, levels = sort(unique(archetype.labels)))
		Annot = levels(archetype.labels)
		if(is.list(CPal)) {
			Pal = CPal[1:length(Annot)]
		} else {
			Pal = ggpubr::get_palette(CPal, length(Annot))
		}
		names(Pal) = Annot
	}		
	if(length(Annot) > 1) {
		vCol = Pal[archetype.labels]
	} else {
		vCol = ACTIONet.out$arch.vis.out$colors
	}		
		
	core.col = vCol[core.index]

	cs = as.numeric(table(ACTIONet.out$core.out$arch.membership))
	core.scale.factor = cs / max(cs)
	core.size = cex*(1 + node.scale.factor*core.scale.factor)


	backbone.graph = construct.sparse.backbone(ACTIONet.out, reduction.slot, stretch.factor)


	w.scaled = E(backbone.graph)$weight; w.scaled = (w.scaled / max(w.scaled))

	E(backbone.graph)$width = 0.5 + 2.5*w.scaled
	E(backbone.graph)$color = ggplot2::alpha('black', 0.25 + 0.35*w.scaled)
	E(backbone.graph)$arrow.size = 0.1

	V(backbone.graph)$color = core.col
	V(backbone.graph)$size = core.size

	HSV = rgb2hsv(col2rgb(V(backbone.graph)$color))
	HSV[3, ] = HSV[3, ]*0.8
	V(backbone.graph)$frame.color = apply(HSV, 2, function(v) do.call(hsv, as.list(v)))


	plot(backbone.graph, vertex.label=NA, layout=core.coor, add = T)
	if(!is.na(archetype.labels)) {
		legend('bottomright', legend = Annot, fill=Pal, cex = 0.5)
	}	
}



plot.ACTIONet.cell.state.view <- function(ACTIONet.out, archetype.labels = NA, transparency.attr = NA, trans.fact = 1, CPal = "d3", cex = 2, node.scale.factor = 3) {
	# Plot main ACTIONet first
	sketch.graph = ACTIONet.out$ACTIONet
	sketch.graph = delete.edges(sketch.graph, E(sketch.graph))
	nV = length(V(sketch.graph))


	vCol = colorspace::lighten('black', 0.97)
	vCol.border = colorspace::lighten('black', 0.9)
	
	if(is.numeric(transparency.attr)) {
		z = (transparency.attr - median(transparency.attr)) / mad(transparency.attr)		
		beta = 1 / ( 1+ exp(-trans.fact*(z)) )
		beta = beta ^ 2

		vCol = scales::alpha(vCol, beta)				
		vCol.border = scales::alpha(vCol.border, beta)				
	
	}	
	


	V(sketch.graph)$color = vCol
	V(sketch.graph)$frame.color = vCol.border
	
	
	
	
	
	V(sketch.graph)$size = cex

	coor = cbind(V(sketch.graph)$x, V(sketch.graph)$y)

	plot(sketch.graph, vertex.label=NA, layout=coor)
	
	
	# Now overlay the core backbone connectome on top
	core.index = ACTIONet.out$core.out$core.archs 
	core.coor = ACTIONet.out$arch.vis.out$coordinates[core.index, ]	
		
	Annot = NA
	if(is.numeric(archetype.labels)) {
		archetype.labels = as.character(archetype.labels)
		archetype.labels = factor(archetype.labels, levels = sort(unique(archetype.labels)))
		Annot = levels(archetype.labels)
		
		if(is.list(CPal)) {
			Pal = CPal[1:length(Annot)]
		} else {
			Pal = ggpubr::get_palette(CPal, length(Annot))
		}
		names(Pal) = Annot
	}
	else if(is.factor(archetype.labels)) {
		Annot = levels(archetype.labels)
		if(is.list(CPal)) {
			Pal = CPal[1:length(Annot)]
		} else {
			Pal = ggpubr::get_palette(CPal, length(Annot))
		}
		names(Pal) = Annot
	} 
	else if(is.character(archetype.labels)) {
		archetype.labels = factor(archetype.labels, levels = sort(unique(archetype.labels)))
		Annot = levels(archetype.labels)
		if(is.list(CPal)) {
			Pal = CPal[1:length(Annot)]
		} else {
			Pal = ggpubr::get_palette(CPal, length(Annot))
		}
		names(Pal) = Annot
	}		
	if(length(Annot) > 1) {
		vCol = Pal[archetype.labels]
	} else {
		vCol = ACTIONet.out$arch.vis.out$colors
	}		
		
	core.col = vCol[core.index]

	cs = as.numeric(table(ACTIONet.out$core.out$arch.membership))
	core.scale.factor = cs / max(cs)
	core.size = cex*(1 + node.scale.factor*core.scale.factor)


	backbone.graph = graph.empty(n = length(core.size), directed = FALSE)	
	V(backbone.graph)$color = core.col
	V(backbone.graph)$size = core.size

	HSV = rgb2hsv(col2rgb(V(backbone.graph)$color))
	HSV[3, ] = HSV[3, ]*0.8
	V(backbone.graph)$frame.color = apply(HSV, 2, function(v) do.call(hsv, as.list(v)))


	plot(backbone.graph, vertex.label=NA, layout=core.coor, add = T)
	if(!is.na(archetype.labels)) {
		legend('bottomright', legend = Annot, fill=Pal, cex = 0.5)
	}
}

plot.ACTIONet.gene.view <- function(ACTIONet.out, sce, top.gene.count = 10, blacklist.pattern = '\\.|^RPL|^RPS|^MRP|^MT-|^MT|^RP|MALAT1|B2M') {
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
	p <- ggplot(genes.df, aes(x, y, label = gene, color=gene)) + scale_colour_manual(values=gene.colors) + geom_point(show.legend = FALSE) + geom_label_repel(show.legend = FALSE, force = 5) + theme_void()

	plot(p)	
}


plot.ACTIONet.phenotype.view <- function(ACTIONet.out, phenotype.scores) {
	phenptype.coors = W %*% ACTIONet.out$arch.vis.out$coordinates[ACTIONet.out$core.out$core.archs, ]
	
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
	p <- ggplot(genes.df, aes(x, y, label = gene, color=gene)) + scale_colour_manual(values=gene.colors) + geom_point(show.legend = FALSE) + geom_label_repel(show.legend = FALSE, force = 5) + theme_void()

	plot(p)	
}



# define utility function to adjust fill-opacity using css
fillOpacity <- function(., alpha = 0.5) {
  css <- sprintf("<style> .js-fill { fill-opacity: %s !important; } </style>", alpha)
  prependContent(., HTML(css))
}

plot.ACTIONet.interactive <- function(ACTIONet.out, sce, labels = NULL, top.arch.genes = 10, blacklist.pattern = '\\.|^RPL|^RPS|^MRP|^MT-|^MT|^RP|MALAT1|B2M', marker.per.cell = 5, node.size = 2, CPal = "Spectral", show.legend = TRUE, annotate.cells = TRUE, opacity = 1.0, title = 'ACTIONet') {
	
	require(plotly)
	require(ACTIONet)

	signature.profile = ACTIONet.out$signature.profile[, ACTIONet.out$core.out$core.archs]

	filtered.row.mask = grepl(blacklist.pattern, toupper(rownames(sce)))
	signature.profile = signature.profile[!filtered.row.mask, ]

	sorted.top.genes = apply(signature.profile, 2, function(x) rownames(signature.profile)[order(x, decreasing = TRUE)[1:top.arch.genes]])

	selected.genes = sort(unique(as.character(sorted.top.genes)))

	if(annotate.cells == TRUE) {
		#imputed.markers = impute.genes.using.ACTIONet(ACTIONet.out, sce, selected.genes, prune = TRUE, rescale = TRUE)
		imputed.markers = t(ACTIONet.out$signature.profile[selected.genes, ACTIONet.out$core.out$core.archs] %*% ACTIONet.out$core.out$H)
		node.annotations = apply(imputed.markers, 1, function(x) {
			top.genes = colnames(imputed.markers)[order(x, decreasing = TRUE)[1:marker.per.cell]]
			label = paste(top.genes, collapse = ';')
			return(label)
		})
	} else {
		node.annotations = ''
	}
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

	Nv <- dim(node.data)[1]
	Ne <- dim(edge.data)[1]

	edge_shapes <- list()

	# Adjust parameters
	node.data$size = node.size

	axis <- list(title = "", showgrid = FALSE, showticklabels = FALSE, zeroline = FALSE)


	if(is.numeric(labels)) {
		labels = as.character(labels)
		labels = factor(labels, levels = sort(unique(labels)))
		Annot = levels(labels)
		
		if(is.list(CPal)) {
			Pal = CPal[1:length(Annot)]
		} else {
			Pal = ggpubr::get_palette(CPal, length(Annot))
		}
		names(Pal) = Annot
	}
	else if(is.factor(labels)) {
		Annot = levels(labels)
		if(is.list(CPal)) {
			Pal = CPal[1:length(Annot)]
		} else {
			Pal = ggpubr::get_palette(CPal, length(Annot))
		}
		names(Pal) = Annot
	} 
	else {
		labels = factor(labels, levels = sort(unique(labels)))
		Annot = levels(labels)
		if(is.list(CPal)) {
			Pal = CPal[1:length(Annot)]
		} else {
			Pal = ggpubr::get_palette(CPal, length(Annot))
		}
		names(Pal) = Annot
	}		
	#vCol = Pal[labels]
	
	node.data$type = labels
	network <- plot_ly(node.data, x = ~x, y = ~y, opacity = opacity, color = ~type, colors = Pal, marker = list(size = ~size, opacity = 1.0, alpha = 1, line = list(width = 0.1*node.size, alpha = 0.5, color = 'rgb(0, 0, 0)')), text = node.annotations, mode = "markers", type = 'scatter', hoverinfo = "text")
	
	p <- plotly::layout(
	  network,
	  title = title,
	  shapes = edge_shapes,
	  xaxis = axis,
	  yaxis = axis,
	  showlegend=show.legend
	) 
}

plot.ACTIONet.interactive.3D <- function(ACTIONet.out, sce, labels = NULL, top.arch.genes = 10, blacklist.pattern = '\\.|^RPL|^RPS|^MRP|^MT-|^MT|^RP|MALAT1|B2M', marker.per.cell = 5, node.size = 2, CPal = "Spectral", show.legend = TRUE, annotate.cells = TRUE, opacity = 1.0, title = 'ACTIONet') {
	
	require(plotly)
	require(ACTIONet)

	signature.profile = ACTIONet.out$signature.profile[, ACTIONet.out$core.out$core.archs]

	filtered.row.mask = grepl(blacklist.pattern, toupper(rownames(sce)))
	signature.profile = signature.profile[!filtered.row.mask, ]

	sorted.top.genes = apply(signature.profile, 2, function(x) rownames(signature.profile)[order(x, decreasing = TRUE)[1:top.arch.genes]])

	selected.genes = sort(unique(as.character(sorted.top.genes)))

	if(annotate.cells == TRUE) {
		#imputed.markers = impute.genes.using.ACTIONet(ACTIONet.out, sce, selected.genes, prune = TRUE, rescale = TRUE)
		imputed.markers = t(ACTIONet.out$signature.profile[selected.genes, ACTIONet.out$core.out$core.archs] %*% ACTIONet.out$core.out$H)
		node.annotations = apply(imputed.markers, 1, function(x) {
			top.genes = colnames(imputed.markers)[order(x, decreasing = TRUE)[1:marker.per.cell]]
			label = paste(top.genes, collapse = ';')
			return(label)
		})
	} else {
		node.annotations = ''
	}
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
	
	Nv <- dim(node.data)[1]
	Ne <- dim(edge.data)[1]

	edge_shapes <- list()

	# Adjust parameters
	node.data$size = node.size


	if(is.numeric(labels)) {
		labels = as.character(labels)
		labels = factor(labels, levels = sort(unique(labels)))
		Annot = levels(labels)
		
		if(is.list(CPal)) {
			Pal = CPal[1:length(Annot)]
		} else {
			Pal = ggpubr::get_palette(CPal, length(Annot))
		}
		names(Pal) = Annot
	}
	else if(is.factor(labels)) {
		Annot = levels(labels)
		if(is.list(CPal)) {
			Pal = CPal[1:length(Annot)]
		} else {
			Pal = ggpubr::get_palette(CPal, length(Annot))
		}
		names(Pal) = Annot
	} 
	else {
		labels = factor(labels, levels = sort(unique(labels)))
		Annot = levels(labels)
		if(is.list(CPal)) {
			Pal = CPal[1:length(Annot)]
		} else {
			Pal = ggpubr::get_palette(CPal, length(Annot))
		}
		names(Pal) = Annot
	}		
	#vCol = Pal[labels]

	
	node.data$type = labels
	
	
	network <- plot_ly(node.data, x = ~x3D, y = ~y3D, z = ~z3D, opacity = opacity, color = ~type, colors = Pal, marker = list(size = ~size, opacity = 1.0, alpha = 1, line = list(width = 0.1*node.size, alpha = 0.5, color = 'rgb(0, 0, 0)')), text = node.annotations, mode = "markers", hoverinfo = "text", type = "scatter3d")
	
	axis <- list(title = "", showgrid = FALSE, showticklabels = FALSE, zeroline = FALSE)
	p <- plotly::layout(
	  network,
	  title = title,
	  shapes = edge_shapes,
	  scene = list(
	  xaxis = axis,
	  yaxis = axis,
	  zaxis = axis),
	  showlegend=show.legend
	) 
}


plot.marker.boxplot <- function(ACTIONet.out, sce, marker.genes, Labels, node.size = 3, CPal = "d3", export_path = NA, thread_no = 8, prune = FALSE, scale.factor = 2) {
	require(igraph)
	require(ACTIONet)
	require(ggpubr)
	
	if(!is.factor(Labels)) {
		Labels = factor(Labels, levels = sort(unique(Labels)))
	}
		
	if(!sum(sapply(marker.genes, length) != 1) & is.null(names(marker.genes))) {
		names(marker.genes) = marker.genes
	}
	 
	gg = unique(unlist(marker.genes))
	all.marker.genes = sort(intersect(gg, rownames(sce)))

	imputed.markers = impute.genes.using.ACTIONet(ACTIONet.out, sce, all.marker.genes, prune = prune)

	imputed.markers.df = as.data.frame(imputed.markers * nrow(imputed.markers))
	imputed.markers.df$Celltype = Labels
	
	
	if(is.list(CPal)) {
		Pal = CPal[1:length(levels(Labels))]
	} else {
		Pal = ggpubr::get_palette(CPal, length(levels(Labels)))
	}
	names(Pal) = levels(Labels)

	sapply(colnames(imputed.markers), function(gene.name) {
		gp = ggboxplot(imputed.markers.df, x = "Celltype", y = gene.name, fill = "Celltype", palette = Pal) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
		print(gp)
				
		if(!is.na(export_path)) {
			fname = sprintf('%s/%s.pdf', export_path, gene.name);
			pdf(fname)
			print(gp)			
			dev.off()
		}		
		
	})

}


plot.marker.dotplot <- function(ACTIONet.out, sce, marker.genes, Labels, CPal = "YlOrRd", thread_no = 8, prune = FALSE) {
	require(ACTIONet)
	library(corrplot)
	library(seriation)
	
	if(!is.factor(Labels)) {
		Labels = factor(Labels, levels = sort(unique(Labels)))
	}
		
	if(!sum(sapply(marker.genes, length) != 1) & is.null(names(marker.genes))) {
		names(marker.genes) = marker.genes
	}
	 
	gg = unique(unlist(marker.genes))
	all.marker.genes = sort(intersect(gg, rownames(sce)))

	imputed.markers = impute.genes.using.ACTIONet(ACTIONet.out, sce, all.marker.genes, prune = prune)


	X = apply(imputed.markers, 2, function(x) return( (x - min(x)) / (max(x) - min(x))))

	IDX = split(1:nrow(imputed.markers), Labels)
	mean.expr = sapply(IDX, function(idx) as.numeric(Matrix::colMeans(X[idx, ])))
	rownames(mean.expr) = colnames(imputed.markers)

	set.seed(0)
	perm = seriation::get_order(seriate(mean.expr, "BEA_TSP"))


	Pal = ggpubr::get_palette(CPal, 11)

	corrplot(mean.expr[perm,], is.corr = FALSE, method = "circle", tl.col="black", cl.pos="n", col = Pal)
}

plot.ACTIONet.gradient <- function(ACTIONet.out, x, max.update.iter = 3, CPal = "tomato", node.size = 3) {

	ACTIONet = ACTIONet.out$ACTIONet


	sketch.graph = ACTIONet
	sketch.graph = delete.edges(sketch.graph, E(sketch.graph))

	V(sketch.graph)$name = ''
	V(sketch.graph)$shape = 'fcircle'


	eps = 1e-16
	A = as(get.adjacency(ACTIONet, attr = "weight"), 'dgTMatrix')
	rs = Matrix::rowSums(A)
	P = sparseMatrix(i = A@i+1, j = A@j+1, x = A@x/rs[A@i+1], dims = dim(A))  


	# Normalize between [0, 1]
	x = (x - min(x)) / (max(x) - min(x))

	# Find effective # of nonzeros using participation ratio
	nnz = round(sum(abs(x))^2 / sum(x^2))
	l = array(0, length(x))
	l[order(x, decreasing = TRUE)[1:nnz]] = 1

	# Use label propagation to smooth nonzero estimates
	for(it in 1:max.update.iter) {
		p = mean(l)

		X = sapply(unique(l), function(i) {
			x = as.numeric(Matrix::sparseVector(x = 1, i = which(l == i), length = length(l)))
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
		l = apply(logPval, 1, which.max)-1
	}


	V(sketch.graph)$color = rgb(0.95, 0.95, 0.95)
	V(sketch.graph)$frame.color = rgb(0.9, 0.9, 0.9)


	if(CPal %in% c("inferno", "magma", "viridis", "BlGrRd")) {
		Pal_grad = switch(CPal, 
		"inferno" = inferno(500, alpha = 0.8),
		"magma" = magma(500, alpha = 0.8),
		"viridis" = viridis(500, alpha = 0.8),
		"BlGrRd" = colorRampPalette(c('blue', 'grey', 'red'))(500))
		
	} else {
		lg = rgb(0.95, 0.95, 0.95)
		Pal_grad = colorRampPalette(c(lg, CPal))(500)
	}

	V(sketch.graph)$color[l == 1] = darken(scales::col_numeric(Pal_grad, domain = NULL)(scale(x[l == 1])), 0.25)
	V(sketch.graph)$frame.color[l == 1] = darken(V(sketch.graph)$color[l == 1], 0.2)


	V(sketch.graph)$size = node.size
	V(sketch.graph)$frame.width = node.size*0.1

	plot(sketch.graph)
}

visualize.markers <- function(ACTIONet.out, sce, marker.genes, max.update.iter = 3, CPal = "d3", node.size = 3, adjust.node.size = TRUE, alpha_val = 0.85, export_path = NA) {
	if(!sum(sapply(marker.genes, length) != 1) & is.null(names(marker.genes))) {
		names(marker.genes) = marker.genes
	}
	Pal = ggpubr::get_palette(CPal, length(names(marker.genes)))
	names(Pal) = names(marker.genes)
	
	gg = unique(unlist(marker.genes))
	all.marker.genes = sort(intersect(gg, rownames(sce)))
	
	imputed.marker.expression = impute.genes.using.ACTIONet(ACTIONet.out, sce, all.marker.genes, prune = prune, alpha_val = alpha_val)


	ACTIONet = ACTIONet.out$ACTIONet
	sketch.graph = ACTIONet
	sketch.graph = delete.edges(sketch.graph, E(sketch.graph))

	V(sketch.graph)$name = ''
	V(sketch.graph)$shape = 'fcircle'


	eps = 1e-16
	A = as(get.adjacency(ACTIONet, attr = "weight"), 'dgTMatrix')
	rs = Matrix::rowSums(A)
	P = sparseMatrix(i = A@i+1, j = A@j+1, x = A@x/rs[A@i+1], dims = dim(A))  

	gradPal = grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu")))(100)
	
	lapply(all.marker.genes, function(gene) {
		print(gene)
		if(! (gene %in% colnames(imputed.marker.expression)) )
			return()
			
		idx = which(sapply(marker.genes, function(gs) gene %in% gs))[1]
		celltype.name = names(marker.genes)[idx]
		
		x = imputed.marker.expression[, gene]


		# Normalize between [0, 1]
		x = (x - min(x)) / (max(x) - min(x))

		# Find effective # of nonzeros using participation ratio
		nnz = round(sum(abs(x^2))^2 / sum(x^4))
		l = array(0, length(x))
		l[order(x, decreasing = TRUE)[1:nnz]] = 1

		# Use label propagation to smooth nonzero estimates
		if(max.update.iter > 0) {
			for(it in 1:max.update.iter) {			
				p = mean(l)

				X = sapply(0:1, function(i) {
					x = as.numeric(Matrix::sparseVector(x = 1, i = which(l == i), length = length(l)))
				})


				Exp = array(1, nrow(A)) %*% t(p)
				Obs = P %*% X 

				Lambda = Obs - Exp


				w2 = Matrix::rowSums(P^2)
				Nu = w2 %*% t(p)

				a = as.numeric(qlcMatrix::rowMax(P)) %*% t(array(1, length(p)))


				logPval = (Lambda^2) / (2 * (Nu + (a*Lambda)/3))
				logPval[Lambda <= 0] = 0
				logPval[is.na(logPval)] = 0
				l = apply(logPval, 1, which.max)-1
			}
		}

		V(sketch.graph)$color = rgb(0.95, 0.95, 0.95)
		V(sketch.graph)$frame.color = rgb(0.9, 0.9, 0.9)

		lg = rgb(0.95, 0.95, 0.95)
		Pal_grad = colorRampPalette(c(lg, Pal[[celltype.name]]))(500)


		V(sketch.graph)$color[l == 1] = darken(scales::col_numeric(Pal_grad, domain = NULL)(scale(x[l == 1])), 0.35)
		V(sketch.graph)$frame.color[l == 1] = darken(V(sketch.graph)$color[l == 1], 0.35)


		if(adjust.node.size == TRUE)
			V(sketch.graph)$size = 3*node.size*x
		else
			V(sketch.graph)$size = node.size
					
		V(sketch.graph)$frame.width = node.size*0.1

		plot(sketch.graph, main = gene)
		
		
		if(!is.na(export_path)) {
			fname = sprintf('%s/%s.pdf', export_path, ifelse(celltype.name == gene, gene, sprintf('%s_%s', celltype.name, gene)));
			pdf(fname)
			plot(sketch.graph, main = gene)
			dev.off()
		}		
	})
}
