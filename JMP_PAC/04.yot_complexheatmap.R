#!/usr/bin/env Rscript


library(ComplexHeatmap)

# mat <- read.table('~/Downloads/lihc2.tsv',header = T, stringsAsFactors = F, sep='\t')
setwd('~/Dropbox/NCI/myc signature/oncoprint')
mat <- read.table('lihc_27genes.tsv',header = T, stringsAsFactors = F, sep='\t')

mat[,1:2] <- NULL
mat[,2] <- NULL

colnames(mat) <- gsub('..AMP.HOMDEL.MUT.FUSION.','',colnames(mat))
rownames(mat) <- mat[,1]
mat[,1] <- NULL

for (i in 1:length(mat[1,])){	
	mat[,i] <- gsub('CNA: AMP;','AMP',mat[,i])	
	mat[,i] <- gsub('MUT\\:.*;','MUT;',mat[,i])
	mat[,i] <- gsub('FUSION\\:.*;','FUSION;',mat[,i])
	mat[,i] <- gsub('CNA: HOMDEL;','HOMDEL;',mat[,i])
		
}

for (i in 1:length(mat[1,])){
	mat[,i] <- gsub('AMP','AMP;',mat[,i])
	mat[,i] <- gsub(' ','',mat[,i])
}

matt <- t(as.matrix(mat))

table(matt)

# rem <- c('AGO2','NELFE')
matt <- matt[!rownames(matt) %in% rem,] 

alter_fun <- list(
	background = function(x,y,w,h) {
		grid.rect(x, y, w-unit(0.5, 'mm'), h-unit(0.5, 'mm'), gp = gpar(fill = '#cccccc', col = NA))
	},
	HOMDEL = function(x,y,w,h) {
		grid.rect(x, y, w-unit(0.5, 'mm'), h-unit(0.5, 'mm'), gp = gpar(fill = 'blue', col = NA))
	},
	AMP = function(x,y,w,h) {
		grid.rect(x, y, w-unit(0.5, 'mm'), h-unit(0.5, 'mm'), gp = gpar(fill = 'red', col = NA))
	},
	MUT = function(x,y,w,h) {
		grid.rect(x, y, w-unit(0.5, 'mm'), h*0.33, gp = gpar(fill = '#008000', col = NA))
	},
	FUSION = function(x,y,w,h) {
		grid.rect(x, y, w-unit(0.5, 'mm'), h*0.50, gp = gpar(fill = 'yellow', col = NA))
	}
	
)

col = c('MUT' = '#008000', 'AMP' = 'red', 'HOMDEL' = 'blue', 'FUSION' = 'yellow')


# ht_list <- oncoPrint(	matt, 
			# get_type = function(x) strsplit(x, ';')[[1]],
			# alter_fun = alter_fun, 
			# col = col, 
			# column_title = 'TCGA LIHC 27 genes',
			# heatmap_legend_param = list(title = 'Alternations', at = c('AMP', 'HOMDEL', 'MUT','FUSION'), labels = c('Amplification', 'Deep deletion', 'Mutation','Fusion')),
			# split = sample(letters[1:2], nrow(matt), replace = T))
			
# draw(ht_list, row_sub_title_side = 'left')

# setwd('~/Dropbox/NCI/myc signature/oncoprint')
# pdf(file = 'oncoprint_27genes.pdf', width = 25, height = 7)

# oncoPrint(matt, 
			# get_type = function(x) strsplit(x, ';')[[1]],
			# alter_fun = alter_fun, 
			# col = col, 
			# column_title = 'TCGA LIHC 27 genes',
			# heatmap_legend_param = list(title = 'Alternations', at = c('AMP', 'HOMDEL', 'MUT','FUSION'), labels = c('Amplification', 'Deep deletion', 'Mutation','Fusion')))

# dev.off()

source('~/Dropbox/NCI/myc signature/prognostic.index/prognostic.index.calculation.R', chdir = TRUE)

sample.hi <- rownames(sample.score)[which(sample.score[,2] == 'high')]
sample.lo <- rownames(sample.score)[which(sample.score[,2] == 'low')]

all.sample <- c(sample.hi,sample.lo)

common.sample <- intersect(all.sample,colnames(matt))
common.hi <- intersect(sample.hi,colnames(matt))
common.lo <- intersect(sample.lo,colnames(matt))

matt.common <- matt[,common.sample]
matt.hi <- matt[,common.hi]
matt.lo <- matt[,common.lo]

setwd('~/Dropbox/NCI/myc signature/oncoprint')
pdf(file = 'oncoprint_27genes.hi.pdf', width = 9, height = 7)
onco.hi <- oncoPrint(matt.hi,
			get_type = function(x) strsplit(x, ';')[[1]],
			alter_fun = alter_fun, 
			col = col, 
			column_title = 'TCGA LIHC 27 genes of NDMT',
			heatmap_legend_param = list(title = 'Alternations', at = c('AMP', 'HOMDEL', 'MUT','FUSION'), labels = c('Amplification', 'Deep deletion', 'Mutation','Fusion')))
onco.hi
dev.off()

hi.onco.order <- colnames(matt.hi)[column_order(onco.hi)[[2]]]

pdf(file = 'oncoprint_27genes.lo.pdf', width = 25, height = 7)
onco.lo <- oncoPrint(matt.lo, 
			get_type = function(x) strsplit(x, ';')[[1]],
			alter_fun = alter_fun, 
			col = col, 
			column_title = 'TCGA LIHC 27 genes of Others',
			heatmap_legend_param = list(title = 'Alternations', at = c('AMP', 'HOMDEL', 'MUT','FUSION'), labels = c('Amplification', 'Deep deletion', 'Mutation','Fusion')))
onco.lo
dev.off()

lo.onco.order <- colnames(matt.lo)[column_order(onco.lo)[[2]]]

sample.order <- c(hi.onco.order,lo.onco.order)

pdf(file = 'oncoprint_27genes.all.pdf', width = 25, height = 7)
oncoPrint(matt[,sample.order], 
			get_type = function(x) strsplit(x, ';')[[1]],
			column_order = NULL,
			alter_fun = alter_fun, 
			col = col, 
			column_title = 'TCGA LIHC 27 genes of All samples\nNDMT first, then Others',
			heatmap_legend_param = list(title = 'Alternations', at = c('AMP', 'HOMDEL', 'MUT','FUSION'), labels = c('Amplification', 'Deep deletion', 'Mutation','Fusion')))
dev.off()