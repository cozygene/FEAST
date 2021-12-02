
# install.packages("C:/path to folder with the package", repos = NULL, type = "source")
# install.packages("/u/home/u/ulzee/FEAST", repos = NULL, type = "source")

library(ggplot2)
library(ggpubr)
library(glmnet)
library(CVXR)
source('R/FEAST.R')
source('R/STENSL.R')
source('R/utils.R')
source('R/Infer_LatentVariables.R')
Rcpp::sourceCpp('src/rcppSchur.cpp')

MAX_ITERS = 20
COVERAGE_DEPTH = 10000
nEnv = 50

meta = read.table('Data_files/metadata_example_stensl.txt', sep='\t', header=T)
rownames(meta) = meta$SampleID
otus = read.table('Data_files/otu_example_stensl.txt', sep='\t', header=T)
print(meta)
# print(dim(otus))
# print(rownames(otus)[1:10])
# print(meta$SampleID[1:10])

feast.result = FEAST(
	C=as.matrix(otus),
	metadata=meta,
	EM_iterations=MAX_ITERS,
	COVERAGE=COVERAGE_DEPTH,
	different_sources_flag=0,
)

stensl.result <- STENSL(
	C=as.matrix(otus),
	metadata=meta,
	EM_iterations=MAX_ITERS,
	COVERAGE=COVERAGE_DEPTH,
	l.range=c(0.1,1,10)
)

ps = unlist(stensl.result$proportions_mat)
g1 = ggplot(data.frame(
	p=c(sum(ps[1:10]), sum(ps[11:50]), ps[51]),
	label=c('Src', 'FP', 'Unk')),
	aes(x = "", y = p, fill = label)
	) +
	geom_col() +
	coord_polar(theta = "y") +
	ggtitle('STENSL') +
	geom_hline(yintercept=0.4,linetype = 'dotted', col = 'blue')+
	geom_hline(yintercept=0,linetype = 'dotted', col = 'blue')

ps = unlist(feast.result$proportions_mat)
g2 = ggplot(data.frame(
	p=c(sum(ps[1:10]), sum(ps[11:50]), ps[51]),
	label=c('Src', 'FP', 'Unk')),
	aes(x = "", y = p, fill = label)
	) +
	geom_col() +
	coord_polar(theta = "y") +
	ggtitle('FEAST')+
	geom_hline(yintercept=0.4,linetype = 'dotted', col = 'blue')+
	geom_hline(yintercept=0,linetype = 'dotted', col = 'blue')

g = ggarrange(g1, g2, common.legend=T, legend='bottom')
g = annotate_figure(g, top = text_grob('Simulated unknown: 40%'))
ggsave('vignettes/comparison.jpg', g, width=8, height=5)
