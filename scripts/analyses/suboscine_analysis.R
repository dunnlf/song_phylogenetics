##########################
#
#
# script for producing suboscine results and figures
#
#
##########################

library(pheatmap)
library(grid)
library(gridExtra)
library(corrplot)

motif_names <- c("Flat whistles",
                 "Slow trills",
                  "Fast trills",
                  "Chaotic songs",
                  "Ultrafast trills",
                  "Slow mod. songs",
                  "Fast mod. songs",
                  "Harmonic stacks")


# load utility functions
utils <- list.files('utils')
for (file in utils) {
  if (!grepl('.Rev', file)){
    source(paste('utils/', file, sep=''))
  }
}

# load data and tree
data_path <- 'data/final/suboscine_uncapped_thinned.csv'
tree_path <- 'data/final/tree_suboscine_thinned.tree'

data <- read.csv(data_path)
tree <- read.tree(tree_path)

# convert gmm clusters binary motif presence, and align rows with tree tip labels
motif_data <- motif_data_to_binary(tree, data)[[2]]


### 1 - identify significant motifs via a permutation test with Moran's I
set.seed(2026)

p_vals <- numeric(8)
moran_res_list <- vector('list', 8)

for (i in 0:7) {
  moran_res_i <- moran_I_permutation_test(tree, motif_data[, i+1], n_samps=10000)
  moran_res_list[[i+1]] <- moran_res_i
  print(paste('Motif ', i, '  I: ', round(moran_res_i[[1]],3),
              '   p: ', round(moran_res_i[[2]],3), sep=''))
  p_vals[i+1] <- moran_res_i[[2]]
}

# visualise selected null distributions
vis_motifs <- c(4,0)
for (ind in vis_motifs){
  moran_res <- moran_res_list[[ind+1]]
  
  hist(moran_res[[3]], main=(paste("Null distribution of Moran's I, motif", ind)),
       xlab=substitute(italic('I')), breaks=30, col='lightblue')
  abline(v=moran_res[[1]], col='red', lty='dotted', lwd=3.5,
         cex.axis=3, cex.main=3)
  legend('topright', legend='Observed', lty='dotted', col='red', lwd=3.5)
}


# stable motifs
significant_motifs <- which(p_vals<0.05) - 1


### plotting recordings of motifs in PCA space

pca_inds <- 17:(17+37)  # columns of PCA data and gmm cluster
plot_PCA_centroids(data[,pca_inds], significant_motifs, main='Suboscine Motifs')


### plotting correlation estimates and PCA distance matrices

# reading log of RevBayes
runs <- read.csv('output/runs.csv', colClasses = c(motif="character"))

cov_mat <- read_cor_matrix(list(model_type='MV_threshold', motif=paste(significant_motifs, sep='', collapse=''), data=data_path), runs, join=TRUE)
pca_mat <- get_PCA_distance_matrix(data[,c(4, pca_inds)], significant_motifs)

rownames(cov_mat) <- significant_motifs
colnames(cov_mat) <- significant_motifs

rownames(pca_mat) <- significant_motifs
colnames(pca_mat) <- significant_motifs


# plot cov and pca matrices via corrplot
signif_motif_names <-  motif_names[significant_motifs+1]
rownames(cov_mat) <- signif_motif_names
colnames(cov_mat) <- signif_motif_names
rownames(pca_mat) <- signif_motif_names
colnames(pca_mat) <- signif_motif_names
corrplot(cov_mat, method='circle', type='upper', diag=FALSE, addCoef.col = 'black',
         number.cex=1.75, cl.cex=1.2, tl.cex=1.75, tl.col='black', tl.srt=30,
         title='Correlation estimates', mar=c(1,1,3.5,1), cex.main=2, col=COL2("RdBu", 200))

corrplot(pca_mat, is.corr=FALSE, method='circle', type='upper', diag=FALSE, addCoef.col = 'black',
         number.cex=1.75, cl.cex=1.2, tl.cex=1.75, tl.col='black', tl.srt=30,
         col=colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu"))[c(2,3,5,9)])(200), col.lim=c(0,70),
         title='PCA distances', mar=c(1,1,3.5,1), cex.main=2)


# confirm results with correlation and significance from Mk model of Pagel 1994

n_motifs <- length(significant_motifs)
signif_Mk <- matrix(nrow=n_motifs,
                    ncol=n_motifs)
corr_Mk <- matrix(nrow=n_motifs,
                  ncol=n_motifs)
diag(corr_Mk) <- 1.0

for (i in 1:(n_motifs-1)) {
  for (j in (i+1):n_motifs){
    res <- get_correlation_from_Mk_model(tree, data,
                                         motif_inds=significant_motifs[c(i,j)])
    corr_ij <- res[[1]]
    signif_ij <- as.integer(res[[2]] < 0.05) * sign(corr_ij)
    
    signif_Mk[i,j] <- signif_ij
    signif_Mk[j,i] <- signif_ij
    
    corr_Mk[i,j] <- corr_ij
    corr_Mk[j,i] <- corr_ij
    
  }
}

print(corr_Mk)
print(signif_Mk)
