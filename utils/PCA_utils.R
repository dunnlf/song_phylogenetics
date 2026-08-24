##########################
#
#
# PCA reconstruction utilities
#
#
##########################


get_PCA_distance_matrix <- function(PCA_df, motif_inds, average_within_species=TRUE){
  # assumed PCA_df has final column for gmm cluster, and remaining are PCA components
  # if average_within_species is TRUE, then distances are computed within each species, then averaged
  # else distance taken between centroids of each motif within species having both motifs
  
  species_motifs <- get_species_motif_presence(PCA_df)

  dist_PCA <- matrix(NA, nrow=length(motif_inds), ncol=length(motif_inds))
  for (i in 1:length(motif_inds)){
    for (j in 1:length(motif_inds)){
      if (i == j){
        dist_PCA[i,j] <- 0
      } else {
        # subset to only species having both motif i and motif j
        species_ij = species_motifs[species_motifs[, paste('motif_', motif_inds[i], sep='')] == 1 &
                                      species_motifs[, paste('motif_', motif_inds[j], sep='')] == 1,]
        species_ij = species_ij$species

        PCA_ij <- PCA_df[PCA_df$species %in% species_ij &
                           PCA_df$gmm_cluster %in% c(motif_inds[i], motif_inds[j]),]
        if (average_within_species){
          dists_ij <- c()
          for (s in species_ij){
            PCA_s <- PCA_ij[PCA_ij$species == s,]
            PCA_means_s <- PCA_s[, 2:ncol(PCA_s)] %>%
              group_by(gmm_cluster) %>%
              summarise_all(mean)

            dist_ij <- dist(PCA_means_s[,2:ncol(PCA_means_s)])

            dists_ij <- c(dists_ij, dist_ij)
          }

          dist_PCA[i,j] <- mean(dists_ij)
          
        } else {
          PCA_means <- PCA_ij[, 2:ncol(PCA_ij)] %>%
            group_by(gmm_cluster) %>%
            summarise_all(mean)

          dist_ij <- dist(PCA_means[,2:ncol(PCA_means)])

          dist_PCA[i,j] <- dist_ij
        }
      }
    }
  }
  
  return(dist_PCA)
}


plot_PCA_centroids <- function(PCA_df, motif_inds,
                               main='Distribution of Motifs', cols=NULL,
                               alpha=0.05, names=NULL){
  
  kdes = list()
  
  if (is.null(cols)){
    cols <- rainbow(length(motif_inds))
  }
  
  PCA_means <- PCA_df[,c(1:(ncol(PCA_df)))] %>%
    group_by(gmm_cluster) %>%
    summarise_all(mean)
  print(PCA_means)
  
  for (i in 1:length(motif_inds)){
    ind = motif_inds[i]
    
    PCA_ind <- PCA_df[PCA_df$gmm_cluster == ind, ]
    
    PCA_ind <- PCA_ind[, 1:2]

    kdes[[i]] <- MASS::kde2d(PCA_ind$PC1, PCA_ind$PC2)
  }
  
  xlim <- range(unlist(lapply(kdes, `[[`, "x")))
  ylim <- range(unlist(lapply(kdes, `[[`, "y")))
  
  contour(kdes[[1]]$x, kdes[[1]]$y, kdes[[1]]$z, col=cols[1],
          xlim = xlim, ylim = ylim, cex.axis=1.3, cex.lab = 1.3,
          xlab='PC1', ylab='PC2')
  PCA_ind <- PCA_df[PCA_df$gmm_cluster == motif_inds[1], ]
  points(PCA_ind$PC1, PCA_ind$PC2,bg=alpha(cols[[1]], alpha), col=alpha(cols[[1]], alpha), pch=21)
  
  for (i in 2:length(kdes)){
    
    contour(kdes[[i]]$x, kdes[[i]]$y, kdes[[i]]$z, col=cols[i],
            xlim = xlim, ylim = ylim, add=TRUE)
    
    PCA_ind <- PCA_df[PCA_df$gmm_cluster == motif_inds[i], ]
    points(PCA_ind$PC1, PCA_ind$PC2,bg=alpha(cols[[i]], alpha), pch=21, col=alpha(cols[[i]], alpha))
    
  }
  
  
  # add centroids
  for (i in 1:length(kdes)) {
    points(PCA_means[PCA_means$gmm_cluster == motif_inds[i],]$PC1,
           PCA_means[PCA_means$gmm_cluster == motif_inds[i],]$PC2, pch=23,
           bg=cols[[i]], col='black', cex=3)
  }
  
  legend("topright", legend = names,
         col = cols, lty = 1, title='Motif', cex=1.5)
  title(main, cex.main=1.85)
}


