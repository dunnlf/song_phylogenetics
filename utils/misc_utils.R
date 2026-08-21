##########################
#
#
# misc. utilities
#
#
##########################

library(RcppAlgos)
library(phytools)
library(ape)
library(cauphy)
library(matrixcalc)
library(PNC)
library(expm)

source('utils/plot.R')

simulate_threshold_model <- function(tree, thresholds, lambdas=NULL) {
  # simulates from (uncorrelated) threshold model with number of dimensions
  # specified by number of thresholds provided
  
  if (is.null(lambdas)) {
    lambdas=rep(1, length(thresholds))
  } else {
    if (length(lambdas) != length(thresholds)) {
      warning('incompatible parameters provided')
      return(NULL)
    }
  }
  
  out_mat <- matrix(NA, nrow=Ntip(tree), ncol=length(thresholds))

  for (i in 1:length(thresholds)){
    lam_i <- lambdas[i]
    th_i <- thresholds[i]
    
    sig <- vcv.phylo(tree)
    diag_sig <- diag(diag(sig))
    sig <- lam_i*(sig-diag_sig) + diag_sig
    eps <- rnorm(Ntip(tree))
    liab <- sqrtm(sig)%*%eps  # liability via transform of std normal
    
    out_mat[,i] <- as.integer(liab>th_i)
  }
  
  return(data.frame(out_mat))
}

compute_unbounded_lambda <- function(tree, x, test=TRUE, estimate=FALSE) {
  # computed unbounded estimate of pagel's lambda via maximum likelihood,
  # if estimate=FALSE, assumes mean of phylogenetic brownian motion is 0, and variance is 1, both known
  
  
  vcv_phylo <- vcv.phylo(tree, model='Brownian')
  vcv_diag <- diag(diag(vcv_phylo))

  
  # likelihood of data given lambda
  lhood <- function(lambda) {
    # define covariance matrix
    sig <- lambda*(vcv_phylo - vcv_diag)
    sig <- sig + vcv_diag
    invSig <- solve(sig)
    n <- ncol(sig)
    
    if(estimate){ # estimate mean and variance
      mu <- as.numeric(sum(invSig%*%x)/sum(invSig))
      sig2 <- as.numeric(t(x-mu)%*%invSig%*%(x-mu)/n)
    } else {
      mu <- 0.0
      sig2 <- 1.0
    }
    
    if (!is.positive.definite(sig)){
      return(-Inf)
    }
    
    # compute likelihood
    return(- t(x - mu)%*%((1/sig2)*invSig)%*%(x-mu)/2 - determinant(sig2 * sig, logarithm = TRUE)$modulus[1]/2 - n*log(2*pi)/2)
  }
  
  H <- nodeHeights(tree)
  max_lambda <- max(H[,2]/max(H[,1]))  # maximum possible scale factor
  
  
  # optimise likelihood over lambda
  # TODO: what is the theoretical minimum lambda?
  res <- optimize(f=lhood, interval=c(-0.05,max_lambda), maximum=TRUE)
  lam_ml <- res$maximum
  
  
  # test via chi-squared test
  if (test){
    stat <- 2*(res$objective - lhood(0))
    p_val <- pchisq(stat, df=1, lower.tail=FALSE, log.p=FALSE)
    
    return(list(lambda=lam_ml,
                logLik = res$objective,
                p=p_val))
  }
  
  return(list(lambda=lam_ml,
              logLik=res$objective))
}


moran_I_permutation_test <- function(tree, X, n_samps=2500, replace=TRUE) {
  # performs a one-sided randomisation test for signifcance of Moran I
  # on phylogenetic data X assumed binary
  
  w <- 1/cophenetic(tree)
  diag(w) <- 0
  
  obs <- Moran.I(X, w)$observed
  
  null_samps <- numeric(n_samps)
  
  for (n in 1:n_samps){
    X_n <- sample(X, length(X), replace=replace)

    if (var(X_n) > 0){
      null_samps[n] <- Moran.I(X_n, w)$observed
    }
    else {
      null_samps[n] <- NA
    }
  }
  
  return(list(observed=obs,
    p=mean(null_samps[!is.na(null_samps)] > obs),
    null_samples=null_samps
  ))
}



plot_lambda_against_moran_I <- function(tree, lambda_range=NULL,
                                           n_points = 400,
                                           n_samps_per_point=250,
                                        quiet=TRUE){
  # compares underlying lambda value in threshold model to induced Moran I value
  
  w <- 1/cophenetic(tree)
  diag(w) <- 0
  
  if (is.null(lambda_range)){
    H <- nodeHeights(tree)
    max_lambda <- max(H[,2]/max(H[,1])) 
    lambda_range <- c(0.0, 1.0)
  }
  lamb_grid <- seq(lambda_range[1],lambda_range[2],
                   (lambda_range[2]-lambda_range[1])/n_points)
  
  out <- numeric(length(lamb_grid))
  
  for (n in 1:length(out)){
    if (!quiet & n%%10==0){
      print(paste('progress:', 100*n/length(out), '/100'))
    }
    
    lamb_n <- lamb_grid[n]
    runs_n <- numeric(n_samps_per_point)
    for (m in 1:n_samps_per_point){
      sig <- vcv.phylo(tree)
      diag_sig <- diag(diag(sig))
      sig <- lamb_n*(sig-diag_sig) + diag_sig
      eps <- rnorm(Ntip(tree))
      liab <- sqrtm(sig)%*%eps  # liability via transform of std normal

      runs_n[m] <- Moran.I(as.integer(liab > 0), w)$observed
    }
    out[n] <- mean(runs_n)
  }
  
  plot(lamb_grid, out, type='l', xlab='Lambda', ylab='I', main="Moran's I, threshold model")
  
  return(list(MoranI=out,
              grid=lamb_grid))
  }



read_mv_ancestral_values <- function(tree, df, motif_inds, raw=FALSE, burnin=0.2,
                                     match_rb_order_to_ape = TRUE){
  # reads multivariate ancestral values
  # returns list of data frames for each motif, if raw: returns continuous values
  
  n <- length(motif_inds)
  m <- ncol(df)
  l <- nrow(df)
  n_taxa <- Ntip(tree)
  
  df <- df[floor(burnin*l):l,]
  
  thresholds <- df[, (m-n+1):m]
  
  df <- df[,5:m]
  
  perm <- match_ape_rb_nodes(tree)
  
  out <- vector('list', n)
  for (i in 1:n) {
    inds <- seq(n*n_taxa + i, n*(2*n_taxa - 3) + i, n)
    
    
    
    ace_i <- df[, inds]
    
    # add root node, which is not recorded and fixed at 0
    ace_i[paste('liability.', 2*n_taxa - 1, '..',i,'.', sep='')] <- rep(0.0, nrow(ace_i))
    
    if (!raw) {
      for (j in 1:nrow(ace_i)) {
        ace_i[j,] <- as.integer(ace_i[j,] > thresholds[j,i])
      
      }
    }
    
    if (match_rb_order_to_ape) {  # reorder to align 
      ace_i <- ace_i[,perm]
    }
    
    out[[i]] <- ace_i
  }
  return(out)
}


read_cor_matrix <- function(filters, runs_df, burnin=0.2, join=TRUE,
                            significance=FALSE){
  # if singifiance==TRUE, returns a matrix of 1/0/-1 for whether the distribution
  # of samples is significantly (0.05) positive (1) negative (-1) or neither (0)
  
  
  keep_inds <- rep(TRUE, nrow(runs_df))
  
  for (col in names(filters)){
    keep_inds <- keep_inds & runs_df[[col]] == filters[[col]]
  }
  
  folders <- (runs_df$save_loc)[keep_inds]
  print(paste('found results at:', folders))
  folders <- sapply(folders, \(x) paste(x, '/rho_entries', sep=''))

  n_entries <- length(list.files(folders[1]))

  n_motifs <- sqrt(2*n_entries + 1/4) + 1/2  # number of motifs
  
  # if significance and join, compute empirical p-value on all sample
  signif_all <- matrix(nrow=n_motifs, ncol=n_motifs)
  diag(signif_all) <- 0
  
  # else, return a list of matrices indicating signficiance for each set of samples
  signif_list <- vector('list', length=length(folders))

  
  out <- vector('list', length=length(folders))
  for(i in 1:length(folders)) {
    out[[i]] <- matrix(nrow=n_motifs, ncol=n_motifs)
    signif_list[[i]] <- matrix(nrow=n_motifs, ncol=n_motifs)
  }
  

  for (j in 1:(n_motifs-1)){
    for (k in (j+1):n_motifs){
      
      all_samps <- numeric(0)
      
      for (i in 1:length(folders)){
        mat_i <- out[[i]]
        signif_mat_i <- signif_list[[i]]
        diag(mat_i) <- rep(1, n_motifs)

        rho_jk_log <- read.csv(paste(folders[i], '/rho',j,k,'.log', sep=''), sep='\t')
        
        m <- nrow(rho_jk_log)
        
        all_samps <- c(all_samps, rho_jk_log[floor(burnin*m):m,5])
          
        rho_jk <- colMeans(rho_jk_log[floor(burnin*m):m,])[5]
          
        mat_i[j,k] <- rho_jk
        mat_i[k,j] <- rho_jk
          
        prop_pos <- mean(rho_jk_log[seq(floor(burnin*m), m, 25),5] > 0)
          
        if (rho_jk > 0){
          signif_mat_i[j,k] <- as.integer((1-prop_pos) < 0.05)
          signif_mat_i[k,j] <- as.integer((1-prop_pos) < 0.05)
        } else {
          signif_mat_i[j,k] <- -as.integer(prop_pos < 0.05)
          signif_mat_i[k,j] <- -as.integer(prop_pos < 0.05)
        }
          
        out[[i]] <- mat_i
        signif_list[[i]] <- signif_mat_i
      }
      
      rho_jk_all <- mean(all_samps)
      prop_pos_all <- mean(all_samps > 0)
      if (rho_jk_all > 0){
        signif_all[j,k] <- as.integer((1-prop_pos_all) < 0.05)
        signif_all[k,j] <- as.integer((1-prop_pos_all) < 0.05)
      } else {
        signif_all[j,k] <- -as.integer(prop_pos_all < 0.05)
        signif_all[k,j] <- -as.integer(prop_pos_all < 0.05)
      }
    }
  }

  if (join & significance){
    return(list(cov=Reduce('+', out)/length(folders),
           signif=signif_all))
  }
  
  if (join){
    return(Reduce('+', out)/length(folders))
  } else {
    if (significance) {
      return(list(cov=out,
                  signif=signif_list))
    } else {
    return(out)
    }
  }
}


get_correlation_from_Mk_model <- function(tree, dat, motif_inds, stationary=TRUE){
  # fits Mk model, and computes correlation implied by stationary distribution
  
  if (!(length(motif_inds) == 2)) {
    warning('Provide 2 motif indeces')
    return(NULL)
  }
  
  if (!('motif_0' %in% colnames(dat))){
    dat <- motif_data_to_binary(tree, dat)[[2]]
  }
  
  if (!stationary){
    table <- table(dat[, motif_inds+1])
    test <- chisq.test(table)
    return(list(cor=sqrt(test$statistic)/nrow(dat),
                p_val=test$p.value))
  } else {
    
    # fit dependent model
    res <- fitPagel(tree, setNames(dat[,motif_inds[1]+1], rownames(dat)),
                    setNames(dat[,motif_inds[2]+1], rownames(dat)))
    
    Q_mat <- matrix(res$dependent.Q, ncol=4, byrow=FALSE)
    
    # estimate stationary distribution
    dist_estimate_1 <- expm(Q_mat*1000)[1,]
    stat_dist <- expm(Q_mat*10000)[1,]
    
    # check convergence
    if (sum(abs(stat_dist - dist_estimate_1)) >1e-8){
      warning('Convergence failure')
      print(paste('1,000', dist_estimate_1))
      print(paste('10,000', stat_dist))
      return(NULL)
    }
    
    stationary_dist_table <- matrix(stat_dist, nrow=2, byrow=TRUE)
    
    corr_exact <- compute_correlation_from_table(stationary_dist_table)
    
    return(list(corr=corr_exact,
                p_val=res$P))
  }
}

compute_correlation_from_table <- function(tab){
  # computes exact correlation from 2x2 table
  mu_1 <- tab[2,1] + tab[2,2]
  mu_2 <- tab[1,2] + tab[2,2]
  
  var_1 <- mu_1*(1 - mu_1)
  var_2 <- mu_2*(1 - mu_2)
  
  tot <- 0.0
  for (i in 0:1) {
    for (j in 0:1) {
      tot <- tot + tab[i+1,j+1]*(i-mu_1)*(j-mu_2)
    }
  }
  return(tot/sqrt(var_1*var_2))
}


get_species_motif_presence <- function(dat){
  # assumes data has columns 'species' and 'gmm_cluster'
  # returns a dataframe of binary motif presence for each species

  species <- unique(dat$species)
  motif_presence <- data.frame(species = species, stringsAsFactors = FALSE)

  for (i in unique(dat$gmm_cluster)){
    motif_presence[[paste('motif_', i, sep='')]] <- as.integer(species %in% dat$species[dat$gmm_cluster == i])
  }

  return(motif_presence)
}















