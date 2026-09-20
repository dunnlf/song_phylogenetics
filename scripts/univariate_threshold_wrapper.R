##########################
#
#
# run univariate MCMC / marginal likelihood via RevBayes
#
#
##########################

library(ape)
library(stringr)

# read data and tree save to tmp files, choose motif to analyse
source('utils/data_handler.R')
data_path <- 'data/oscine_data.csv'
tree_path <- 'data/oscine_tree.tree'

dat <- read.csv(data_path)
tree <- read.tree(tree_path)

motif_inds <- 2:7
null_dist <- FALSE # whether or not to shuffle tip values

for (motif_ind in rep(motif_inds, 1)){
  for (analysis in c('MCMC')){
  
  res <- motif_data_to_binary(tree, dat)
  
  tree <- res[[1]]
  motif_data <- res[[2]]
  
  if (null_dist){
    # resample data and assign to original species
    motif_data <- setNames(motif_data[sample(1:nrow(motif_data), nrow(motif_data),
                                    replace=TRUE),], tree$tip.label)
  }
  
  write.tree(tree, 'data/tree_tmp.tree')
  write_to_nexus(motif_data, motif_ind)
  
  
  # RevBayes script paths
  model_scripts = data.frame(row.names=c('Mk', 'threshold', 'HRM', 'threshold_alt'),
                       loc=c('rb-scripts/lines/Mk.Rev',
                               'rb-scripts/lines/threshold.Rev',
                               'rb-scripts/lines/HRM.Rev',
                               'rb-scripts/lines/threshold_altmodel.Rev'))
  
  analysis_scripts = data.frame(row.names=c('MCMC', 'ML'),
                                loc=c('rb-scripts/lines/MCMC.Rev',
                                      'rb-scripts/lines/ML.Rev'))
  
  # parameters
  model_type = 'threshold_alt' # available: Mk, threshold, threshold_alt
  n_chains = 1
  n_generations = 8e5
  burnin_generations = 300000 # for ML only
  tuning_interval = 1000 # for ML only
  printgen=500
  sd_prior=10.0
  run_id = NULL ### populate with a run id, else date and time used
  
  if (is.null(run_id)) {
      run_id = date()
  }
  
  # folder of run results
  save_loc = paste('output/',
                   model_type,
                   '/',
                   analysis,
                   '_',
                   run_id,
                   sep=''
  )
  save_loc <- str_replace_all(save_loc, ':', '')
  
  # temporary Rev script
  tmp <- tempfile(fileext = '.Rev')
  
  con <- file(tmp, open='a')
  
  # write parameter values
  lines_0 <- c(
    paste('save_loc ="', save_loc, '"', sep=''),
    paste('n_chains =', n_chains, sep=''),
    paste('burnin_generations =', burnin_generations, sep=''),
    paste('n_generations =', n_generations, sep=''),
    paste('tuning_interval =', tuning_interval, sep=''),
    paste('printgen =', printgen, sep='')
  )
  
  writeLines(lines_0, con)
  
  # add model code
  lines_model <- readLines(model_scripts[model_type,])
  if (model_type %in% c('threshold_alt', 'threshold')){
    lines_model[11] <- str_replace(lines_model[11], 'SD_PRIOR', as.character(sd_prior))
  } else {
    sd_prior <- NULL
  }
  
  writeLines(lines_model, con)
  
  # add code to run analysis
  lines_analysis <- readLines(analysis_scripts[analysis,])
  writeLines(lines_analysis, con)
  
  
  close(con)
  
  
  # run script
  system(paste('rb', tmp))
  
  # writing details of run to results list
  results <- read.csv('output/runs.csv', row.names = 1, header=TRUE)
  run_details <- c(date(),
                   model_type,
                   analysis,
                   run_id,
                   n_chains,
                   n_generations,
                   printgen,
                   tuning_interval,
                   save_loc,
                   data_path,
                   tree_path,
                   null_dist,
                   as.character(sd_prior),
                   motif_ind)
  
  results[nrow(results)+1,] <- run_details
  write.csv(results, file='output/runs.csv')
  
  # remove script
  unlink(tmp)
}
}

