# Code for phylogenetic analysis of passerine songs

Tree and data from: https://github.com/quentinbacquele/phylogenetic_analysis


Code for analysing binary trait co-occurence using the threshold model in RevBayes.

### Data:

Trait data comprise motif presence/absence for the 8 Motifs identified in Bacquelé, Q. et al. (2026). The global biogeography of passerine songs. Science 393(6809), 381–386. https://doi.org/10.1126/science.aee6239

**Oscines** - 100 species thinned from QB tree
- max_genera_100_osc_dat.csv
- max_genera_100_tree.tree

**Suboscines** - 133 species thinned from QB tree
- suboscine_uncapped_thinned.csv
- tree_suboscine_thinned.tree

### Scripts

`rb-scripts` contains RevBayes scripts for running the threshold model, either to draw posterior samples via MCMC, or to compute the marginal likelihood via stepping stone sampling.

`output/` contains outputs of RevBayes, as well as `runs.csv` which contains metadata and file locations for previous runs.

`scripts` contains R wrappers for running the RevBayes scripts on given data.

Finally, scripts in `scripts/analyses` produce the core figures for the two data sets, reading from `runs.csv` for MCMC results.




