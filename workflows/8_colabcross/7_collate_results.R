#!/home/momaksimov/anaconda3/envs/r4/bin/Rscript

# about: collate results from analysis using different minprob thresholds

# options
options(stringsAsFactors = FALSE, dplyr.summarise.inform = FALSE)

# libraries
library(tidyverse)
library(cowplot)
library(fs)

# config
data_dir = '../../data/colab_cross/denovo_info'

# load denovo genotypes
denovo_gts = map_df(dir_ls(data_dir, regexp = '\\d+_.*.rds'), readRDS)

# check completeness
denovo_gts %>% pull(chr) %>% unique
all(str_c('chr', 1:19) %in% (denovo_gts %>% pull(chr) %>% unique))

# collate into tibble
denovo_gts = denovo_gts %>% map_df(~.x, .id = 'minprob')

# check number of denovo loci
denovo_gts %>% distinct(chr, pos, end) %>% nrow

# save combined data
saveRDS(denovo_gts, path(data_dir, 'denovo_gts.rds'))
saveRDS(denovo_gts %>% distinct(chr, pos, end), path(data_dir, 'denovo_loci.rds'))
