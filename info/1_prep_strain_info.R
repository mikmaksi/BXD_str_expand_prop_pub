#!/home/momaksimov/anaconda3/envs/r4.2/bin/Rscript

# about: prep information about strains for convenient loading in other scripts

# options
options(stringsAsFactors = FALSE)

# libraries
library(tidyverse)
library(cowplot)
library(fs)

# load strain info
strain_info = read_tsv('bxd_strain_names_plus.tsv',
    col_types = cols(
      long_name      = 'c',
      short_name     = 'c',
      bxd_id         = 'c',
      batch          = 'c',
      off_epoch      = 'c',
      type           = 'c',
      gen_inbreeding = 'i',
      gen_breeding   = 'i',
      note           = 'c',
      cross_type     = 'c'
    )
)

# remove redundant columns
strain_info = strain_info %>% select(!c(batch, note, type))

# get list of strains for which sequencing data is available (either have strain in snp vcf or bams for str genotyping)
seqed_strains = list(snp = 'snp_strain_list', 
		     str = 'str_strain_list') %>%
    map_df(~read_tsv(.x, 
		     col_names = 'short_name', 
		     col_types = cols('c')), .id = 'seq_data') %>%
    mutate(is_seq = TRUE) %>%
    pivot_wider(id_cols = short_name, 
		names_from = seq_data, names_prefix = 'is_seq_',
		values_from = is_seq, values_fill = list(is_seq = FALSE))
strain_info = strain_info %>% 
    left_join(seqed_strains, by = 'short_name') %>%
    mutate(across(matches('is_seq_'), ~replace_na(.x, FALSE)))

# save object
saveRDS(strain_info, 'strain_info.rds')

