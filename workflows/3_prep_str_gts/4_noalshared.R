#!/home/momaksimov/anaconda3/envs/r4.2/bin/Rscript

# about: find str loci where parent-of-origin can be found (no alleles shared) 

# options
options(stringsAsFactors = FALSE)

# libraries
library(tidyverse)
library(fs)
library(DBI)
library(dbplyr)
library(progress)

# for parallel
# library(furrr)
# library(progressr)

# config
out_dir = '../../data/str_gts'; dir_create(out_dir)
geno_file = 'all_repcn_proc_nosegdup_nolowcr.rds'

# read processed version of all the genotypes
gen = readRDS(path(out_dir, geno_file))

# find loci where founder strains don't share and allele and are not missing
no_alshared = geno %>% 
    select(chr, pos, end, C57BL, DBA) %>%
    rename(B = C57BL, D = DBA) %>%
    filter(!is.na(B)) %>%
    filter(!is.na(D)) %>%
    filter(!B==D) %>%
    mutate(across(c(B, D), ~.x %>% str_split(','))) %>%
    mutate(cond = map2_lgl(B, D, function(a, b) {a[1] == b[1] || a[1] == b[2] || a[2] == b[1] || a[2] == b[2]})) %>%
    filter(!cond)
   
# write locus list
write_tsv(no_alshared %>% select(chr, pos, end), '../../data/ref/str_noalshared.tsv')

{ # DEBUGGING
    # get previous lists
    data_dir = '/projects/ps-gymreklab/mikhail/090520_unified_workflow/data/ref'
    ori_regions = map_df(
	list(lowcr = 'exclude/lowcr_strs.bed',
	     segdup = 'exclude/segdup_strs.bed',
	     # zv_hetonly = 'exclude/zv_hetonly_strs.bed', # don't need zv het only
	     no_alshared = 'str_regions_founder_no_shared_alleles.bed'),
	~read_tsv(path(data_dir, .x), 
		  col_names = c('chr', 'pos', 'end'), 
		  col_types = cols(chr = 'c', pos = 'i', end = 'i')), 
	.id = 'type')
    ori_noalshared = ori_regions %>%
	filter(type == 'no_alshared') %>%
	anti_join(ori_regions %>% filter(type != 'no_alshared'), by = c('chr', 'pos', 'end'))
    
    # get new lowcr
    lowcr_loci = read_tsv('../../data/ref/str_lowcr.tsv')

    # compare counts
    comp = list(ori = ori_noalshared %>% select(chr, pos, end), new = no_alshared %>% select(chr, pos, end)); comp %>% map(nrow)
    
    # no diff once new lowcr loci are accounted for
    miss = setdiff(comp$ori, comp$new) %>% anti_join(lowcr_loci)

    # some extra loci
    extra = setdiff(comp$new, comp$ori)
}

# subset list of genotypes
no_alshared = geno %>% semi_join(no_alshared, by = c('chr', 'pos', 'end'))

# convert to B/D labels
no_alshared_recode = no_alshared %>% 
    mutate(across(matches('BXD'), ~case_when(.x == C57BL ~ 'B', .x == DBA ~ 'D', TRUE ~ NA_character_))) %>%
    select(!c(C57BL, DBA))

# check counts
no_alshared_recode %>% 
    pivot_longer(matches('BXD'), names_to = 'strain', values_to = 'gt') %>%
    count(strain, gt) %>%
    pivot_wider(id_cols = strain, names_from = gt, values_from = n)

# save genotype lists
# write_tsv(no_alshared, path(out_dir, 'noalshared_repcn.tsv.gz'))
saveRDS(no_alshared, path(out_dir, 'all_repcn_proc_nosegdup_nolowcr_noalshared.rds'))
saveRDS(no_alshared_recode, path(out_dir, 'all_foulab_nosegdup_nolowcr_noalshared.rds'))

# load list of segregating strs
segreg = read_tsv('../../data/ref/str_segreg.tsv', col_types = cols(chr = 'c', pos = 'i', end = 'i'))

# padd the no_alsharead founder genotypes with segreg loci
padded = no_alshared_recode %>%
    bind_rows(segreg %>% anti_join(no_alshared_recode, by = c('chr', 'pos', 'end'))) %>%
    mutate(chr = fct_relevel(chr, str_c('chr', 1:19))) %>%
    arrange(chr, pos, end)

# write the padded founder genotypes
saveRDS(padded, path(out_dir, 'all_foulab_nosegdup_nolowcr_noalshared_padded.rds'))
