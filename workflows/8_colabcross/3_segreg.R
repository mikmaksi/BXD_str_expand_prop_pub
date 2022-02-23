#!/home/momaksimov/anaconda3/envs/r4.2/bin/Rscript

# about: find segregating str loci - not monoallelic and not only 

# options
options(stringsAsFactors = FALSE)

# libraries
library(tidyverse)
library(fs)
library(DBI)
library(dbplyr)
library(progress)

# config
out_dir   = '../../data/colab_cross/str_gts'; dir_create(out_dir)
ref_dir   = '../../data/colab_cross/ref'; dir_create(ref_dir)
geno_file = 'all_repcn_proc_nosegdup_nolowcr.rds'

# read processed version of all the genotypes
geno = readRDS(path(out_dir, geno_file))

# blank out hets (this is fast column-wise)
# NAs are already blanked out
geno_het_clean = geno %>%
    mutate(across(c(matches('^CC'), LETTERS[1:8]), function(.x) {
	.y = .x %>% str_split(',', simplify = TRUE)
	.x[.y[,1] != .y[,2]] = NA
	.x
    }))

# count number of unique genotypes per row
# CHECK: CC_variants_MM_strs_highQ_calls_strain_filt.vcf.gz with 74 strains (number of strains checks out)
# CHECK: founders are counted in number of unique genotypes

cache_file = path(out_dir, 'unq_gt_per_loc.tsv'); redo = FALSE
if (!file_exists(cache_file) | redo) {
    pb <- progress_bar$new(total = dim(geno_het_clean)[1])
    n_unq_gt = geno_het_clean %>% 
	select(matches('^CC'), LETTERS[1:8]) %>%
	apply(MARGIN = 1, FUN = function(.x) {pb$tick(); .x %>% discard(is.na) %>% unique %>% length})
    unq_gt_per_loc = geno %>% select(chr, pos, end) %>% mutate(n_unq_gt = !!n_unq_gt)

    # save unq_gt_per_loc
    write_tsv(unq_gt_per_loc, cache_file)
} else { unq_gt_per_loc = read_tsv(cache_file, col_types = cols(chr = 'c', pos = 'i', end = 'i', n_unq_gt = 'i')) }

# spot check some of these
if (0) {
    geno_het_clean %>%
	semi_join(unq_gt_per_loc %>% filter(n_unq_gt == 2) %>% slice_sample(n = 1)) %>%
	pivot_longer(matches('^CC'), 'var', 'val') %>%
	count(value)
}

# make list of segregating strs
segreg = geno %>%
    semi_join(unq_gt_per_loc %>% filter(n_unq_gt >= 2), by = c("chr", "pos", "end"))

# make a list of monoallelic loci
monoal_loci = unq_gt_per_loc %>% filter(n_unq_gt < 2) %>% select(chr, pos, end)

# sort list by locus
segreg = segreg %>% 
    mutate(chr = fct_relevel(chr, str_c('chr', 1:19))) %>%
    arrange(chr, pos, end)

# write segregating genotypes and locus list
saveRDS(segreg, path(out_dir, 'all_repcn_proc_nosegdup_nolowcr_segreg.rds'))
write_tsv(segreg %>% select(chr, pos, end), path(ref_dir, 'str_segreg.tsv'))
write_tsv(monoal_loci %>% select(chr, pos, end), path(ref_dir, 'str_zv_hetonly.tsv'))
