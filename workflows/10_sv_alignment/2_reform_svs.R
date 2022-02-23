#!/home/momaksimov/anaconda3/envs/r4.2/bin/Rscript

# about: make a new SV data file in a simple format where start,stop,st_type are indicated
# along with strains that contain this variant
# use REF/ALT alignments from previous script

# clean variables
rm(list = ls())

# options
options(stringsAsFactors = FALSE, dplyr.summarise.inform = FALSE)

# libraries
library(tidyverse)
library(cowplot)
library(fs)

# config
vcf_dir = '/projects/ps-gymreklab/resources/datasets/BXD/pangenome/symlinks'
# vcf_file = 'chr13.pang.5mb.large.vcf.gz'
vcf_file = 'chr13.pang.10mb.large.vcf.gz'
alig_files_dir = '../../data/sv_alig/refalt_alig'
out_dir = '../../data/sv_alig'; dir_create(out_dir)

# load REF/ALT alignment data
alig_data = dir_ls(alig_files_dir) %>% map_df(readRDS)

# check that sv_size matches calculated total_sv
# alig_data %>% count(sv_size == total_sv)

# load SV genotypes
cmd = sprintf("bcftools query -H -f '%%CHROM\t%%POS\t%%REF\t%%ALT[\t%%GT]\n' %s", path(vcf_dir, vcf_file))
gt_data = read_tsv(pipe(cmd), col_types = cols(.default = 'c'))

# rename columns
gt_data = gt_data %>% 
    rename_all(~str_replace(.x, '\\[.*\\]', '') %>% str_replace('# ', '') %>% str_replace(':GT', '')) %>%
    mutate(across(c(POS), as.integer)) %>%
    rename(chr = CHROM, pos = POS)

# remove duplicates if any
gt_data = gt_data %>% distinct

# check for overlap of loci b/w alig_data and gt_data
dkeys = list(
    alig_data = alig_data %>% distinct(chr, pos, REF, ALT),
    gt_data = gt_data %>% distinct(chr, pos, REF, ALT)
)

# there are some for which not alignment data
# NOTE: need to investigate, for now just filter out
setdiff(dkeys$alig_data, dkeys$gt_data)
setdiff(dkeys$gt_data, dkeys$alig_data)
dkeys = intersect(dkeys$gt_data, dkeys$alig_data)

# add a unique identifier
dkeys = dkeys %>%
    unite('uid', c('chr', 'pos', 'REF', 'ALT'), remove = FALSE) %>%
    mutate(uid = uid %>% as.factor %>% as.integer)
join_keys = function(data, dkeys) {
    data %>% inner_join(dkeys, by = c('chr', 'pos', 'REF', 'ALT')) %>% select(!c(chr, pos, REF, ALT)) %>% relocate(uid, .before = 1)
}
alig_data = alig_data %>% join_keys(dkeys)
gt_data = gt_data %>% join_keys(dkeys)

# recalculate largest_sv (was wrong in previous script)
# NOTE: this is simply a check
if (0) {
	largest_sv_per_alig = alig_data %>%
		select(uid, alig_tbls) %>%
		mutate(indels = map(alig_tbls, ~.x$indels)) %>%
		select(!alig_tbls) %>%
		unnest(indels) %>% 
		arrange(desc(width)) %>%
		distinct(uid, .keep_all = TRUE) %>%
		select(uid, largest_sv = width)
		alig_data = alig_data %>%
		select(!largest_sv) %>%
		left_join(largest_sv_per_alig, by = 'uid')
}

# there are multiple rows per uid with different number of missing genotypes
# this shouldn't happen
# calculate the number of missing genotypes per locus and keep one with fewest missing values
gt_data = gt_data %>%
    rowwise %>%
    mutate(n_miss = sum(c_across(!uid) == '.')) %>%
    ungroup %>%
    arrange(n_miss) %>%
    distinct(uid, .keep_all = TRUE) %>%
    select(!n_miss)

# format long
gt_data_long = gt_data %>%
    pivot_longer(!uid, names_to = 'strain', values_to = 'gt')

# calculate unique genotypes per locus 
ngt_per_loc = gt_data %>%
    pivot_longer(!uid, names_to = 'strain', values_to = 'gt') %>%
    group_by(uid) %>%
    summarise(n_unq_gt = gt[gt != '.'] %>% unique %>% length, .groups = 'drop')

# check numb of variants with X unique genotypes
ngt_per_loc %>% count(n_unq_gt, name = 'n_variants') %>% print

# some mono-allelic variants, but are these all 0 or all 1?
# mostly al zeros
monoal_vars = gt_data_long %>%
    semi_join(ngt_per_loc %>% filter(n_unq_gt == 1), by = 'uid') %>%
    filter(gt != '.') %>%
    distinct(uid, gt) %>%
    rename(unq_genotype = gt)
monoal_vars %>% count(unq_genotype, name = 'n_variants') %>% print

# filter these out
gt_data   = gt_data %>% anti_join(monoal_vars %>% filter(unq_genotype   == 0), by = c('uid'))
alig_data = alig_data %>% anti_join(monoal_vars %>% filter(unq_genotype == 0), by = c('uid'))

# |                   adjust POS and END positions for SVs                   |
# |   ---------------------------------------------------------------------  |
# |                        1. take REF/ALT alignments                        |
# | 2. get table of indels for each with pos relative to the reference start |
# |     3. add the REF start to each relative start to get adjusted start    |

# extract table of indels from alignment
locs_data = alig_data %>%
    select(uid, alig_tbls) %>%
    mutate(indels = map(alig_tbls, ~.x$indels)) %>%
    select(!alig_tbls) %>%
    unnest(indels) %>% 
    select(uid, rel_pos = start, rel_end = end, rel_width = width)

# one-to-many possible b/c we've split one REF/ALT alignment into multiple SVs in some cases
locs_data = alig_data %>% 
    # select(uid, chr, pos, REF, ALT, sv_type, ref_alt_diff = sv_size) %>%
    select(uid, sv_type, ref_alt_diff = sv_size) %>%
    left_join(locs_data, by = 'uid')

# adjust position
locs_data = locs_data %>% 
    left_join(dkeys %>% select(uid, chr, pos), by = 'uid') %>%
    mutate(pos_adj = pos + rel_pos,
	   end_adj = pos + rel_end) %>%
    select(uid, chr, pos, sv_type, ref_alt_diff, pos_adj, end_adj, rel_width)

# check how the sv size compares when calculating by taking difference between REF/ALT vs by alignment (rel_width)
locs_data %>% count(ref_alt_diff == rel_width, name = 'n_vars') %>% print

# determine strains with variant for each locus
strains_per_al = gt_data_long %>%
    filter(gt == 1) %>%
    group_by(uid) %>%
    summarise(strains = list(strain %>% unique), n_strains = strain %>% unique %>% length, .groups = 'drop')

# add the number of strains that each variant is present in and the identity of the strains as a list
# now we drop the original chr,pos,REF,ALT identifier and consolidate by chr, pos_adj, end_adj, sv_type
locs_data = locs_data %>% 
    left_join(strains_per_al, by = 'uid') %>%
    select(uid, chr, pos = pos_adj, end = end_adj, sv_type, sv_size = rel_width, strains, n_strains) %>%
    arrange(chr, pos, end)

# save object
sv_data = list(locs_data = locs_data, dkeys = dkeys)
out_file = path(out_dir, 'sv_data.rds')
saveRDS(sv_data, out_file)
