#!/home/momaksimov/anaconda3/envs/r4.2/bin/Rscript

# about: process raw str repeat counts
# remove loci that overlap with segmental duplications
# remove loci with low call rates
# drop BXD210 from analysis since there is no snp information for it

# options
options(stringsAsFactors = FALSE)

# libraries
library(tidyverse)
library(fs)
library(DBI)
library(dbplyr)
library(progress)

# config
geno_file = '../../data/str_gts/all_repcn_raw.tsv.gz'
out_dir = '../../data/str_gts'; dir_create(out_dir)

# load strain info
strain_info = readRDS('../../info/strain_info.rds')

# connect to ucsc sql and download segdup regions
conn = dbConnect(RMySQL::MySQL(), 
		 host = 'genome-mysql.soe.ucsc.edu', 
		 user = 'genome', 
		 port = 3306,
		 dbname = 'mm10') 
segdup_reg = tbl(conn, sql('select chrom, chromStart, chromEnd from genomicSuperDups')) %>% 
    collect() %>% set_names(c('chr', 'pos', 'end'))
dbDisconnect(conn)

# read the str regions file
geno = read_tsv(geno_file, col_types = cols(.default = 'c'))

# reformat header
header = names(geno)
header = header %>% 
    str_replace('#', '') %>% 
    str_replace(':REPCN', '') %>% 
    str_replace('\\[\\d+\\]', '') %>% 
    str_trim()
geno = geno %>% set_names(header)

# recode strains
recodes = strain_info %>% filter(is_seq_str) %>% pull(short_name) 
names(recodes) = strain_info %>% filter(is_seq_str) %>% pull(bxd_id) 
geno = geno %>% select(chr = CHROM, pos = POS, end = END, all_of(recodes))

# drop BXD210
geno = geno %>% select(!BXD210)

# convert types
geno = geno %>% mutate(across(c(pos, end), as.integer))

# convert to '.' to NA
geno = geno %>% mutate(across(c(matches('BXD'), C57BL, DBA), ~if_else(.x == '.', NA_character_, .x)))

# save processed version
saveRDS(geno, path(out_dir, 'all_repcn_proc.rds'))

# combine segdup and str regions into a list
bed_regions = list(str = geno %>% select(chr, pos, end), segdup = segdup_reg %>% select(chr, pos, end))

# keep only releveant chromsomes
bed_regions = bed_regions %>% 
    map(~.x %>% filter(chr %in% str_c('chr', 1:19)))

# sort loci in lists
bed_regions = bed_regions %>% 
    map(~.x %>% mutate(chr = fct_relevel(chr, str_c('chr', 1:19))) %>% arrange(chr, pos, end))

# write temporary .bed files and run `betools intersect` to find loci that overlap with segdup regions
tmp_files = list(a = tempfile(), b = tempfile())
walk2(bed_regions, tmp_files, ~write_tsv(.x, .y, col_names = FALSE))
segdup_strs = read_tsv(
    pipe(sprintf("bedtools intersect -a %s -b %s -u -sorted", tmp_files$a, tmp_files$b)),
    col_names = c('chr', 'pos', 'end'), 
    col_types = cols('c', 'i', 'i')
)
map(tmp_files, file_delete)

# save list of segdup strs
segdup_strs %>% 
    mutate(chr = fct_relevel(chr, str_c('chr', 1:19))) %>%
    arrange(chr, pos, end) %>%
    write_tsv('../../data/ref/str_segdup.tsv')

# calculate per-locus call rate
# CHECK: 154 strains in original analysis (BXD210 already removed) - same as here
# CHECK: filter if > 50, not >= 50
miss_strains_per_loc = geno %>%
    mutate(miss_strains = rowSums(is.na(across(c(matches('BXD'), C57BL, DBA))))) %>%
    select(chr, pos, end, miss_strains)
lowcr_loci = miss_strains_per_loc %>% filter(miss_strains > 50)

# save list of segdup strs
lowcr_loci %>% 
    mutate(chr = fct_relevel(chr, str_c('chr', 1:19))) %>%
    arrange(chr, pos, end) %>%
    select(chr, pos, end) %>%
    write_tsv('../../data/ref/str_lowcr.tsv')

# save a version of genotypes with segdup and lowcr loci filtered
geno_nosegdup_nolowcr = geno %>% 
    anti_join(segdup_strs, by = c('chr', 'pos', 'end')) %>%
    anti_join(lowcr_loci, by = c('chr', 'pos', 'end'))
saveRDS(geno_nosegdup_nolowcr, path(out_dir, 'all_repcn_proc_nosegdup_nolowcr.rds'))

{ # DEBUGGING
    # check with previous lists
    data_dir = '/projects/ps-gymreklab/mikhail/090520_unified_workflow/data/ref/exclude'
    ori_regions = map_df(
	list(lowcr = 'lowcr_strs.bed',
	     segdup = 'segdup_strs.bed',
	     zv_hetonly = 'zv_hetonly_strs.bed'),
	~read_tsv(path(data_dir, .x), 
		  col_names = c('chr', 'pos', 'end'), 
		  col_types = cols(chr = 'c', pos = 'i', end = 'i')), 
	.id = 'type')
    ori_regions %>%
	bind_rows(extra %>% mutate(type = 'extra')) %>%
	group_by(chr, pos, end) %>%
	summarise(lists = str_c(type %>% unique, collapse = ', '), .groups = 'drop') %>% 
	count(lists)
    # lists                          n
    # <chr>                      <int>
    # extra                         93 only these are of concern
    # lowcr                     109915
    # lowcr, segdup               2047
    # lowcr, segdup, zv_hetonly  11734
    # lowcr, zv_hetonly         280399
    # segdup                      1334
    # segdup, extra                 35
    # segdup, zv_hetonly         17369
    # segdup, zv_hetonly, extra   1838
    # zv_hetonly                671634
    # zv_hetonly, extra           2984

    # ones of concerns are those that aren't filtered in other lists
    real_extra = extra %>% anti_join(ori_regions)
    lowcr_loci %>% semi_join(real_extra)
    geno %>% 
	semi_join(real_extra) %>%
	pivot_longer(c(matches('BXD'), C57BL, DBA), 'strain', 'gt') %>%
	group_by(chr, pos, end) %>%
	summarise(sum(is.na(value)))
    # NOTE: this version is correct; found a bug in previous code that didn't distinguish b/w "." and "./."
}
