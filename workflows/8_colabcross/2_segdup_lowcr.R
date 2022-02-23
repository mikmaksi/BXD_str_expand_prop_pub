#!/home/momaksimov/anaconda3/envs/r4.2/bin/Rscript

# about: process raw str repeat counts
# remove loci that overlap with segmental duplications
# remove loci with low call rates

# options
options(stringsAsFactors = FALSE)

# libraries
library(tidyverse)
library(fs)
library(DBI)
library(dbplyr)
library(progress)

# config
geno_file = '../../data/colab_cross/str_gts/all_repcn_raw.tsv.gz'
out_dir   = '../../data/colab_cross/str_gts'; dir_create(out_dir)
info_dir  = '../../info'; dir_create(info_dir)
ref_dir   = '../../data/colab_cross/ref'; dir_create(ref_dir)

# load strain info
strain_info = read_tsv(path(info_dir, 'cc_strain_info.tsv'), comment = '#',
    col_types = cols(
      strain         = 'c',
      long_name      = 'c',
      m_chrom        = 'c',
      y_chrom        = 'c',
      funnel_code    = 'c',
      gen_inbreeding = 'i',
      n_founders     = 'i',
      origin         = 'c',
      origin_id      = 'i',
      long_name      = 'c')
)

# filter strain info (these strains don't have long names)
strain_info = strain_info %>% 
    filter(!is.na(long_name)) %>%
    filter(strain != 'CC071') # CC071 is not complete; leave out of analysis for now

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

# load strain list 
strain_list = read_tsv(path(info_dir, 'cc_strain_list'), 
		       col_names = 'strain', col_types = cols(strain = 'c'), comment = '#')$strain

# only keep these
geno = geno %>% select(chr = CHROM, pos = POS, end = END, all_of(strain_list))

# recode strains
recodes = strain_info %>% pull(long_name) 
names(recodes) = strain_info %>% pull(strain) 
geno = geno %>% select(chr, pos, end, all_of(recodes))

# convert types
geno = geno %>% mutate(across(c(pos, end), as.integer))

# convert to '.' to NA
geno = geno %>% mutate(across(c(matches('^CC'), LETTERS[1:8]), ~if_else(.x == '.', NA_character_, .x)))

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
    write_tsv(path(ref_dir, 'str_segdup.tsv'))

# calculate per-locus call rate
# CHECK: 154 strains in original analysis (BXD210 already removed) - same as here
# CHECK: filter if > 50, not >= 50
miss_strains_per_loc = geno %>%
    mutate(miss_strains = rowSums(is.na(across(c(matches('^CC'), LETTERS[1:8]))))) %>%
    select(chr, pos, end, miss_strains)
lowcr_loci = miss_strains_per_loc %>% filter(miss_strains > 25)

# save list of segdup strs
lowcr_loci %>% 
    mutate(chr = fct_relevel(chr, str_c('chr', 1:19))) %>%
    arrange(chr, pos, end) %>%
    select(chr, pos, end) %>%
    write_tsv(path(ref_dir, 'str_lowcr.tsv'))

# save a version of genotypes with segdup and lowcr loci filtered
geno_nosegdup_nolowcr = geno %>% 
    anti_join(segdup_strs, by = c('chr', 'pos', 'end')) %>%
    anti_join(lowcr_loci, by = c('chr', 'pos', 'end'))
saveRDS(geno_nosegdup_nolowcr, path(out_dir, 'all_repcn_proc_nosegdup_nolowcr.rds'))
