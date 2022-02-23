#!/home/momaksimov/anaconda3/envs/r4.2/bin/Rscript

# about: remove any loci that overlap with segdup regions in mm10 

# options
options(stringsAsFactors = FALSE)

# libraries
library(tidyverse)
library(fs)
library(DBI)
library(dbplyr)

# config
snp_regions_file = '../../data/snp_gts/nohetfou_nomiss_noident.tsv.gz'
out_dir = '../../data/snp_gts'; dir_create(out_dir)

# connect to ucsc sql and download segdup regions
conn = dbConnect(RMySQL::MySQL(), 
		 host = 'genome-mysql.soe.ucsc.edu', 
		 user = 'genome', 
		 port = 3306,
		 dbname = 'mm10') 
segdup_reg = tbl(conn, sql('select chrom, chromStart, chromEnd from genomicSuperDups')) %>% 
    collect() %>% set_names(c('chr', 'pos', 'end'))
dbDisconnect(conn)

# read founder genotypes
snp_regions = read_tsv(snp_regions_file, 
		       col_names = c('chr', 'pos', 'end', 'C57', 'DBA'), 
		       col_types = cols('c', 'i', 'i', 'c', 'c'), 
		       comment = '#')

# check distribution of founder genotypes to confirm previous script is working
# snp_regions %>% count(C57 == DBA)
# snp_regions %>% count(C57, DBA)

# combine segdup and snp regions into a list
bed_regions = list(snp = snp_regions, segdup = segdup_reg)

# keep only releveant chromsomes
bed_regions = bed_regions %>% 
    map(~.x %>% filter(chr %in% str_c('chr', 1:19)))

# sort loci in each list
bed_regions = bed_regions %>% 
    map(~.x %>% mutate(chr = fct_relevel(chr, str_c('chr', 1:19))) %>% arrange(chr, pos, end))

# run bed files to temporary location for `bedtools intersect`
# get list of loci not intersecting with segdup regions
tmp_files = list(a = tempfile(), b = tempfile())
walk2(bed_regions, tmp_files, ~write_tsv(.x, .y, col_names = FALSE))
nonsegdup_snps = read_tsv(
    pipe(sprintf("bedtools intersect -a %s -b %s -v -sorted", tmp_files$a, tmp_files$b)),
    col_names = c('chr', 'pos', 'end', 'C57', 'DBA'), 
    col_types = cols('c', 'i', 'i', 'c', 'c')
)
map(tmp_files, file_delete)

# sort by locus position
nonsegdup_snps = nonsegdup_snps %>% 
    mutate(chr = fct_relevel(chr, str_c('chr', 1:19))) %>% 
    arrange(chr, pos, end)

# check counts
nonsegdup_snps %>% count(chr)

# write output
# no header b/c needs to be .bed like for querying
write_tsv(nonsegdup_snps, path(out_dir, 'nohetfou_nomiss_noident_nosegdup.tsv.gz'), col_names = FALSE)
