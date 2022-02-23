#!/home/momaksimov/anaconda3/envs/r4.2/bin/Rscript
#PBS -q hotel
#PBS -N my_script
#PBS -l nodes=1:ppn=1
#PBS -l walltime=00:30:00
#PBS -o job.o
#PBS -e job.e 
#PBS -V
#PBS -M mikmaksi@gmail.com
#PBS -m abe
#PBS -t 1-19
#PBS -A gymreklab-group

# about: make a locus-by-strain file with founder haplotypes for snps
# match each RI genotype to C57 and to DBA and assign B or D respectively
# if not matching then assign NA
# fast column-wise implementation

# cd into main analysis directory
setwd('/projects/ps-gymreklab/mikhail/072321_bxd_mutator_paper/workflows/1_prep_snp_gts')

# get job environmental variable
PBS_ARRAYID = as.integer(Sys.getenv('PBS_ARRAYID'))

# debug
# PBS_ARRAYID=1

# options
options(stringsAsFactors = FALSE)

# libraries
library(tidyverse)
library(fs)
library(optparse)

# config
data_dir = '../../data/snp_gts/raw'
out_dir = '../../data/snp_gts/fou_lab'; dir_create(out_dir)
chrom = str_c('chr', PBS_ARRAYID)

# load strain info
strain_info = readRDS('../../info/strain_info.rds')

# user output
print(sprintf('Processing %s', chrom))

# read genotypes
geno = read_tsv(path(data_dir, str_c(chrom, '.tsv.gz')))

# reformat header
header = names(geno)
header = header %>% 
    str_replace('#', '') %>% 
    str_replace(':GT', '') %>% 
    str_replace('\\[\\d+\\]', '') %>% 
    str_trim()
header[header == c('POS')] = c('POS', 'END')
geno = geno %>% set_names(header)

# recode strains
recodes = strain_info %>% filter(is_seq_snp) %>% pull(short_name) 
names(recodes) = strain_info %>% filter(is_seq_snp) %>% pull(bxd_id) 
geno = geno %>% select(chr = CHROM, pos = POS, end = END, all_of(recodes))

# convert to B/D labels
geno_recode = geno %>% 
    mutate(across(matches('BXD'), ~case_when(.x == C57BL ~ 'B', .x == DBA ~ 'D', TRUE ~ NA_character_)))

# check counts
geno_recode %>% 
    pivot_longer(matches('BXD'), names_to = 'strain', values_to = 'gt') %>%
    count(strain, gt) %>%
    pivot_wider(id_cols = strain, names_from = gt, values_from = n)

# check
# geno_recode %>% names

# remove founders
geno_recode = geno_recode %>% select(!c(C57BL, DBA))

# combine into locus
geno_recode = geno_recode %>%
    mutate(across(c(pos, end), as.integer)) %>%
    arrange(pos, end) %>%
    unite('locus', c('chr', 'pos', 'end'))

# write output
write_tsv(geno_recode, path(out_dir, str_c(chrom, '.tsv.gz')))
