#!/home/momaksimov/anaconda3/envs/r4.2/bin/Rscript
#PBS -q hotel
#PBS -N my_script
#PBS -l nodes=1:ppn=1
#PBS -l walltime=01:00:00
#PBS -V
#PBS -M mikmaksi@gmail.com
#PBS -m abe
#PBS -t 1-19
#PBS -A gymreklab-group

# about: 
# 1) impute founder labels for STRs
# 2) join founder genotype to ri genotype by founder label
# 3) find denovo variants
# 4) calculate `delta_founder` metrics
# 5) save results
# one job per chromosome

# set current directory
setwd('/projects/ps-gymreklab/mikhail/072321_bxd_mutator_paper/workflows/8_colabcross')

# options
options(stringsAsFactors = FALSE)

# libraries
library(tidyverse)
library(fs)
library(qtl2)

# for debug
# chrom = '13'

# config
snp_dir = '../../data/colab_cross/snp_qtl2'
str_markers_file = '../../data/colab_cross/ref/str_segreg.tsv'
chrom = Sys.getenv('PBS_ARRAYID')
minprob = 0.99 # more conservative threshold
qtl_dir = '../../data/colab_cross/str_qtl2/'; dir_create(qtl_dir)
gts_dir = '../../data/colab_cross/str_gts/'; dir_create(gts_dir)
gt_vals_file = '../../data/colab_cross/str_gts/all_repcn_proc_nosegdup_nolowcr_segreg.rds'
founder_names  = 'A,B,C,D,E,F,G,H'
founder_labels = 'A,B,C,D,E,F,G,H'
out_dir = '../../data/colab_cross/denovo_info'; dir_create(out_dir)
print(sprintf('Processing %s, %s', chrom, minprob))

# read list of segregating str loci; NOTE: we don't need STR genotypes because we have no way to assign founders
str_markers = read_tsv(str_markers_file, col_types = cols(chr = 'c', pos = 'i', end = 'i')) %>%
    mutate(chr = str_replace(chr, 'chr', ''))
    
# select just a single chromosome
str_markers = str_markers %>% filter(chr == chrom)

# load snp qtl data
snp_data = list(cross = 'cross.rds', gmap = 'gmap.rds', pmap = 'pmap.rds', probs = 'probs.rds', phys_map = 'phys_map.rds') %>% map(~readRDS(path(snp_dir, .x)))

# calculate interpolation function for genomic coordinates
interp_fns = snp_data$phys_map %>%
    nest(data = !chr) %>%
    group_by(chr) %>%
    summarise(fn = map(data, ~approxfun(x = .x$pos, y = .x$cM, rule = 2)))

# make physical and genotype maps for strs
str_markers = str_markers %>%
    select(chr, pos, end) %>%
    nest(data = c(pos, end)) %>%
    left_join(interp_fns, by = 'chr') %>%
    mutate(data = map2(data, fn, ~.x %>% mutate(cM = .y(pos)))) %>%
    select(chr, data) %>%
    unnest(data) %>%
    mutate(Mb = pos/1e6) %>%
    mutate(loc_id = sprintf('chr%s_%s_%s', chr, pos, end))

# check interpolation
'Check interpolation of genetic coordinates' %>% print
str_markers %>% skimr::skim(cM, Mb) %>% print

# format str markers as qtl2 pseudomarkers
str_pseudomarkers = str_markers %>%
    mutate(chr = str_replace(chr, 'chr', '')) %>%
    split(.$chr) %>%
    map(~.x$cM %>% set_names(.x$loc_id))

# insert str markers into snp markers as pseudomarkers; tol = 0 to make sure all are added
padded_markers = insert_pseudomarkers(snp_data$gmap[chrom], pseudomarker_map = str_pseudomarkers, tol = 0)

# check dimensions or arrays
'Check number of markers' %>% print
list(str_pseudomarkers = str_pseudomarkers, snp_gmap = snp_data$gmap[chrom], padded_markers = padded_markers) %>% 
    map(~map_int(.x, length)) %>% print

# calculate genotype probabilities with pseudomarkers and save the result
probs = calc_genoprob(snp_data$cross[,chrom], map = padded_markers[chrom], quiet = FALSE)

# save probabilities
saveRDS(probs, path(qtl_dir, str_c(chrom, '_probs.rds')))

# pick the max marginal genotype based on probabilites
fou_lab = maxmarg(probs, return_char = TRUE, minprob = minprob, quiet = FALSE)

# form a tibble from founder labels
fou_lab = fou_lab[[chrom]] %>% 
    t() %>% as.data.frame %>% 
    rownames_to_column(var = 'marker') %>% 
    as_tibble

# subset fou_lab objects down to only strs
fou_lab = fou_lab %>% filter(marker %in% str_markers$loc_id)

# recode BB->B etc ...
fou_lab = fou_lab %>% 
    mutate(across(!marker, ~recode(.x, !!!( LETTERS[1:8] %>% set_names(str_c(LETTERS[1:8], LETTERS[1:8])) ) )))

# write outputs
saveRDS(fou_lab, path(qtl_dir, str_c(chrom, '_all_foulab_nosegdup_nolowcr_segreg.rds')))

# load numeric genotype values and select a single chromosome
gt_vals = readRDS(gt_vals_file) %>% filter(chr == str_c('chr', chrom))

# split out marker
fou_lab = fou_lab %>% separate('marker', c('chr', 'pos', 'end'), convert = TRUE)

# gather by strain
fou_lab = fou_lab %>% 
    pivot_longer(cols = !c(chr, pos, end), names_to = 'strain', values_to = 'gt')

# split comma separated genotype into integers
# `dplyr::separate` is too slow for this
gt_vals = gt_vals %>%
    pivot_longer(cols = !c(chr, pos, end), names_to = 'strain', values_to = 'gt')
split_gt = str_split(gt_vals$gt, pattern = ',', simplify = TRUE)
gt_vals = gt_vals %>% 
    mutate(GT_A = split_gt[,1], GT_B = split_gt[,2]) %>%
    mutate(across(c(GT_A, GT_B), as.integer)) %>%
    select(!gt)

# split founder names and labels
founder_strains = strsplit(founder_names, ',')[[1]]
founder_labs    = strsplit(founder_labels, ',')[[1]]

# split gt values for RI strains from gt values for founder strains
ri_gts      = gt_vals %>% filter(!strain %in% founder_strains)
founder_gts = gt_vals %>% filter(strain %in% founder_strains)

# join original and imputed genotypes to known genotypes of RI strains
ri_gts = ri_gts %>%
    left_join(fou_lab %>% rename(imp_gt = gt), by = c('chr', 'pos', 'end', 'strain'))

# join founder genotype for the expected founder
denovo_gts = ri_gts %>%
    rename(founder = imp_gt) %>%
    left_join(founder_gts, by = c('chr', 'pos', 'end', 'founder' = 'strain'), suffix = c('.ri', '.fou'))
    
# find denovos
denovo_gts = denovo_gts %>%
    filter(!is.na(GT_A.ri)) %>%
    filter(!is.na(GT_A.fou)) %>%
    filter(GT_A.ri  == GT_B.ri) %>% # hom ri
    filter(GT_A.fou == GT_B.fou) %>% # hom fou
    filter(GT_A.ri  != GT_B.fou)

# calculate delta founder
denovo_gts = denovo_gts %>%
    mutate(fou_rn = GT_A.fou) %>%
    mutate(delta_fou = GT_A.ri - fou_rn) %>%
    mutate(expand_sign = sign(delta_fou)) %>%
    mutate(expand_type = case_when(expand_sign == -1 ~'contr', 
				   expand_sign == 1 ~ 'expan',
				   expand_sign == 0 ~ 'equal')) %>%
    mutate(delta_fou = abs(delta_fou))

# final formatting
denovo_gts = denovo_gts %>% 
    select(chr, pos, end, strain, RN_A = GT_A.ri, RN_B = GT_B.ri, founder, fou_rn, delta_fou, expand_sign, expand_type)

# check expand vs contract
denovo_gts %>% count(expand_type) %>% mutate(perc = n*100/sum(n)) %>% print

# write outputs
saveRDS(denovo_gts, path(out_dir, str_c(chrom, '_denovo_gts.rds')))
