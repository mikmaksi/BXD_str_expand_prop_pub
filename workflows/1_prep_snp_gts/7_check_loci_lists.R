#!/home/momaksimov/anaconda3/envs/r4.2/bin/Rscript

# about: check that no snps of interest have been missed

# options
options(stringsAsFactors = FALSE, dplyr.summarise.inform = FALSE)

# libraries
library(tidyverse)
library(cowplot)
library(fs)

# read list of snp loci of interest
snp_regions_file = '../../data/snp_gts/nohetfou_nomiss_noident_nosegdup.tsv.gz'
snp_regions = read_tsv(snp_regions_file, 
		       col_names = c('chr', 'pos', 'end', 'C57', 'DBA'), 
		       col_types = cols('c', 'i', 'i', 'c', 'c'), 
		       comment = '#')

# get list of loci for which raw genotype values have been pulled
cmd = "zcat ../../data/snp_gts/raw/* | awk '$0!~\"CHROM\" {print $1,$2}' OFS='\t'"
rawgts_loci = read_tsv(pipe(cmd), col_names = c('chr', 'pos'), col_types = cols(chr = 'c', pos = 'i'))
rawgts_loci %>% count(chr)

# summarise requested (req) vs. returned (ret) loci sets
comp = list(
    req = snp_regions %>% select(chr, pos),
    ret = pulled_loci %>% select(chr, pos)
)
setdiff(comp$req, comp$ret)
setdiff(comp$ret, comp$req)

# get list of loci from for-PLINK vcfs
vcf_loci = map_df(str_c('chr', 1:19), function(.x) {
    print(.x)
    cmd = sprintf("bcftools query -f '%%CHROM\t%%POS\n' ../../data/snp_gts/vcf/%s.vcf.gz", .x)
    read_tsv(pipe(cmd), col_names = c('chr', 'pos'), col_types = cols(chr = 'c', pos = 'i'))
})

# summarise requested (req) vs. returned (ret) loci sets
comp = list(
    req = snp_regions %>% select(chr, pos),
    ret = vcf_loci %>% select(chr, pos)
)
setdiff(comp$req, comp$ret)
setdiff(comp$ret, comp$req)
