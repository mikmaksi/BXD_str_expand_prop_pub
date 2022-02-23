#!/home/momaksimov/anaconda3/envs/r4/bin/Rscript

# about: Rdata binary file from dryad archive for Attie dataset contains multiple pieces of information
# extract the useful ones

# options
options(stringsAsFactors = FALSE, dplyr.summarise.inform = FALSE)

# libraries
library(tidyverse)
library(cowplot)
library(fs)
library(qtl2)

# load Attie dataset
load('../../data/diversity_outbred/Attie_DO378_eQTL_viewer_v1.Rdata')

# check objects
# ls()
# [1] "dataset.exvivo"         "dataset.invivo"         "dataset.islet.proteins" "dataset.islet.rnaseq"  
# [5] "dataset.liver.lipids"   "dataset.plasma.lipids"  "ensembl.version"        "genoprobs"             
# [9] "K"                      "map"                    "markers"               

# ensebml version
ensembl.version
# 75 used GRCm38 which is mm10

# markers (some of these are pseudo markers
# marker.id chr     pos    cM      bp
# <chr>     <chr> <dbl> <dbl>   <dbl>
# 1_3000000 1      3     1.45 3000000
# 1_3041392 1      3.04  1.47 3041392
# 1_3346528 1      3.35  1.48 3346528

# genotype probabilities (already in r/qtl2 format
# genoprobs[[13]] %>% dimnames %>% map(head)
# [[1]]
# [1] "DO021" "DO022" "DO023" "DO024" "DO025" "DO026"
# 
# [[2]]
# [1] "A" "B" "C" "D" "E" "F"
# 
# [[3]]
# [1] "13_3000000" "13_3396518" "13_3793036" "13_3803732" "13_3807812" "13_3811891"

# convert genotype probabilities to founder labels using qtl2
fou_lab = maxmarg(genoprobs, minprob = 0)
object.size(fou_lab) %>% print(unit = 'Mb')

# convert to a tibble
fou_lab_tib = fou_lab %>% 
    map_df(~.x %>% t() %>% as.data.frame %>% rownames_to_column(var = 'marker.id') %>% as_tibble, .id = 'chr')

# check on islet RNA data
# dataset.islet.rnaseq %>% names
data_ele = c("annot.mrna", "annot.samples", "covar.matrix", "covar.info", "data", "datatype", "display.name", "lod.peaks")
rnaseq = dataset.islet.rnaseq
rnaseq[[data_ele[1]]]
rnaseq[[data_ele[2]]]
rnaseq[[data_ele[3]]] %>% as.data.frame %>% head
rnaseq[[data_ele[4]]]

# raw and normalized? data
rnaseq[[data_ele[5]]] %>% map(dim)
# $raw
# [1]   378 21771
# 
# $rz
# [1]   378 21771

# try qtl mapping with Msh3
embl_id = rnaseq$annot.mrna %>% filter(symbol == 'Msh3') %>% pull(gene.id)
qtl_res = scan1(genoprobs[,13],
      pheno = rnaseq$data$rz[,embl_id, drop = FALSE],
      # pheno = rnaseq$data$raw[,embl_id, drop = FALSE],
      kinship = K['13'],
      addcovar = rnaseq$covar.matrix)
qtl_res = qtl_res %>% as.data.frame %>% rownames_to_column(var = 'marker') %>% as_tibble  %>% rename(pheno = embl_id)

# check the result against precomputed values values
# here result: matches when "rz" is used, but not when "raw" is used
qtl_res %>% arrange(desc(pheno)) %>% slice(1)

# pre-comp values
( rnaseq$lod %>% map(~.x %>% filter(gene.id == embl_id)) )$additive
# # A tibble: 1 x 2
#   marker      pheno
#   <chr>       <dbl>
# 1 13_92367481  39.6
# # A tibble: 4 x 3
#   gene.id            marker.id     lod
#   <chr>              <chr>       <dbl>
# 1 ENSMUSG00000014850 12_31064532  7.57
# 2 ENSMUSG00000014850 13_92367481 39.6 
# 3 ENSMUSG00000014850 16_59515531  7.18
# 4 ENSMUSG00000014850 17_79487020  6.82

# combine into a single named list
# only save chr13
chrom = 13
DO_gt = list(
    genoprobs = genoprobs[,chrom],
    pmap = map[chrom],
    phys_map = markers %>% filter(chr == chrom),
    kinship = K[chrom],
    fou_lab = fou_lab_tib %>% filter(chr == chrom))
DO_gene_expr = list(
    annot.gene = rnaseq$annot.mrna,
    annot.smp = rnaseq$annot.samples,
    pheno = rnaseq$data$rz,
    covar = rnaseq$covar.matrix)
# list(DO_gt = DO_gt, DO_gene_expr = DO_gene_expr) %>% walk(~object.size(.x) %>% print(unit = 'Mb'))

# save data
out_dir = '../../data/diversity_outbred'; dir_create(out_dir)
saveRDS(DO_gt, path(out_dir, str_c('DO_gt_', chrom, '.rds')))
saveRDS(DO_gene_expr, path(out_dir, 'DO_gene_expr.rds'))
