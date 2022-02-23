#!/home/momaksimov/anaconda3/envs/r4/bin/Rscript

# about: calculate various qtl2 objects from snp genotypes for CC strains

# clean variables
rm(list = ls())

# libraries
library(qtl2)
library(tidyverse)
library(fs)

# update marker map and geno object from cross
update_objs = function(cross_obj) {
    # convert qtl2 objects to tibbles
    pmaps    = map_df(cross_obj$pmap, ~tibble(marker = names(.x), pos = as.integer(as.numeric(.x)*1e6)), .id = 'chr')
    ri_genos = map_df(cross_obj$geno, ~.x %>% t() %>% as.data.frame %>% rownames_to_column(var = 'marker') %>% as_tibble)
    list(pmaps = pmaps, ri_genos = ri_genos)
}

# config
in_dir = '../../data/colab_cross/unc_data'; dir_create(in_dir)
out_dir = '../../data/colab_cross/snp_qtl2'; dir_create(out_dir)

# make list to store similar data at different points of analysis
cross_data = list()
pmaps      = list()
ri_genos   = list()
probs      = list()
phys_maps  = list()

# load cross data
cross_data$full = read_cross2(str_c(in_dir, '/', 'cc.zip'), quiet=FALSE)

# convert qtl2 objects to tibbles
ret = update_objs(cross_data$full)
pmaps$full    = ret$pmaps
ri_genos$full = ret$ri_genos
probs$full    = NA # no point to generate

# check unique values in genotype table
# ri_genos$full %>% select(!marker) %>% unlist %>% table(useNA = 'always')
#      0       1       3    <NA> 
# 205256 3696928 3691542       0 

# pivot ri genotypes
ri_genos_t = ri_genos$full %>% 
    pivot_longer(matches('^CC'), names_to = 'strain', values_to = 'gt') %>%
    pivot_wider(id_cols = strain, names_from = marker, values_from = gt)

# calculate number of genotypes per locus
n_gt_per_loc = apply(ri_genos_t %>% select(-strain), 2, 
		     function(.x) .x %>% discard(~.x == 0) %>% unique %>% length)
# n_gt_per_loc %>% table(useNA = 'always')
# there are some monoallelic loci (not counting missing value)
#     1      2  <NA>
#   688 109366     0

# spot check monoallelic loci
# ri_genos$full %>% 
#     filter(marker %in% (n_gt_per_loc %>% keep(~.x == 1) %>% names)) %>% 
#     pivot_longer(matches('^CC'), names_to = 'strain', values_to = 'gt') %>%
#     filter(marker == sample(unique(marker), 1)) %>%
#     count(gt)

# remove monoallelic loci and generate a new cross objecte
cross_data$nozv = drop_markers(cross_data$full, (n_gt_per_loc %>% keep(~.x == 1) %>% names))

# update qtl2 objects
ret = update_objs(cross_data$nozv)
pmaps$nozv    = ret$pmaps
ri_genos$nozv = ret$ri_genos
probs$nozv    = NA # no point to generate

# make list of markers that overlap exacly
dup_marks = pmaps$nozv %>%
    semi_join(pmaps$nozv %>% 
	      count(chr, pos, name = 'n_mark') %>%
	      filter(n_mark > 1), by = c('chr', 'pos')) %>%
    arrange(chr, pos)

# join number of missing genotype per marker to duplicated markers
dup_marks = dup_marks %>%
    left_join(n_missing(cross_data$nozv, by = 'marker') %>% tibble(marker = names(.), n_miss = .), by = 'marker')

# for each set of duplicated markers drop the one with most missing values
drop_marks = dup_marks %>% arrange(n_miss) %>% distinct(chr, pos, .keep_all = TRUE)
cross_data$nozv_nodup = drop_markers(cross_data$nozv, (drop_marks %>% pull(marker)))

# update qtl2 objects
ret = update_objs(cross_data$nozv_nodup)
pmaps$nozv_nodup    = ret$pmaps
ri_genos$nozv_nodup = ret$ri_genos

# calculate genotype probabilities
probs$nozv_nodup = calc_genoprob(cross_data$nozv_nodup, quiet = FALSE, cores = 1)

# construct readable physical/genetic map files
phys_maps = map(names(cross_data) %>% set_names(.), function(.x) {
    map_df(cross_data[[.x]]$pmap, ~tibble(marker = names(.x), pos = as.integer(as.numeric(.x)*1e6), Mb = as.numeric(.x)), .id = 'chr') %>%
	left_join(map_df(cross_data[[.x]]$gmap, ~tibble(marker = names(.x), cM = as.numeric(.x)), .id = 'chr'), by = c('chr', 'marker'))
})

# calculate kinship
kinships = map(probs, function(.x) if(length(.x) == 1 && is.na(.x)) { NA } else { calc_kinship(.x, type = 'loco') })

# check
list(cross = cross_data, pmap = pmaps, ri_genos = ri_genos, 
     probs = probs, phys_map = phys_maps, kinship = kinships) %>% map(~.x %>% names)

# collate by analysis type
# .x = 'nozv_nodup'
master_list = map(c('nozv_nodup') %>% set_names(.), function(.x) {
    list(
	cross    = cross_data[[.x]],
	gmap     = cross_data[[.x]]$gmap,
	pmap     = cross_data[[.x]]$pmap,
	probs    = probs[[.x]],
	kinship  = kinships[[.x]],
	phys_map = phys_maps[[.x]],
	out_dir  = out_dir
    )
})

# save data objects
walk(master_list, function(.x) {
    saveRDS(.x$cross,    path(.x$out_dir, "cross.rds"))
    saveRDS(.x$gmap,     path(.x$out_dir, "gmap.rds"))
    saveRDS(.x$pmap,     path(.x$out_dir, "pmap.rds"))
    saveRDS(.x$probs,    path(.x$out_dir, "probs.rds"))
    saveRDS(.x$kinship,  path(.x$out_dir, "kinship.rds"))
    saveRDS(.x$phys_map, path(.x$out_dir, 'phys_map.rds'))
    mutate(fname = path %>% path_file) %>%
    print(n = nrow(.))
