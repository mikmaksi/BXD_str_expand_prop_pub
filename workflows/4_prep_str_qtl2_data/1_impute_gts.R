#!/home/momaksimov/anaconda3/envs/r4.2/bin/Rscript

# about: impute missing founder labels for strs

# options
options(stringsAsFactors = FALSE)

# libraries
library(tidyverse)
library(fs)
library(qtl2)

## functions

# for imputing genotypes from cross
impute_gt_func = function(cross_obj, ncores = 1) {
    # calculate genotype probabilities
    gprobs = calc_genoprob(cross_obj, cores = ncores, quiet = FALSE)

    # pick the max marginal genotype based on probabilites
    impute_gt = maxmarg(gprobs, cores = ncores, quiet = FALSE)

    # replace na values with zeros 
    for (chrom in 1:length(impute_gt)) {
	mp = impute_gt[[chrom]]
	mp[is.na(mp)] = 0
	impute_gt[[chrom]] = mp
    }

    # output
    return(impute_gt)
}

# function for converting from genotype matrix to tidy format
geno_to_tbl = function(geno_list) {
    geno_tbl = lapply(geno_list, function(x) {
		ret = as.data.frame(x) %>% rownames_to_column(var = 'strain')
		ret = ret %>% gather(locus, gt, -strain)
	     }) %>% bind_rows %>% as_tibble
    return(geno_tbl)
}

# fix imputation result
fix_imp_tbl = function(imp_tbl) {
    imp_tbl = imp_tbl %>%
	mutate_at(vars(matches('gt\\.')), as.integer) %>%
	mutate(gt.imp.fixed = gt.imp) %>%
	mutate(gt.imp.fixed = if_else(gt.ori != 0, gt.ori, gt.imp.fixed))
    
    # check
    imp_tbl %>% 
	filter(gt.ori != 0) %>%
	filter(gt.ori != gt.imp.fixed) %>% nrow %>% print

    # output
    return(imp_tbl)
}

# function for summarizing imputation results
imputation_smry = function(ori_geno_list, imp_geno_list) {
    # convert from geno to tbl
    ori_geno_tbl = geno_to_tbl(ori_geno_list)
    imp_geno_tbl = geno_to_tbl(imp_geno_list)

    # check values
    # ori_geno_tbl %>% distinct(gt)
    # imp_geno_tbl %>% distinct(gt)
    
    # check number of missing genotypes 
    # ori_geno_tbl %>% summarise(frac_miss = mean(gt == 0))
    # imp_geno_tbl %>% summarise(frac_miss = mean(gt == 0))

    # join the original and imputed genotypes together
    comb = ori_geno_tbl %>% left_join(imp_geno_tbl, by = c('strain', 'locus'), suffix = c('.ori', '.imp'))
    
    # summarize by category
    imp_smry = comb %>% mutate(label = case_when(
	    (gt.ori == 0) & (gt.imp != 0) ~ 'impute_success',# successfully imputed genotypes
	    (gt.ori == 0) & (gt.imp == 0) ~ 'impute_fail', # not successfully imputed genotypes (still missing)
	    (gt.ori != 0) & (gt.imp == 0) ~ 'impute_wrong', # poorly imputed (were note missing, but now are)
	    (gt.ori != 0) & (gt.imp != 0) ~ 'impute_not_needed', # 
	    TRUE ~ 'other'
	)) %>% count(label) %>% mutate(perc = n*100/sum(n))

    # factorize output
    imp_smry = imp_smry %>% 
	mutate(label = factor(label, 
	       levels = c('impute_not_needed', 'impute_success', 'impute_fail', 'impute_wrong'))) %>% 
	arrange(as.numeric(label))

    # check if there were any changes to genotypes where genotypes were already known
    known_changes_smry = comb %>% 
	filter((gt.ori != 0) & (gt.imp != 0)) %>%
	mutate(status = if_else(gt.ori == gt.imp, 'not_changed', 'changed')) %>%
	count(status) %>% mutate(perc = n*100/sum(n)) 

    # print results
    cat('Summary of genotypes by type;
   "impute_not_needed": gt known before and after impute
   "impute_success": gt unknown before, known after impute
   "impute_fail": gt unknown before and after impute
   "impute_wrong": gt known before, unkown after\n')
   print(imp_smry)

   #
   cat('\n')

   cat('Summary of changes to known loci before imputation ("impute_not_needed" group);\n')
   print(known_changes_smry)

   # outputs
   return(list(imp_smry = imp_smry, known_changes_smry = known_changes_smry, ori_imp_joined = comb))
}

# config
geno_file = '../../data/str_gts/all_foulab_nosegdup_nolowcr_noalshared_padded.rds'
out_dir = '../../data/str_qtl2'; dir_create(out_dir)

# read geno
geno = readRDS(geno_file)

# check the number of NAs per strain
nas_per_strain = apply(geno %>% select(matches('BXD')), 2, function(x) 100*mean(is.na(x))) %>% 
    as.data.frame %>%
    rownames_to_column(var = 'strain') %>% as_tibble %>%
    set_names(c('strain', 'perc_na')) %>%
    arrange(desc(perc_na))
nas_per_strain %>% skimr::skim(perc_na)

# check unique values
# geno %>% select(matches('BXD')) %>% unlist %>% unique

# load snp coordinates for interpolation
gw_dir = '../../data/snp_qtl2/gw'
gw_pmap = readRDS(path(gw_dir, 'pmap.rds')) %>% map_df(~tibble(marker = names(.x), pos = .x), .id = 'chr')
gw_gmap = readRDS(path(gw_dir, 'gmap.rds')) %>% map_df(~tibble(marker = names(.x), cM = .x), .id = 'chr')
interp_fns = gw_gmap %>%
    left_join(gw_pmap, by = c('chr', 'marker')) %>%
    mutate(chr = str_c('chr', chr)) %>%
    nest(data = !chr) %>%
    group_by(chr) %>%
    summarise(fn = map(data, ~approxfun(x = .x$pos, y = .x$cM, rule = 2)))

# make physical and genotype maps
phys_map = geno %>% select(chr, pos, end)
geno_map = geno %>%
    select(chr, pos, end) %>%
    nest(data = c(pos, end)) %>%
    left_join(interp_fns, by = 'chr') %>%
    mutate(data = map2(data, fn, ~.x %>% mutate(cM = .y(pos/1e6)))) %>%
    select(chr, data) %>%
    unnest(data)
# geno_map %>% group_by(chr) %>% skimr::skim(cM)
	
# write outputs in qtl2 format
phys_map %>% 
    unite('marker', c('chr', 'pos', 'end'), remove = FALSE) %>%
    select(marker, chr, pos) %>%
    mutate(chr = str_replace(chr, 'chr', ''),
	   pos = pos/1e6) %>% 
    write_csv(path(out_dir, 'pmap.csv'))
geno_map %>% 
    unite('marker', c('chr', 'pos', 'end'), remove = FALSE) %>%
    select(marker, chr, pos = cM) %>%
    mutate(chr = str_replace(chr, 'chr', '')) %>% 
    write_csv(path(out_dir, 'gmap.csv'))
geno %>%
    unite('marker', c('chr', 'pos', 'end')) %>%
    write_csv(path(out_dir, 'geno.csv'))

# load the cross object
config_file = path(out_dir, 'config.json')
write_control_file(config_file, 
    description = 'BXD strs gw',
    crosstype = 'risib',
    sep = ",",
    na.strings = c("-", "NA"),
    comment.char = "#",
    geno_file = '../../data/str_qtl2/geno.csv',
    gmap_file = '../../data/str_qtl2/gmap.csv',
    pmap_file = '../../data/str_qtl2/pmap.csv',
    alleles = c("B", "D"),
    geno_codes = c(B = 1, D = 2),
    geno_transposed = TRUE, 
    overwrite = TRUE)
bxd_cross = read_cross2(config_file, quiet = FALSE)
# file_delete(config_file)

# impute genotypes
impute_gt = impute_gt_func(bxd_cross, 1)

# check matrices
impute_gt$`1` %>% .[1:3,1:3]

# check distributions of values in the gt matrices
# gts_by_chrom = function(geno_list) {
#     map_df(1:length(geno_list), function(chrom) {data.frame(table(geno_list[[chrom]]))}, .id = 'chrom') %>%
# 	rename(gt = Var1, n = Freq) %>% spread(gt, n)
# }
# gts_by_chrom(impute_gt)

# summarise changes made during imputation 
smry = imputation_smry(bxd_cross$geno, impute_gt)

# extract the joined genotype tbl
imp_tbl = smry$ori_imp_joined

# we want to revert genotypes for "wrong imputations" where genotype was not originally missing but is missing after imputations
# also want to revert cases when imputation was not needed in the first place, i.e. the genotype was already known
# this is an extra layer of protection from keeping wrong calls
# another similar check is applied during inference of numeric genotypes, which also brings back het and newvariant genotypes
imp_tbl = fix_imp_tbl(imp_tbl)

# keep the fixed genotypes and convert codes
imp_tbl = imp_tbl %>% 
    select(strain, locus, gt.imp.fixed) %>% 
    rename(gt = gt.imp.fixed) %>%
    mutate(gt = recode(gt, '0' = NA_character_, '1' = 'B', '2' = 'D'))

# make a genotype file with imputed founder labels
geno_imp = imp_tbl %>% 
    pivot_wider(id_cols = locus, names_from = strain, values_from = gt) %>%
    separate('locus', c('chr', 'pos', 'end')) %>% 
    mutate_at(c('pos', 'end'), as.integer) %>%
    mutate(chr = fct_relevel(chr, str_c('chr', 1:19))) %>%
    arrange(chr, pos, end)

# blank this map to keep the summary files light
# smry$ori_imp_joined = NULL

# write outputs
saveRDS(geno_imp, '../../data/str_gts/all_foulab_nosegdup_nolowcr_noalshared_padded_imp.rds')
saveRDS(smry, file = '../../data/str_gts/imp_smry.rds')
