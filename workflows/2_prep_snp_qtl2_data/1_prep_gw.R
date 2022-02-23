#!/home/momaksimov/anaconda3/envs/r4.2/bin/Rscript

# about: prep qtl2 compatible data (probabilties and maps) using a gw LD-indep snp marker set

# options
options(stringsAsFactors = FALSE, dplyr.summarise.inform = FALSE)

# libraries
library(tidyverse)
library(fs)
library(readxl)
library(qtl2)

# clean variables
rm(list = ls())

# config
bxd_geno_url = 'http://datafiles.genenetwork.org/download/GN600/BXD_Geno_2017.xlsx'
snp_geno_dir = '../../data/snp_gts/fou_lab'
out_dir = '../../data/snp_qtl2/gw'; dir_create(out_dir)

# download marker data from GeneNetwork
bxd_geno_file = 'BXD_Geno_2017.xlsx'
if (!file_exists(bxd_geno_file)) download.file(url = bxd_geno_url, destfile = 'BXD_Geno_2017.xlsx')

# read genotype file
bxd_geno = read_xlsx(bxd_geno_file, skip = 20, sheet = '2017',
		     col_types = c('text', 'text', 'numeric', 'numeric', 'numeric', rep('text', 202)),
		     .name_repair = 'universal')

# keep only relevant chromosomes
bxd_geno = bxd_geno %>% filter(Chr %in% 1:19)

# rename columns and nest the data
bxd_markers = bxd_geno %>% 
    select(!Mb_mm9) %>%
    rename(chr = Chr, marker = Locus, pos = Mb_mm10, cM = cM_BXD) %>%
    # nest(gts = !c(chr, marker, pos, cM)) %>%
    mutate(chr = str_c('chr', chr), pos = as.integer(pos*1e6), end = pos) %>%
    select(chr, pos, end, marker, cM) %>%
    unite('locus', c(chr, pos, end), sep = '_', remove = FALSE)

# switch genomic coordinates for markers that appear to be out of order
marker_switches = list(
    list(from = c("UNC5348732", "rs6301139"), to = c("rs6301139", "UNC5348732")),
    list(from = c("UNCHS045465", "rs6377403", "UNCHS045463"), to = c("UNCHS045463", "UNCHS045465", "rs6377403")) # all of these
)
bxd_markers_sw = bxd_markers
# i = marker_switches[[1]]
for (i in marker_switches) {
    # bxd_markers_sw[match(i$to, bxd_markers_sw$marker),] %>% arrange(pos) %>% print
    bxd_markers_sw[match(i$to, bxd_markers_sw$marker),c('cM')] = 
	bxd_markers_sw[match(i$from, bxd_markers_sw$marker),c('cM')]
    # bxd_markers_sw[match(i$to, bxd_markers_sw$marker),] %>% arrange(pos) %>% print
}
bxd_markers = bxd_markers_sw

# load founder genotypes for snps
# .x = bxd_markers %>% nest(req_marks = !chr) %>% pull(req_marks) %>% .[[11]]; .y = 'chr11'
geno = bxd_markers %>%
    nest(req_marks = !chr) %>%
    mutate(ret_snps = map2(req_marks, chr, function(.x, .y) {
	# debug
	print(.y)

	# load snp genotype
    	gt_vals = read_tsv(path(snp_geno_dir, str_c(.y, '.tsv.gz')), col_types = cols(.default = 'c'))

	# select those SNP that match identically with marker location
	match_snps = gt_vals %>% semi_join(.x, by = 'locus')

	# make list of markers for which a snp was not found
	nonmatch_mark = .x %>% anti_join(match_snps, by = 'locus')

	# extract genomic position of all loci
	snp_locs = gt_vals %>% select(locus) %>% separate('locus', c('chr', 'pos', 'end'), convert = TRUE)

	# run `bedtools closest` to find nearest SNP for those markers where exact match was not found
	tmp_files = list(a = tempfile(), b = tempfile())
	bed_files = list(a = nonmatch_mark %>% mutate(chr = .y) %>% select(chr, pos ,end) %>% arrange(pos, end),
			 b = snp_locs %>% arrange(pos, end))
	walk2(bed_files, tmp_files, ~write_tsv(.x, .y))
	closest_snps = read_tsv(pipe(sprintf("bedtools closest -a %s -b %s -t first -D ref", tmp_files$a, tmp_files$b)),
				col_names = c('chr', 'pos', 'end', 'chr_near', 'pos_near', 'end_near', 'dist'), 
				col_types = c(chr = 'c', pos = 'i', end = 'i', chr_near = 'c', pos_near = 'i', end_near = 'i', dist = 'i'), 
	)
	walk(tmp_files, file_delete)

	# unite loci definitions
	closest_snps = closest_snps %>% 
	    unite('locus_near', c('chr_near', 'pos_near', 'end_near')) %>%
	    unite('locus', c('chr', 'pos', 'end'))

	# check for missing
	# closest_snps %>% count(is.na(locus), is.na(locus_near))

	# might find the same closest locus for multiple requested loci
	# in this case just keep one
	# closest_snps %>% filter(locus_near %in% locus_near[duplicated(locus_near)])
	closest_snps = closest_snps %>% group_by(locus_near) %>% slice_min(n = 1, order_by = 'dist', with_ties = FALSE) %>% ungroup

	# want to make sure closest_snp is not already used as an exact match snp to avoid duplicates
	closest_snps = closest_snps %>% anti_join(match_snps, by = c('locus_near' = 'locus'))

	# check distances
	# closest_snps %>% skimr::skim(dist)

	# filter out closest snps that are too far
	closest_snps = closest_snps %>% filter(abs(dist) <= 5e5)

	# subset genotypes for the nearest loci
	nonmatch_snps = closest_snps %>% 
	    select(!dist) %>%
	    inner_join(gt_vals, by = c('locus_near' = 'locus'))
    
	# output
	ret = bind_rows(
	    match_snps %>% mutate(locus_near = locus),
	    nonmatch_snps
	) %>% separate('locus', c('chr', 'pos', 'end'), convert = TRUE, remove = FALSE) %>%
	    arrange(pos, end)

	# final check
	# setdiff(ret %>% select(locus), .x %>% select(locus))
	# setdiff(.x %>% select(locus), ret %>% select(locus))
	# ret %>% count(locus == locus_near)

	# add missing markers as an attribute
	attr(ret, 'missing') = .x %>% anti_join(ret, by = 'locus') %>% select(locus, marker)
	# attributes(ret)$missing

	# output
	return(ret)
    }))

# check requested markers vs. returned snps
# we attempt to find the closest snp to a marker if one isn't found at exact location
# in some cases one closest snp is shared between multiple markers and we keep just one of these marker/snp combinations
loc_smry = pmap_df(geno, function(chr, req_marks, ret_snps) {
    tibble(chr = chr, req_loci = nrow(req_marks), ret_loci = nrow(ret_snps))
}) %>% 
    mutate(miss = req_loci - ret_loci)
loc_smry %>% summarise(across(c(req_loci, ret_loci, miss), sum))
# req_loci ret_loci  miss
#     7124     7101   23

# get dropped markers
dropped_marks = geno %>%
    pull(ret_snps) %>%
    map_df(~attributes(.x)$missing) %>%
    separate('locus', c('chr'), extra = 'drop', remove = FALSE)
dropped_marks %>% count(chr)

# unnest the genotypes
geno_mat = geno %>% 
    select(ret_snps) %>% 
    unnest(ret_snps) %>%
    distinct

# check names
# geno_mat %>% names

# check number of loci where nearest neighbor is used
# "locus": location of marker in LD-indep list
# "locus_near": either same as marker location or locatin of nearest snp
geno_mat %>% count(locus == locus_near)

# check
# (geno_mat$locus %in% bxd_markers$locus) %>% mean

if (0) { # check for an alternative to using `interp_map` from qtl2
    # make list of loci for which cM coordinates need to be interpolated
    to_interp = geno_mat %>%
	filter(locus != locus_near) %>%
	select(locus_near) %>%
	separate('locus_near', c('chr', 'pos'), extra = 'drop', convert = TRUE, remove = FALSE)
    interp_comp = list(
	# interp_map from qtl2
	qtl2 = to_interp %>%
	    mutate(pos = pos/1e6, chr = str_replace(chr, 'chr', '')) %>%
	    split(.$chr) %>%
	    map(~(.x$pos %>% set_names(.x$locus_near))) %>%
	    interp_map(bxd_cross$pmap, bxd_cross$gmap) %>%
	    map_df(~tibble(locus = names(.x), cM = .x)),

	# linear fits by chromosome
	lin_fit = { 
	    lm_fits = bxd_markers %>%
		nest(data = !chr) %>%
		mutate(lm_fit = map(data, ~lm(cM ~ pos + 0, data = .x))) %>%
		select(chr, lm_fit)
	    to_interp %>%
		nest(data = !chr) %>%
		left_join(lm_fits, by = 'chr') %>%
		mutate(.pred = map2(data, lm_fit, function(.x, .y) {
		    parsnip::augment(.y, newdata = .x) %>% select(locus = locus_near, cM = .fitted)
		})) %>%
		select(chr, .pred) %>%
		unnest(.pred)
	},

	# interpolation using approxfun
	lin_interp = {
	    interp_fns = bxd_markers %>%
		nest(data = !chr) %>%
		mutate(fn = map(data, ~approxfun(x = .x$pos, y = .x$cM, rule = 2))) %>%
		select(chr, fn)
	    to_interp %>%
		nest(data = !chr) %>%
		left_join(interp_fns, by = 'chr') %>%
		mutate(.pred = map2(data, fn, function(.x, .y) {
		    tibble(locus = .x$locus_near, cM = .y(.x$pos))
		})) %>%
		select(chr, .pred) %>%
		unnest(.pred)
	}
    )
    interp_comp %>% 
	map_df(~.x, .id = 'type') %>%
	pivot_wider(id_cols = locus, names_from = type, values_from = cM) %>%
	column_to_rownames(var = 'locus') %>%
	cor(method = 'pearson', use = 'pairwise.complete.obs')
    interp_comp %>% 
	map_df(~.x, .id = 'type') %>%
	pivot_wider(id_cols = locus, names_from = type, values_from = cM) %>%
	count(is.na(qtl2), is.na(lin_fit), is.na(lin_interp))
    interp_comp %>%
	map_df(~.x, .id = 'type') %>%
	pivot_wider(id_cols = locus, names_from = type, values_from = cM) %>%
	# filter(is.na(lin_interp))
	filter(qtl2 != lin_interp) %>%
	mutate(diff = abs(lin_interp - qtl2)) %>%
	arrange(desc(diff))
    # use `approxfun` for interpolation of maps, not lm fit
}

# interpolate the locations of the closest loci
interp_fns = bxd_markers %>%
    nest(data = !chr) %>%
    mutate(fn = map(data, ~approxfun(x = .x$pos, y = .x$cM, rule = 2))) %>%
    select(chr, fn)
closest_interp = geno_mat %>%
    filter(locus != locus_near) %>%
    select(locus_near) %>%
    separate('locus_near', c('chr', 'pos'), extra = 'drop', convert = TRUE, remove = FALSE) %>%
    nest(data = !chr) %>%
    left_join(interp_fns, by = 'chr') %>%
    mutate(.pred = map2(data, fn, function(.x, .y) {
	tibble(locus = .x$locus_near, cM = .y(.x$pos))
    })) %>%
    select(chr, .pred) %>%
    unnest(.pred)

# update bxd marker list
bxd_markers_up = bxd_markers %>%
    # get rid of the few markers which were dropped
    semi_join(geno_mat, by = c('locus' = 'locus')) %>% 
    left_join(geno_mat %>% select(locus, locus_near), by = 'locus') %>%
    left_join(closest_interp %>% select(locus, cM_near = cM), by = c('locus_near' = 'locus')) %>%
    mutate(cM = if_else(is.na(cM_near), cM, cM_near))

# make physical and genetic maps
geno_map = bxd_markers_up %>% select(chr, marker = locus_near, cM)
phys_map = bxd_markers_up %>% select(chr, marker = locus_near, pos)

# write outputs in qtl2 format
phys_map %>% 
    select(marker, chr, pos) %>%
    mutate(chr = str_replace(chr, 'chr', ''),
	   pos = pos/1e6) %>% 
    write_csv(path(out_dir, 'pmap.csv'))
geno_map %>% 
    select(marker, chr, pos = cM) %>%
    mutate(chr = str_replace(chr, 'chr', '')) %>% 
    write_csv(path(out_dir, 'gmap.csv'))
geno_mat %>%
    rename(marker = locus_near) %>%
    select(!c(locus, chr, pos, end)) %>%
    relocate(marker, .before = 1) %>%
    write_csv(path(out_dir, 'geno.csv'))

# also save the marker to locus position which doesn't get captured within qtl2
write_tsv(bxd_markers_up %>% select(chr, marker, locus, locus_near), path(out_dir, 'snp_phys_map_w_loc.tsv'))

# load the cross object
bxd_cross = read_cross2('gw.json', quiet = FALSE)

# calculate genotype probabilities
# defaults error_prob = 1e-4 and haldane map function
snp_probs = calc_genoprob(bxd_cross, quiet = FALSE)

# calculate kinship leave-one-chromosome-out kinship matrices
kinship = calc_kinship(snp_probs, "loco")

# save rds objects
saveRDS(snp_probs,      file = path(out_dir, 'probs.rds'))
saveRDS(bxd_cross$gmap, file = path(out_dir, 'gmap.rds'))
saveRDS(bxd_cross$pmap, file = path(out_dir, 'pmap.rds'))
saveRDS(kinship,        file = path(out_dir, 'kinship.rds'))
