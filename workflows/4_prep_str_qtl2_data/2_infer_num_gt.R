#!/home/momaksimov/anaconda3/envs/r4.2/bin/Rscript

# about: infer founder genotypes
# NOTE: including old code for inferring RI genotypes as well, but this is no longer necessary for this workflow
# NOTE: we are still using imputed RI labels, not just RI labels that can be calculated from comparing fou and RI directly

# clean vars
rm(list = ls())

# options
options(stringsAsFactors = FALSE)

# libraries
library(tidyverse)
library(fs)
library(optparse)

# for parallel (for deprecated code only)
# library(furrr)
# library(progressr)

# function for inferring missing founder genotypes (non-ununanimous RI genotypes)
infer_founder_gts_nonunan = function(ri_gts, founder_gts) {
    # make list of loci where at least one founder genotype is unkown
    found_gt_miss_loci = founder_gts %>%
	filter(is.na(GT_A)) %>%
	distinct(chr, pos, end, strain) %>%
	mutate(strain = recode(strain, C57BL = 'B', DBA = 'D'))

    # check how many missing per founder
    # found_gt_miss_loci %>% count(strain)

    # subselect RI genotypes for the loci where founder genotypes are missing (only strains with inheritance of the founder that is missing)
    found_miss_gts = ri_gts %>%
	semi_join(found_gt_miss_loci, by = c('chr', 'pos' ,'end', 'imp_gt' = 'strain'))

    # for each locus and founder label:
    # count the frequency of observed RI genotypes
    ri_gt_freq = found_miss_gts %>%
	# remove missing genotypes values (uninformative)
	filter(!is.na(GT_A)) %>%
	# filter out het RI genotypes (these are untrustworthy)
	filter(GT_A == GT_B) %>%
	# remove cases when imputed genotype is unkown (i.e. no founders label)
	filter(!is.na(imp_gt)) %>%
	group_by(chr, pos, end, imp_gt) %>%
	count(GT_A, GT_B, name = 'n_strains') %>%
	mutate(n_gts = length(imp_gt),
	       min_af = min(n_strains),
	       maj_af = max(n_strains)) %>%
	ungroup
    
    # check examples
    if (0) {
	td = ri_gt_freq %>% nest(data = !c(chr, pos, end, imp_gt, n_gts, min_af, maj_af))

	# unanimous (easy keep but doesn't discover new denovos)
	td %>% filter(n_gts == 1)

	# two genotypes; single denovo strain; multiple non-denovo strains (can infer)
	# convervative option; allowing for possibility of a singleton denovo
	td %>% filter(n_gts == 2 & min_af == 1 & maj_af > 1) %>% 
	    slice_sample(n = 1) %>% unnest(data)

	# two genotypes; single denovo strain, but also single non-denovo strain (cannot infer)
	# there are only a few of these anyway
	# maj_af should be bigger than min_af, otherwise unreliable
	td %>% filter(n_gts == 2 & min_af == 1 & maj_af == 1) %>% 
	    slice_sample(n = 1) %>% unnest(data)
   
	# two genotypes; more than one denovo strain; single non-denovo strain
	# this doesn't exist
	td %>% filter(n_gts == 2 & min_af > 1 & maj_af == 1)

	# two genotypes; one or more than one denovo strains; multiple non-denovo strains (can infer)
	# less conservative approach
	td %>% filter(n_gts == 2 & min_af >= 1 & maj_af > min_af) %>% 
	    filter(min_af > 1) %>% # for debug
	    arrange(desc(min_af)) %>% # if there is a close to even split on min and maj_al fraction then more likely to be error
	    slice(1:10) %>% 
	    slice_sample(n = 1) %>% unnest(data)

	# allow for more than one denovo genotype; least convervative approach
	# might be too error prone if we're getting 2 different hom denovos
	# this this case maj_af > min_af comparison doesn't make sense
	td %>% filter(n_gts >= 2 & min_af >= 1 & maj_af > min_af) %>% 
	    filter(n_gts == 3) %>% # for debug
	    arrange(desc(min_af)) %>% # if there is a close to even split on min and maj_al fraction then more likely to be error
	    slice(1:10) %>% 
	    slice_sample(n = 1) %>% unnest(data)

	# hot to write this: "if at most one additional repeat genotypes was present differing from the modeal repeat genotype, 
	# the founder was inferred to have the modal genotype." 
    }
    # v0: 0 denovo genotypes; unanimous; can discover new denovo
    # v1: only 1 denovo genotype; only 1 denovo strain; more than one founder strain
    # v2: only 1 denovo genotype; 1 or more denovo strains; more founder strains than denovo strains
    # v3: more than 1 denovo genotype; more than one denovo strain; more founder strains than smallest number of denovo strains
    miss_found_inferred = ri_gt_freq %>% 
	# filter(n_gts == 1) # v0
	# filter(n_gts == 1 | (n_gts == 2 & min_af == 1 & maj_af > 1)) # v1
	filter(n_gts == 1 | (n_gts == 2 & min_af >= 1 & maj_af > min_af)) # v2: selected
	# filter(n_gts == 1 | (n_gts >= 2 & min_af >= 1 & maj_af > min_af)) # v3

    # different way to compute the same thing
    check_compute = miss_found_inferred %>% 
	arrange(chr, pos, end, imp_gt, desc(n_strains)) %>%
	# unite('id', c('chr', 'pos', 'end', 'imp_gt'), remove = FALSE) %>%
	# filter(!duplicated(id)) %>%
	# select(-id)
	distinct(chr, pos, end, imp_gt, .keep_all = TRUE)

    # now pick a single founder genotype with most "votes" from RIL strians
    miss_found_inferred = miss_found_inferred %>%
	# filter(n_gts > 1) %>%
	group_by(chr, pos, end, imp_gt) %>%
	top_n(1, n_strains) %>%
	ungroup
   
    # run the check
    # top_n will keep duplicates, which are a sign of problems and check below will fail
    if (!all_equal(miss_found_inferred, check_compute)) stop('Error detected in founder inference')

    # check
    founder_gts_inf = founder_gts %>%
	pivot_longer(cols = contains('GT_'),
		     names_to = 'allele',
		     values_to = 'gt.ori') %>%
	left_join(miss_found_inferred %>% 
		    mutate(strain = recode(imp_gt, B = 'C57BL', D = 'DBA')) %>%
		    select(-imp_gt) %>%
		    select(chr, pos, end, strain, contains('GT_')) %>%
		    pivot_longer(cols = contains('GT_'),
				 names_to = 'allele',
				 values_to = 'gt.inf'),
		  by = c('chr', 'pos', 'end', 'strain', 'allele'))

    # create a new "fixed" label
    founder_gts_inf = founder_gts_inf %>%
	# if original genotype is not NA, then keep it, otherwise use inferred
	mutate(gt.fix = if_else(!is.na(gt.ori), gt.ori, gt.inf)) 

    # checks
    # founder_gts_inf %>% 
    #     mutate(has_ori = !is.na(gt.ori), has_inf = !is.na(gt.inf), has_fix = !is.na(gt.fix)) %>%
    #     count(has_ori, has_inf, has_fix)
    # should not have any of the below
    # founder_gts_inf %>% 
    #     mutate(has_ori = !is.na(gt.ori)) %>%
    #     filter(has_ori) %>%
    #     filter(gt.ori != gt.fix)

    # cast founder_gts_inf wide to recreate founder_gts
    founder_gts_inf = founder_gts_inf %>%
	select(-gt.ori, -gt.inf) %>%
	pivot_wider(id_cols = c('chr', 'pos', 'end', 'strain'),
		    names_from = 'allele',
		    values_from = 'gt.fix')

    # double check that known founder values didn't get inferred and the perc missing founder values before/after
    if (1) {
	print('Looking for TRUE -> FALSE; FALSE -> TRUE should not happen')
	founder_gts %>%
	    gather(allele, gt, matches('GT_')) %>%
	    left_join(founder_gts_inf %>%
			gather(allele, gt, matches('GT_')),
		      by = c('chr', 'pos', 'end', 'strain', 'allele'), 
		      suffix = c('.ori', '.inf')) %>%
	    # filter(!is.na(gt.ori) & (gt.ori != gt.inf)) # should be none
	    count(is.na(gt.ori), is.na(gt.inf)) %>% 
	    mutate(perc = n*100/sum(n)) %>% print # TRUE -> FALSE should not happen, FALSE -> TRUE hopefully happends
    }

    # check that number of rows didn't get inflated
    if (nrow(founder_gts) != nrow(founder_gts_inf)) {
	stop('Number of rows in founder_gts not constant after founder genotype inference')
    }

    # check that there is a single genotype value per locus/strain combination
    check_mult_values_per_locstrain(founder_gts_inf)

    # output
    return(founder_gts_inf)
}

# in_df = ri_gts
check_mult_values_per_locstrain = function(in_df) {
    # calculate the number of unique genotypes per locus/strain
    vals_per_locstrain = in_df %>%
	unite('gt', GT_A, GT_B, sep = '/') %>%
	group_by(chr, pos, end, strain) %>%
	summarise(n_unq_rnt = length(unique(gt)), .groups = 'drop')
    mult_gt_per_locstrain = vals_per_locstrain %>% filter(n_unq_rnt > 1)
    if (nrow(mult_gt_per_locstrain) != 0) {
	print("Multiple genotypes found for following locus/strain combinations")
	in_df %>%
	    semi_join(mult_gt_per_locstrain, by = c('chr', 'pos', 'end', 'strain')) %>%
	    arrange(chr, pos, end, strain) %>% print
	stop()
    }
}

# config
geno_dir = '../../data/str_gts/'
imp_geno_file = path(geno_dir, 'all_foulab_nosegdup_nolowcr_noalshared_padded_imp.rds')
gt_vals_file = path(geno_dir, 'all_repcn_proc_nosegdup_nolowcr_segreg.rds')
founder_names  = 'C57BL,DBA'
founder_labels = 'B,D'
out_dir = '../../data/str_gts/'; dir_create(out_dir)

# load imputed B/D labels
imp_geno = readRDS(imp_geno_file)

# load numeric genotype values
gt_vals = readRDS(gt_vals_file)

# subset to only segregating strs
imp_geno = imp_geno %>% semi_join(gt_vals, by = c('chr', 'pos', 'end'))

# gather by strain
imp_geno = imp_geno %>% 
    pivot_longer(cols = !c(chr, pos, end), names_to = 'strain', values_to = 'gt')

# check genotype counts
imp_geno %>% count(gt)

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
    left_join(imp_geno %>% rename(imp_gt = gt), by = c('chr', 'pos', 'end', 'strain'))

# infer missing founder genotypes
founder_gts_inf = infer_founder_gts_nonunan(ri_gts, founder_gts)

# cast founder_gts wide (make sure to use the inferred ones)
founder_gts_wide = founder_gts_inf %>% 
    gather(var, val, GT_A, GT_B) %>%
    mutate(var = sub('GT_', '', var)) %>%
    unite('var', c('strain', 'var')) %>% 
    spread(var, val)

# save founder_gts
saveRDS(founder_gts_wide, path(out_dir, 'fou_repcn_inf.rds'))

if (0) { # DEPRECATED
    # old code for inferring RI genotypes, which is not necessary anymore

    # merge founder gts back
    ri_gts = ri_gts %>% left_join(founder_gts_wide, by = c('chr', 'pos', 'end'))

    # debug
    # ri_gts %>% count(!is.na(C57BL_A), !is.na(DBA_A))

    # assign total genotype dose value based on imputed genotypes
    # NOTE: we only want to overwrite values for which original genotype was unknown
    ri_gts = ri_gts %>% 
	mutate(r = row_number()) %>%
	mutate(to_replace = is.na(GT_A) & !is.na(imp_gt)) %>%
	mutate(replace_with = case_when(imp_gt == founder_labs[1] ~ founder_strains[1], 
					imp_gt == founder_labs[2] ~ founder_strains[2],
					TRUE ~ NA_character_)) %>%
	mutate(replace_gt_known = case_when(
		    replace_with == founder_strains[1] ~ !is.na(.data[[paste0(founder_strains[1], '_A')]]),
		    replace_with == founder_strains[2] ~ !is.na(.data[[paste0(founder_strains[2], '_A')]]),
		    TRUE ~ NA))

    # user out
    print('Summary of fraction of loci to replace with imputed values')
    ri_gts %>% count(to_replace, replace_with) %>% 
	mutate(perc = n*100/sum(n))
    # ri_gts %>% filter(is.na(replace_with))

    # how many of these values can be successfully replaced (i.e. D and B genotypes are known)
    print('For what fraction are replacement genotypes known?')
    ri_gts %>% count(to_replace, replace_with, replace_gt_known) %>% 
	mutate(perc = n*100/sum(n)) # %>%
	# filter(to_replace) %>% print

    # identify denovo loci where either A or B allele is unlike all of the founder alleles
    # compute intensive to do, so stick to column operations instead of pivoting
    # split loci into batches and use parallellization
    plan(multisession, workers = 10)
    # plan(sequential)
    denovo_gts = ri_gts %>%
	select(chr, pos, end, strain, matches('_[ab]')) %>%
	mutate(batch = cut_width(row_number(), width = 1e4, labels = FALSE)) %>%
	nest(data = !batch)
    with_progress({
	p <- progressor(steps = nrow(denovo_gts))
	# bonus points for using both `across` and `c_across` with rowwise :)
	# essentially this is: compare GT_A to C57BL_[AB],DBA_[AB] and GT_B to C57BL_[AB],DBA_[AB]
	denovo_gts = denovo_gts %>%
	    mutate(res = future_map(data, function(.x) {
		p()
		.x %>%
		# slice(1:00) %>%
		rowwise %>%
		summarise(across(c(GT_A, GT_B), ~any(.x %in% c_across(matches('(C57BL|DBA)'))), .names = 'match_{col}'), 
			  .groups = 'drop')
	}))
    })

    # unnest columns
    # label strain/locus comb where wither ri or one of the founders is no-call (can't tell if these are denovo)
    # denovos are where neither ri or founder are  missing and either alleles doesn't match any of the founder alleles
    denovo_gts = denovo_gts %>% 
	unnest(c(data, res)) %>%
	mutate(ri_gt_miss = is.na(GT_A),
	       fou_gt_miss = is.na(C57BL_A) | is.na(DBA_A), 
	       is_het = !ri_gt_miss & GT_A != GT_B) %>%
	# count(ri_gt_miss, fou_gt_miss, across(contains('match_GT')))
	filter(!ri_gt_miss & !fou_gt_miss & (!match_GT_A | !match_GT_B)) %>%
	# count(is_het)
	filter(!is_het)

    # make list of denovo loci
    denovo_loci = denovo_gts %>% distinct(chr, pos, end)

    # prevent imputation at denovo loci
    ri_gts = ri_gts %>% 
	left_join(denovo_loci %>% mutate(allow_imp = FALSE),
		  by = c('chr', 'pos', 'end')) %>%
	mutate(allow_imp = replace_na(allow_imp, TRUE))

    # check
    # ri_gts %>% count(allow_imp)
    # ri_gts %>% select(to_replace, replace_with, replace_gt_known)

    # easier to replace values using indeces
    which_replace_A = ri_gts %>% 
	filter(to_replace & 
	       allow_imp &
	       (replace_with == founder_strains[1]) & 
	       !is.na(replace_gt_known) & 
	       replace_gt_known) %>%
	pull(r)
    which_replace_B = ri_gts %>% 
	filter(to_replace & 
	       allow_imp &
	       (replace_with == founder_strains[2]) & 
	       !is.na(replace_gt_known) & 
	       replace_gt_known) %>%
	pull(r)

    # check
    # identical(which_replace_A, which_replace_B)
    # length(which_replace_A)*100/nrow(ri_gts) # % loci to replace
    # length(which_replace_B)*100/nrow(ri_gts) # % loci to replace

    # make new version of data.frame and then replace values
    # for B
    ri_gts_imp = ri_gts
    ri_gts_imp[which_replace_A, 'GT_A'] = ri_gts_imp[which_replace_A, paste0(founder_strains[1], '_A')]
    ri_gts_imp[which_replace_A, 'GT_B'] = ri_gts_imp[which_replace_A, paste0(founder_strains[1], '_B')]
    ri_gts_imp[which_replace_A, 'GT_T'] = ri_gts_imp[which_replace_A, paste0(founder_strains[1], '_T')]

    # for D
    ri_gts_imp[which_replace_B, 'GT_A'] = ri_gts_imp[which_replace_B, paste0(founder_strains[2], '_A')]
    ri_gts_imp[which_replace_B, 'GT_B'] = ri_gts_imp[which_replace_B, paste0(founder_strains[2], '_B')]
    ri_gts_imp[which_replace_B, 'GT_T'] = ri_gts_imp[which_replace_B, paste0(founder_strains[2], '_T')]

    # delete unnecessary columsn
    ri_gts_imp = ri_gts_imp %>% select(-r, -to_replace, -replace_with, -replace_gt_known)

    # check
    if (0) {
	ri_gts_imp %>% slice(which_replace_A) %>% 
	    mutate(check1 = GT_A == C57BL_A, check2 = GT_B == C57BL_B) %>%
	    summarise(mean(check1, na.rm = TRUE), mean(check2, na.rm = TRUE))
	ri_gts_imp %>% slice(which_replace_B) %>% 
	    mutate(check1 = GT_A == DBA_A, check2 = GT_B == DBA_B) %>%
	    summarise(mean(check1, na.rm = TRUE), mean(check2, na.rm = TRUE))
	
	# ori_geno %>% filter(chr == 'chr1', pos == 3016261)
	# imp_geno %>% filter(chr == 'chr1', pos == 3016261)
	# ri_gts %>% filter(chr == 'chr1', pos == 3016261)
	# ri_gts_imp %>% filter(chr == 'chr1', pos == 3016261)
	# ri_gts %>% filter(chr == 'chr1', pos == 3016261) %>% filter(to_replace)
	# ri_gts_imp %>% filter(chr == 'chr1', pos == 3016261) %>% filter(to_replace)
    }

    # check
    print(sprintf('Perc NA value in raw genotype table %0.2f%%', mean(is.na(ri_gts[,'GT_T']))*100))
    print(sprintf('Perc NA value in imputed genotype table %0.2f%%', mean(is.na(ri_gts_imp[,'GT_T']))*100))

    # prepare final table
    imputed_gt_vals = bind_rows(
	ri_gts_imp,
	founder_gts_inf
    ) %>% select(chr, pos, end, strain, GT_A, GT_B, GT_T) %>%
	# join strain id back
	left_join(strains, by = 'strain') %>%
	select(-strain)

    # check
    # imputed_gt_vals %>% distinct(strain_id)

    # final check on number of loci
    cat(sprintf('Number of unique loci in regions: %d
    Number of unique loci in gt_vals from out_db: %d
    Number of unique loci in imputed_gt_vals: %d\n', 
	nrow(distinct(regions)),
	gt_vals %>% distinct(chr, pos, end) %>% nrow,
	imputed_gt_vals %>% distinct(chr, pos, end) %>% nrow))

    # rename columns for strs
    if (loc_type == 'str') {
	imputed_gt_vals = imputed_gt_vals %>% rename(RN_A = GT_A, RN_B = GT_B, RN_T = GT_T)
	denovo_gts = denovo_gts %>% rename(RN_A = GT_A, RN_B = GT_B)
    }

    # write output file
    write_tsv(imputed_gt_vals, str_c(prefix, '.tsv'))
}
