#!/home/momaksimov/anaconda3/envs/r4.2/bin/Rscript

# about: find new str variants by comparing RI to founder genotypes

# clean variables
rm(list = ls())

# options
options(stringsAsFactors = FALSE, dplyr.summarise.inform = FALSE)

# libraries
library(tidyverse)
library(DBI)
library(dbplyr)
library(fs)

# for parallel (deprecated)
# library(furrr)
# library(progressr)

{ ### READ RAW GENOTYPES AND MAKE DENOVO RIL GENOTYPE LIST
    # config
    geno_dir = '../../data/str_gts/'
    gtvals_file = path(geno_dir, 'all_repcn_proc_nosegdup_nolowcr_segreg.rds')
    fou_gtvals_file = path(geno_dir, 'fou_repcn_inf.rds') # inferred founder genotypes

    # load numeric genotype values
    gt_vals          = readRDS(gtvals_file)
    founder_gts_wide = readRDS(fou_gtvals_file)

    # split comma separated genotype into integers
    # `dplyr::separate` is too slow for this
    gt_vals = gt_vals %>%
	pivot_longer(cols = !c(chr, pos, end), names_to = 'strain', values_to = 'gt')
    split_gt = str_split(gt_vals$gt, pattern = ',', simplify = TRUE)
    gt_vals = gt_vals %>% 
	mutate(GT_A = split_gt[,1], GT_B = split_gt[,2]) %>%
	mutate(across(c(GT_A, GT_B), as.integer)) %>%
	select(!gt)

    # get ri only genotypes
    ri_gts = gt_vals %>% filter(!strain %in% c('C57BL', 'DBA'))

    # recode column names
    ri_gts = ri_gts %>% rename(RN_A = GT_A, RN_B = GT_B)
    
    # save number of rows to check against after joins
    rows_gt = nrow(ri_gts)
}

{ ### LOAD INHERITANCE HAPLOBLOCKS
    # determine which parent the denovo variant descended from for each strain by looking at the inheritance blocks

    # load major inheritance blocks, so that founder of origin can be ascertained for denovo genotypes
    haplo_blocks = read_tsv('../../data/inherit_smry/snp_haplo_blocks.tsv', 
			    col_types = cols(strain = 'c', chr = 'c', pos = 'i', end = 'i', lab = 'c', n_loci = 'i'))

    # check haplo blocks
    # haplo_blocks %>% count(lab)
    # haplo_blocks %>% skimr::skim(n_loci)
}

{ ### IDENTIFY DENOVO LOCI
    # identify denovo loci where either A or B allele is unlike all of the founder alleles

    # cast founder_gts wide (make sure to use the inferred ones
    # NOTE: already done for inferred founder genotypes
    # founder_gts_wide = founder_gts %>% 
    #     gather(var, val, RN_A, RN_B) %>%
    #     mutate(var = sub('RN_', '', var)) %>%
    #     unite('var', c('strain', 'var')) %>% 
    #     spread(var, val)

    # compute intensive, so stick to column operations instead of pivoting
    # first join founder genotypes
    denovo_gts = ri_gts %>%
	select(chr, pos, end, strain, matches('_[AB]')) %>%
	# join founder genotypes onto ri genotypes
	left_join(founder_gts_wide, by = c('chr', 'pos', 'end'))

    # essentially this is: compare GT_A to C57BL_[AB],DBA_[AB] and GT_B to C57BL_[AB],DBA_[AB]
    denovo_gts = denovo_gts %>% 
	mutate(across(.cols = matches('^(C57BL|DBA)'),
		      .fns = list(RN_A = ~.x == RN_A, RN_B = ~.x == RN_B),
		      .names = 'lgl_{.fn}_{.col}')) %>%
	mutate(matches_founder = rowSums(across( matches('^lgl_') )) != 0) %>%
	select(!matches('^lgl_')) %>%
	mutate(ri_gt_miss = is.na(RN_A),
	       fou_gt_miss = is.na(C57BL_A) | is.na(DBA_A), 
	       is_het = !ri_gt_miss & RN_A != RN_B) %>%
	# count(ri_gt_miss, fou_gt_miss, across(contains('match_GT')))
	filter(!ri_gt_miss & !fou_gt_miss & !matches_founder) %>%
	# count(is_het)
	filter(!is_het)

    if (0) { # slower rowwise parallel way
	denovo_gts = ri_gts %>%
	    select(chr, pos, end, strain, matches('_[AB]')) %>%
	    # join founder genotypes onto ri genotypes
	    left_join(founder_gts_wide, by = c('chr', 'pos', 'end')) %>%
	    mutate(batch = cut_width(row_number(), width = 1e4, labels = FALSE)) %>%
	    nest(data = !batch)

	plan(multicore, workers = 10)
	# plan(sequential)
	with_progress({
	    p <- progressor(steps = nrow(denovo_gts))
	    denovo_gts = denovo_gts %>%
		mutate(res = future_map(data, function(.x) {
		    p()
		    .x %>%
		    # slice(1:00) %>%
		    rowwise %>%
		    summarise(across(c(RN_A, RN_B), ~any(.x %in% c_across(matches('^(C57BL|DBA)'))), .names = 'match_{col}'), 
			      .groups = 'drop')
	    }))
	})

	# unnest columns
	# label strain/locus comb where wither ri or one of the founders is no-call (can't tell if these are denovo)
	# denovos are where neither ri or founder are  missing and either alleles doesn't match any of the founder alleles
	denovo_gts = denovo_gts %>% 
	    unnest(c(data, res)) %>%
	    mutate(ri_gt_miss = is.na(RN_A),
		   fou_gt_miss = is.na(C57BL_A) | is.na(DBA_A), 
		   is_het = !ri_gt_miss & RN_A != RN_B) %>%
	    # count(ri_gt_miss, fou_gt_miss, across(contains('match_GT')))
	    filter(!ri_gt_miss & !fou_gt_miss & (!match_RN_A | !match_RN_B)) %>%
	    # count(is_het)
	    filter(!is_het)
    }

    # make list of denovo loci
    denovo_loci = denovo_gts %>% distinct(chr, pos, end)
 
    # make a list of all possible strain/locus combinations
    denovo_loci_strains = crossing(denovo_loci %>% select(chr, pos, end), denovo_gts %>% distinct(strain))

    # find assignments for loci which are enclosed
    # should only be one block per locus, because blocks are non-overlapping
    encl_loci = full_join(
	denovo_loci_strains,
	haplo_blocks %>% select(strain, chr, block_pos = pos, block_end = end, lab),
	by = c('strain', 'chr')
    ) %>% filter(pos >= block_pos & end <= block_end)

    # assign nearest block to loci which aren't within a block
    gap_loci = full_join(
	denovo_loci_strains %>% anti_join(encl_loci, by = c('chr', 'pos', 'end', 'strain')),
	haplo_blocks %>% select(strain, chr, block_pos = pos, block_end = end, lab),
	by = c('strain', 'chr')
    ) %>% 
    filter(!is.na(pos)) %>%
    mutate(across(c(pos, end), ~pmin(abs(.x - block_pos), abs(.x - block_end)), .names = 'dist_{col}')) %>%
    # rowwise %>%
    # mutate(across(c(pos, end), ~min(abs(.x - c_across(c(block_pos, block_end)))), .names = 'dist_{col}')) %>%
    # ungroup %>%
    mutate(dist = pmin(dist_pos, dist_end)) %>%
    arrange(dist) %>%
    distinct(chr, pos, end, strain, .keep_all = TRUE)

    # combine encl and gap loci
    founder_of_ori = bind_rows(encl_loci, gap_loci %>% select(chr, pos, end, strain, lab, block_pos, block_end))

    # join "founder of origin" label to denovo genotypes
    denovo_gts = denovo_gts %>%
	left_join(founder_of_ori %>% select(chr, pos, end, strain, founder = lab),
		  by = c('strain', 'chr', 'pos', 'end'))

    # check that no dups
    if (nrow(denovo_gts %>% distinct(strain, chr, pos, end, founder)) != nrow(denovo_gts)) stop('Rows lost or duplicated')
    print('Founder or origin for denovo STRs')
    denovo_gts %>% count(founder) %>% mutate(perc = n*100/sum(n)) %>% print

    # check if any outside of block
    if (nrow(denovo_gts %>% filter(is.na(founder))) > 0) stop('Some denovo loci not assigned founder of origin!')
}

{ ### SAVE DATA
    out_dir = '../../data/denovo_info'; dir_create(out_dir)
    saveRDS(denovo_gts,  path(out_dir, 'denovo_gts.rds'))
    saveRDS(denovo_loci, path(out_dir, 'denovo_loci.rds'))
}
