#' Cacluate STR mutator phenotypes
#' 
#' Function for calculating various STR mutator phenotypes which could be used for QTL mapping
#' 
#' @param denovo_ri_gts TODO
#' @param pheno TODO
#' @param min_pts_per_phe TODO
#' @param max_denovo_strains_per_loc TODO
#' @param ml TODO: motif length
#' @param gtloc_per_strain TODO
#' 
#' @return tibble
#' @export 
#' @examples 
#' pheno_vals = calc_pheno_vals(denovo_ri_gts, 'proportion_expanded', 1, 1000, ml = motif_len)
calc_pheno_vals = function(denovo_ri_gts, 
			   pheno = c('denovo_abundance', 'denovo_perc_abundance', 'proportion_expanded', 
				     'expand_delta_ru', 'contract_delta_ru', 
				     'frac_D_block_denovo', 
				     'denovo_ri_len', 'denovo_fou_len', 
				     'denovo_ri_len_expan', 'denovo_fou_len_expan', 
				     'denovo_ri_len_contr', 'denovo_fou_len_contr'), 
			   min_pts_per_phe = 1, 
			   max_denovo_strains_per_loc = 1000, 
			   ml = c('all', 2, 3, 4, 5, 6), 
			   gtloc_per_strain = NULL, verbose = FALSE) {
    # calculate denovo strains per locus
    strains_per_loc = denovo_ri_gts %>% count(chr, pos, end, name = 'n_strains')

    # keep only loci with a maximum number of denovo strains per locus
    lower_lim = 1; upper_lim = max_denovo_strains_per_loc
    strains_per_loc_filt = strains_per_loc %>%
	filter(n_strains >= lower_lim & n_strains <= upper_lim)
    denovo_ri_gts = denovo_ri_gts %>% semi_join(strains_per_loc_filt, by = c('chr', 'pos', 'end'))

    # user out
    if (verbose) {
	print(sprintf('max_denovo_strains_per_loc: %d\nInput loci: %d\nFiltered loci: %d\n', 
		       max_denovo_strains_per_loc, 
		       nrow(strains_per_loc),
		       nrow(distinct(denovo_ri_gts, chr, pos, end))))
    }

    # calculate metrics which depend on motif_len before (potentially) setting it to "all"
    denovo_ri_gts = denovo_ri_gts %>%
	mutate(
	    denovo_ri_len = RN_T*motif_len, # how long the denovo strs are currently in RI strains
	    denovo_fou_len = founder_rn*motif_len # how long were denovo strs in fou before they mutated
	)

    # set motif to "all" if not calculating phenotype separately for each motif
    if (ml == 'all') {
	if (verbose) cat('Using all motif lengths\n')
	denovo_ri_gts = denovo_ri_gts %>% 
	    mutate(motif_len = 'all')
	# print(denovo_ri_gts %>% distinct(motif_len))
    } else {
	if (verbose) cat(sprintf('Filtering in motif_len: %s\n', ml))
	denovo_ri_gts = denovo_ri_gts %>% filter(motif_len == ml)
	if (verbose) cat(sprintf('Unique motif_len in data.frame: %s\n', denovo_ri_gts %>% pull(motif_len) %>% unique %>% str_c(collapse = ',')))
    }

    # switch b/w phenotype cases
    if (pheno == 'denovo_abundance') {
	# calculate the number of de novo loci per strain
	denovo_per_strain = denovo_ri_gts %>% 
	    count(strain, motif_len, name = 'n_loci') %>%
	    complete(strain = denovo_ri_gts %>% distinct(strain) %>% pull(strain),  # some strains might have zero new variants
		     motif_len = denovo_ri_gts %>% filter(motif_len != 1) %>% distinct(motif_len) %>% pull(motif_len),
		     fill = list(n_loci = as.integer(0)))
	
	# check missing values
	# denovo_per_strain %>% filter(is.na(n_loci))
	
	# select mutator phenotype
	mutator_phenos = denovo_per_strain %>%
	    # mutate(frac = n_loci/tot_loci) %>% # note: 04/23/21 using fraction
	    select(strain, motif_len, pheno = n_loci, n = n_loci) %>%
	    mutate(metric = 'denovo_abundance')

    } else if (pheno == 'denovo_perc_abundance') {
	if (is.null(gtloc_per_strain))
	    stop('If pheno is "denovo_perc_abundance", then need "gtloc_per_strain"')

	# calculate the number of de novo loci per strain
	denovo_per_strain = denovo_ri_gts %>% 
	    count(strain, motif_len, name = 'n_loci') %>%
	    complete(strain = denovo_ri_gts %>% distinct(strain) %>% pull(strain),  # some strains might have zero new variants
		     motif_len = denovo_ri_gts %>% filter(motif_len != 1) %>% distinct(motif_len) %>% pull(motif_len),
		     fill = list(n_loci = as.integer(0)))
	
	# check missing values
	# denovo_per_strain %>% filter(is.na(n_loci))

	# normalize
	denovo_per_strain = denovo_per_strain %>%
	    left_join(gtloc_per_strain, by = 'strain') %>%
	    mutate(perc = n_loci/n_gt) %>%
	    select(-n_gt)
	
	# select mutator phenotype
	mutator_phenos = denovo_per_strain %>%
	    # mutate(frac = n_loci/tot_loci) %>% # note: 04/23/21 using fraction
	    select(strain, motif_len, pheno = perc, n = n_loci) %>%
	    mutate(metric = 'denovo_perc_abundance')

    } else if (pheno == 'proportion_expanded') {
	# calculate the number of expansions/contractions relative to founder (INCLUDE DENOVO LOCI FILTERED WITH DIFFERENT THRESHOLDS)
	n_delta_fou_by_strain = denovo_ri_gts %>%
	    filter(!is.na(founder_rn)) %>%
	    filter(delta_fou != 0) %>%
	    count(strain, motif_len, expand_type, name = 'n_loci') %>%
	    complete(strain = denovo_ri_gts %>% distinct(strain) %>% pull(strain),
		     motif_len = denovo_ri_gts %>% filter(motif_len != 1) %>% distinct(motif_len) %>% pull(motif_len),
		     expand_type,
		     fill = list(n_loci = as.integer(0))) %>%
	    group_by(strain, motif_len) %>%
	    mutate(tot_denovo = sum(n_loci)) %>%
	    ungroup %>%
	    mutate(frac = n_loci/tot_denovo) %>%
	    complete(strain = denovo_ri_gts %>% distinct(strain) %>% pull(strain),  # some strains might have zero new variants
		     motif_len = denovo_ri_gts %>% filter(motif_len != 1) %>% distinct(motif_len) %>% pull(motif_len),
		     expand_type,
		     fill = list(n_loci = as.integer(0), tot_denovo = as.integer(0), frac = NA))
	
	# check missing values
	# n_delta_fou_by_strain %>% filter(is.na(frac))
	
	# select mutator phenotype
	mutator_phenos = n_delta_fou_by_strain %>% 
		filter(expand_type == 'expan') %>% 
		select(strain, motif_len, pheno = frac, n = tot_denovo) %>%
		mutate(metric = 'proportion_expanded')

    } else if (pheno %in% c('expand_delta_ru', 'contract_delta_ru')) {
	# calculate difference relative to founder (INCLUDE DENOVO LOCI FILTERED WITH DIFFERENT THRESHOLDS)
	len_delta_fou_by_strain = denovo_ri_gts %>%
	    filter(!is.na(founder_rn)) %>%
	    filter(delta_fou != 0) %>%
	    group_by(strain, motif_len, expand_type) %>%
	    summarise(avg__delta = mean(delta_fou, na.rm = TRUE),
		      stdev__delta = sd(delta_fou, na.rm = TRUE),
		      n_agg_pts = length(delta_fou), .groups = 'drop') %>%
	    complete(strain = denovo_ri_gts %>% distinct(strain) %>% pull(strain),  # some strains might have zero new variants
		     motif_len = denovo_ri_gts %>% filter(motif_len != 1) %>% distinct(motif_len) %>% pull(motif_len),
		     expand_type,
		     fill = list(NA))

	# select mutator phenotype
	if (pheno == 'expand_delta_ru') {
	    mutator_phenos = len_delta_fou_by_strain %>%
		filter(expand_type == 'expan') %>%
		select(strain, motif_len, pheno = avg__delta, n = n_agg_pts) %>%
		mutate(metric = 'expand_delta_ru')
	} else {
	    mutator_phenos = len_delta_fou_by_strain %>%
		filter(expand_type == 'contr') %>%
		select(strain, motif_len, pheno = avg__delta, n = n_agg_pts) %>%
		mutate(metric = 'contract_delta_ru')
	}

	# check missing values
	# len_delta_fou_by_strain %>% filter(is.na(avg__delta))
    } else if (pheno == 'frac_D_block_denovo') { # this is only relevant for BXD 
	frac_denovo_by_block = denovo_ri_gts %>%
	    distinct(chr, pos, end, strain, motif_len, founder) %>%
	    count(strain, motif_len, founder) %>%
	    complete(strain = denovo_ri_gts %>% distinct(strain) %>% pull(strain),  # some strains might have zero new variants
		     motif_len = denovo_ri_gts %>% filter(motif_len != 1) %>% distinct(motif_len) %>% pull(motif_len),
		     founder,
		     fill = list(n = 0)) %>%
	    group_by(strain, motif_len) %>%
	    mutate(frac = n/sum(n)) %>%
	    ungroup

	# select mutator phenotype
	mutator_phenos = frac_denovo_by_block %>%
	    filter(founder == 'D') %>%
	    # mutate(frac = n_loci/tot_loci) %>% # note: 04/23/21 using fraction
	    select(strain, motif_len, pheno = frac, n = n) %>%
	    mutate(metric = 'frac_D_blck_denovo')

    } else if (pheno %in% c('denovo_ri_len', 'denovo_fou_len', 
			    'denovo_ri_len_expan', 'denovo_fou_len_expan',
			    'denovo_ri_len_contr', 'denovo_fou_len_contr')) {
	# filter
	denovo_len = denovo_ri_gts %>%
	    filter(!is.na(founder_rn)) %>%
	    filter(delta_fou != 0)

	# decide if expand vs. contract is important
	if (pheno %in% c('denovo_ri_len', 'denovo_fou_len')) {
	    denovo_len = denovo_len %>% 
		group_by(strain, motif_len) %>%
		summarise(across(c(denovo_ri_len, denovo_fou_len), mean, na.rm = TRUE),
			  n = length(strain),
			  .groups = 'drop') %>%
		complete(strain = denovo_ri_gts %>% distinct(strain) %>% pull(strain),  # some strains might have zero new variants
			 motif_len = denovo_ri_gts %>% filter(motif_len != 1) %>% distinct(motif_len) %>% pull(motif_len),
			 fill = list(NA))
	    
	    # select mutator phenotype
	    mutator_phenos = denovo_len %>%
		select(strain, motif_len, pheno = !!pheno, n = n) %>%
		mutate(metric = !!pheno)
	} else {
	    denovo_len = denovo_len %>% 
		group_by(strain, motif_len, expand_type) %>%
		summarise(across(c(denovo_ri_len, denovo_fou_len), mean, na.rm = TRUE),
			  n = length(strain),
			  .groups = 'drop') %>%
		complete(strain = denovo_ri_gts %>% distinct(strain) %>% pull(strain),  # some strains might have zero new variants
			 motif_len = denovo_ri_gts %>% filter(motif_len != 1) %>% distinct(motif_len) %>% pull(motif_len),
			 expand_type = c('contr', 'expan'),
			 fill = list(NA))
	
	    # select mutator phenotype
	    et = str_replace(pheno, '.*_(expan|contr)$', '\\1')
	    mutator_phenos = denovo_len %>%
		filter(expand_type == !!et) %>%
		select(strain, motif_len, pheno = !!str_replace(pheno, '_(expan|contr)$', ''), n = n) %>%
		mutate(metric = !!pheno)
	}

    } else if (pheno == 'denovo_len_bytype') {
	denovo_ri_gts %>%
	    filter(!is.na(founder_rn)) %>%
	    filter(delta_fou != 0) %>%
	    group_by(strain, motif_len, expand_type) %>%
	    summarise(across(c(denovo_len, denovo_fou_len), mean, na.rm = TRUE, .names = "avg__{.col}"),
		      .groups = 'drop') %>%
	    complete(strain = denovo_ri_gts %>% distinct(strain) %>% pull(strain),  # some strains might have zero new variants
		     motif_len = denovo_ri_gts %>% filter(motif_len != 1) %>% distinct(motif_len) %>% pull(motif_len),
		     expand_type,
		     fill = list(NA))
	    
	
	
    } else {
	stop("Don't know how to compute requested phenotype")
    }

    # check 
    # pheno_vals %>% group_by(metric) %>% 
    #     skimr::skim(n) %>% 
    #     select(metric, numeric.p0, numeric.p100)

    # filter mutator phenotyps by minimum number of points per phenotype
    n_strains = mutator_phenos %>% pull(strain) %>% unique %>% length
    mutator_phenos = mutator_phenos %>% filter(n >= min_pts_per_phe)
    n_strains_filt = mutator_phenos %>% pull(strain) %>% unique %>% length

    # user out
    if (verbose) {
	cat(sprintf('min_pts_per_phe: %d\nInput strains: %d\nFiltered strains: %d\n', 
		     min_pts_per_phe, 
		     n_strains,
		     n_strains_filt))
    }

    # output
    return(mutator_phenos)
}
