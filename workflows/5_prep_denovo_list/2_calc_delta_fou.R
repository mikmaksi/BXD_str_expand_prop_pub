#!/home/momaksimov/anaconda3/envs/r/bin/Rscript

# about: put mutator phenotypes into an expression database

# clean variables
rm(list = ls())

# options
options(stringsAsFactors = FALSE, dplyr.summarise.inform = FALSE)

# libraries
library(tidyverse)
library(fs)

# gts_df = denovo_ri_gts; col_name = 'foun_rn'
find_max_delta = function(gts_df, col_name) {
    # pivot longer by allele
    max_delta = gts_df %>%
	select(chr, pos, end, strain, matches('_[AB]'), .data[[col_name]]) %>%
	pivot_longer(cols = matches('_[AB]'), names_to = 'al', values_to = 'rn') %>%
	mutate(delta_ru = rn - .data[[col_name]])

    # arrange so that max delta_ru is on top
    max_delta = max_delta %>% arrange(chr, pos, end, strain, desc(abs(delta_ru)))

    # make a unique id
    max_delta = max_delta %>% unite('id', c('strain', 'chr', 'pos', 'end'), remove = FALSE) 

    # keep unique
    max_delta = max_delta[!duplicated(max_delta$id),]

    # output
    return(max_delta %>% distinct(strain, chr, pos, end, delta_ru))
}

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
    
    # pivot founder genotypes back
    founder_gts = founder_gts_wide %>%
	pivot_longer(cols = !c(chr, pos, end), names_to = c('strain', 'allele'), names_sep = '_', values_to = 'gt') %>%
	pivot_wider(id_cols = c(chr, pos, end, strain), names_from = 'allele', values_from = 'gt', names_prefix = 'RN_') %>%
	mutate(strain = recode(strain, 'C57BL' = 'B', 'DBA' = 'D'))
}


{ ### READ DENOVO GENOTYPES
    # config
    data_dir = '../../data/denovo_info'
    denovo_gts = readRDS(path(data_dir, 'denovo_gts.rds'))
    denovo_gts = denovo_gts %>% select(chr, pos, end, strain, RN_A, RN_B, founder)
}

{ ### JOIN NUMERIC FOUNDER GENOTYPES TO DENOVO RI GENOTYPES
    # check that there are no hets among the denovo genotypes
    if (!all(denovo_gts %>% mutate(check = RN_A == RN_B) %>% pull(check)))
	stop('Some het denovo variants detected')
	
    # get list of denovo loci where founders are het
    # only one of the founder may be het so keep strain column
    denovo_fou_het_loci = founder_gts %>%
	filter(RN_A != RN_B) %>%
	distinct(chr, pos, end, strain) %>%
	semi_join(denovo_gts %>% distinct(chr, pos, end), by = c('chr', 'pos', 'end'))

    # assign label to denovo loci where founder is het
    denovo_gts = denovo_gts %>%
	left_join(denovo_fou_het_loci %>% mutate(fou_het = TRUE), by = c('chr', 'pos', 'end', 'founder' = 'strain')) %>%
	mutate(fou_het = replace_na(fou_het, FALSE))

    # check hets
    denovo_gts %>% count(fou_het)
   
    # join founder genotypes to RIL genotypes
    denovo_gts = denovo_gts %>%
	# can divide by 2 because founders are hom (otherwise unassig)
	left_join(founder_gts %>% 
		    unite('fou_gt', c('RN_A', 'RN_B'), sep = '/', remove = FALSE) %>%
		    mutate(RN_T = RN_A + RN_B) %>%
		    mutate(fou_rn = ceiling(RN_T/2)) %>%
		    select(chr, pos, end, strain, fou_gt, fou_rn),
		  by = c('founder' = 'strain', 'chr' = 'chr', 'pos' = 'pos', 'end' = 'end'))
}

{ ### CALCULATE DELTA RN RELATIVE TO FOUNDER FOR DENOVO LOCI
    # calculate max difference in repeat number between RIL and founder strains
    max_delta_fou = find_max_delta(denovo_gts, 'fou_rn')

    # join RN delta relative to founder back to denovo RIL genotypes
    denovo_gts = denovo_gts %>% 
	left_join(max_delta_fou, by = c('strain', 'chr', 'pos', 'end')) %>%
	mutate(expand_sign = sign(delta_ru)) %>%
	mutate(expand_type = case_when(expand_sign == -1 ~'contr', 
				       expand_sign == 1 ~ 'expan',
				       expand_sign == 0 ~ 'equal')) %>%
	mutate(delta_ru = abs(delta_ru)) %>%
	rename(delta_fou = delta_ru)

    # NOTE: max_delta_fou produces the same result as taking RN_A - fou_rn for homozygous loci (both RI and fou)
    # using this function is more general since it would also produce a meaningful value for heterozygous new variant if 
    # we ever wanted to consider these
    if (0) {
	denovo_gts %>% 
	    # filter(fou_het) %>%
	    mutate(check = RN_A - fou_rn) %>%
	    # slice_sample(n = 10)
	    count(check == delta_fou*expand_sign)
    }

    # check locus counts
    # denovo_gts %>% count(chr, pos, end)

    # remove loci where founders are het from the "bycomp" list
    denovo_gts = denovo_gts %>% filter(!fou_het)

    # check counts again
    denovo_gts %>% count(fou_het, expand_type) %>% mutate(perc = n*100/sum(n))
    
    # save denovo_ri_gts for reporting scripts
    out_dir = '../../data/denovo_info'; dir_create(out_dir)
    write_tsv(denovo_gts, path(out_dir, 'denovo_ri_gts_hom.tsv'))
}
