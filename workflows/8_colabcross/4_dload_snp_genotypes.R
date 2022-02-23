#!/home/momaksimov/anaconda3/envs/r4/bin/Rscript

# about: following the same steps as CC example directly from R/qtl2 developer: https://github.com/rqtl/qtl2data/trunk/CC
# adapting scripts for simplicity and clarity
# download founder inheritance probabilities for CC strains, CC genotypes and founder genotypes
# create crossinfo, covar, .geno files for CC strains and CC founders compaible with R/qtl2 input
# final create a single cc.zip for import to R/qtl2 as a cross

# clean variables
rm(list = ls())

# options
options(stringsAsFactors = FALSE, dplyr.summarise.inform = FALSE)

# libraries
library(tidyverse)
library(cowplot)
library(fs)

# config
base_dir = '../../data/colab_cross/unc_data'; dir_create(base_dir)
zenodo_url = 'https://zenodo.org/record/377036/files/'
fg_url = "https://ndownloader.figshare.com/files/13623080"
cc_strain_info_file = '../../info/cc_strain_info.tsv'

# read the strain info
# NOTE: some funnel codes are missing, but can be inferred from genotype probabilities
cc_strain_info = read_tsv(cc_strain_info_file, comment = '#', col_types = cols())

## download 36-state genotype probabilities

# set up directories
prob_file = str_c(zenodo_url, 'Prob36.zip')
out_dir = str_c(base_dir, '/', 'Prob36'); dir_create(out_dir)
if (nrow(dir_info(out_dir)) == 0) {
    # download zip file if not done already
    cmd = sprintf('wget -P %s %s', base_dir, prob_file)
    # cat(cmd)
    system(cmd, intern = TRUE)

    # unzip
    unzip(str_c(base_dir, '/', 'Prob36.zip'), exdir = out_dir)

    # clean old file
    file_delete(str_c(base_dir, '/', 'Prob36.zip'))
}

# load genotype probabilities
prob_files = as.character(dir_ls(out_dir))
names(prob_files) = sub('.csv.gz', '', basename(prob_files))
probs = map(prob_files, function(f_path) {
    read_csv(f_path, col_type = cols(marker = 'c', chromosome = 'c', 'position(B38)' = 'i', 
				     .default = 'd'))
})

## determine cross orders; use orders in Table S2 when available
# make a matrix from funnel codes
# columns A,B,C ... H do not indicated founders A-H
# rather then indicate order of cross 1-8
# value indicates the number of the letter in the alphabet
# ...
#              ABCDEFGH
# BECADHFG --> 25314867 
funnel_codes_split = cc_strain_info %>%
    group_by(strain) %>%
    group_modify(.f = function(.x, .y) {
	if (is.na(.x$funnel_code)) return(tibble(NULL))
	fc = str_split(.x$funnel_code, '')[[1]]
	# tibble(ord = 1:length(fc), founder = fc)
	tibble(ord = LETTERS[1:8], founder = match(fc, LETTERS[1:8]))
    }) %>%
    ungroup %>%
    # pivot_wider(id_cols = everything(), names_from = 'founder', values_from = 'ord') %>%
    pivot_wider(id_cols = everything(), names_from = 'ord', values_from = 'founder') %>%
    select(strain, LETTERS[1:8])

# check 
# funnel_codes_split

# join split funnel codes back
cross_info = cc_strain_info %>% left_join(funnel_codes_split, by = 'strain')

# check
# cross_info %>% select(m_chrom, y_chrom, strain, LETTERS[1:8])

# when not available, use the values from the table for A and H 
# A: first in 1st pair in 1st cross
# H: second in 4th pair in 1st cross 
cross_info = cross_info %>%
    mutate(m_chrom = if_else(m_chrom %in% c('A/D', 'D/A'), 'A', m_chrom)) %>%
    rowwise %>%
    mutate(A = replace_na(A, match(m_chrom, LETTERS[1:8])),
	   H = replace_na(H, match(y_chrom, LETTERS[1:8]))) %>%
    ungroup 

# problems: CC013, CC023, CC027
# CC013 M = Y = E [genotypes say Y could be B or C, too] --- change Y to B
cross_info[cross_info$strain == "CC013", "H"] <- 2

# this version of file doesn't need founders
cross_info = cross_info %>% filter(!strain %in% LETTERS[1:8])

# check
# cross_info[cross_info$strain == "CC013",] %>% select(LETTERS[1:8])

# for each CC strain calculate the average probabilities inheritance from A, B, C, D ... founder
# across all chrX, chrY, chrM markers
# from which founder strain did the CC strain inherit the X, Y and M chromosomes
mprob = t(sapply(probs, function(a) colMeans(a[a[,2]=="M", paste0(LETTERS, LETTERS)[1:8]])))
yprob = t(sapply(probs, function(a) colMeans(a[a[,2]=="Y", paste0(LETTERS, LETTERS)[1:8]])))
xprob = t(sapply(probs, function(a) colMeans(a[a[,2]=="X", paste0(LETTERS, LETTERS)[1:8]])))

# check cases where mitochondrial genotype is clear
max_mprob = apply(mprob, 1, max)
wh_mprob  = apply(mprob, 1, which.max)
max_yprob = apply(yprob, 1, max)
wh_yprob  = apply(yprob, 1, which.max)

# mitochondrial DNA should be from first strain in the cross
# at least in cases where probability is >90%
stopifnot(all(cross_info[max_mprob>0.9,"A"] == wh_mprob[max_mprob > 0.9]))

# "A" and "H" strains cannot be the same
stopifnot(all(cross_info[,"A"] != cross_info[,"H"]))

# for some strains where "A" strains does not match mitochondrial DNA, average probability should be <50%
stopifnot(all(sort(max_mprob[wh_mprob != cross_info[,"A"]]) < 0.5))

# same as above for "H" strain
stopifnot(all(sort(max_yprob[wh_yprob != cross_info[,"H"]]) < 0.5))

# if 3rd founder is missing, then replace it with the first founder from max X chromosome genotype
wh_xprob = t(apply(xprob, 1, order, decreasing=TRUE))
# head(wh_xprob)
cross_info[is.na(cross_info[,"C"]),"C"] = wh_xprob[is.na(cross_info[,"C"]),1]
stopifnot(all(cross_info %>% pull("A") != cross_info %>% pull("C")))

# TODO: i don't like sampling here because it is not deterministic. Possible solution is to request this information from Churchill lab
if (0) {
    # sort other X chr probs ... put three most probable in the 2nd, 5th, and 6th slots
    # i = 4
    x_ambig_strains_id = which(is.na(cross_info[,"B"]))
    # cross_info[x_ambig_strains,] %>% pull(strain)
    # x_ambig_strains = c("CC004", "CC013", "CC016", "CC019", "CC026", 
    #  			"CC031", "CC032", "CC037", "CC041", "CC042", 
    #  			"CC051", "CC056", "CC059", "CC072")
    # cc_strain_info %>% filter(strain %in% x_ambig_strains)
    for (i in x_ambig_strains_id) {
	# head(wh_xprob)
	# cross_info[i,LETTERS[1:8]]
	xprobs = wh_xprob[i, ]
	known = cross_info[i, c("A","C","H")]
	z = xprobs[!xprobs %in% known]
	# z = (wh_xprob[i, ] %wnin% cross_info[i, c("A","C","H")])
	cross_info[i, c("B","E","F")] = as.list(sample(z[1:3]))
	cross_info[i, c("D","G")] <- as.list(sample(z[4:5]))
    }
} else {
    # instead simply use the coding from KBroman prepared file
    fou_replace = tribble(~strain,~A,~B,~C,~D,~E,~F,~G,~H,
	    'CC004',3,1,8,7,2,6,5,4,
	    'CC013',5,4,3,6,8,1,7,2,
	    'CC016',8,2,6,5,4,3,7,1,
	    'CC019',5,1,2,6,4,8,7,3,
	    'CC026',6,3,4,7,1,8,2,5,
	    'CC031',1,6,5,3,7,2,4,8,
	    'CC032',1,4,2,3,6,5,8,7,
	    'CC037',3,4,8,7,6,1,5,2,
	    'CC041',1,2,4,7,8,3,6,5,
	    'CC042',7,3,5,1,6,8,4,2,
	    'CC051',1,2,4,7,8,5,6,3,
	    'CC056',3,8,1,7,5,4,6,2,
	    'CC059',1,2,4,6,5,8,7,3,
	    'CC072',1,7,2,8,3,4,6,5)
    # strain = 'CC004'
    for (strain in (fou_replace %>% pull(strain))) {
	cross_info[cross_info$strain == strain,LETTERS[1:8]] = 
	    fou_replace[fou_replace$strain == strain,LETTERS[1:8]] 
    }

    # check
    # cross_info %>% filter(strain %in% fou_replace$strain) %>% select(LETTERS[1:8])
}

# further problems:
# CC031 Y chr is B but this is clearly on the X chromosome
# CC037 Y chr is D but this is clearly on the X chromosome
# CC056 Y chr is E but this is clearly on the X chromosome
# ... swap these with one of 2,5,6
cross_info[which(cross_info$strain == "CC031"), LETTERS[1:8]] = as.list(c(1,6,5,3,7,2,4,8))
cross_info[which(cross_info$strain == "CC037"), LETTERS[1:8]] = as.list(c(3,4,8,7,6,1,5,2))
cross_info[which(cross_info$strain == "CC056"), LETTERS[1:8]] = as.list(c(3,8,1,7,5,4,6,2))

# check for missing values
# cross_info %>% 
#     select(strain, LETTERS[1:8]) %>%
#     pivot_longer(cols = LETTERS[1:8], 
# 		 names_to = 'var',
# 		 values_to = 'val') %>%
#     group_by(strain) %>%
#     summarise(miss_vals = sum(is.na(val))) %>%
#     count(miss_vals)

# make directory for all CC related R/qtl2 files from which a zip will be created
cc_dir = str_c(base_dir, '/', 'CC'); dir_create(cc_dir)

# write cross info in rqtl2 format
out_file = str_c(cc_dir, '/', 'cc_crossinfo.csv')
comment_line = paste("# Cross information for Collaborative Cross (CC) lines inferred from",
                        "genotypes from Srivastava et al. (2017)",
                        "doi:10.1534/genetics.116.198838,",
                        "data at doi:10.5281/zenodo.377036", 
			"\n# nrow", nrow(cross_info),
			"\n# ncol", 9)
write_lines(comment_line, out_file)
write_csv(cross_info %>% select(id = strain, LETTERS[1:8]),
	  out_file, append = TRUE, col_names = TRUE)

# make covar dataframe
# NOTE: use cc_strain_info and not cross info, b/c we replaced "A/D" with "A" in mitochondria column
covar = cc_strain_info %>%
    select(id = strain, 
	   mitochondria = m_chrom,
	   Ychr = y_chrom,
	   n_founders,
	   gen_inbreeding,
	   origin)

# write covar  info in rqtl2 format
out_file = str_c(cc_dir, '/', 'cc_covar.csv')
comment_line = paste("# Covariate information for Collaborative Cross (CC) lines inferred from",
                        "genotypes from Srivastava et al. (2017)",
                        "doi:10.1534/genetics.116.198838,",
                        "data at doi:10.5281/zenodo.377036", 
			"\n# nrow", nrow(covar),
			"\n# ncol", ncol(covar))
write_lines(comment_line, out_file)
write_csv(covar, out_file, append = TRUE, col_names = TRUE)

# download genotypes.zip for CC strains
genotypes_file = 'genotypes.zip'
out_dir = str_c(base_dir, '/', 'cc_genotypes'); dir_create(out_dir)
if (nrow(dir_info(out_dir)) == 0) {
    # download zip file if not done already
    cmd = sprintf('wget -P %s %s/%s', base_dir, zenodo_url, genotypes_file)
    # cat(cmd)
    system(cmd, intern = TRUE)

    # unzip
    unzip(str_c(base_dir, '/', genotypes_file), exdir = out_dir)

    # clean old file
    file_delete(str_c(base_dir, '/', genotypes_file))
}

# download founder genotypes from figshare
# https://doi.org/10.6084/m9.figshare.5404762.v2
genotypes_file = 'MMnGM_processed_files.zip'
out_dir = str_c(base_dir, '/', 'fou_genotypes'); dir_create(out_dir)
if (nrow(dir_info(out_dir)) == 0) {
    # download zip file if not done already
    # cmd = sprintf('wget -P %s %s', base_dir, fg_url)
    cmd = sprintf('wget -O %s/%s %s', base_dir, genotypes_file, fg_url)
    # cat(cmd)
    system(cmd, intern = TRUE)

    # unzip
    unzip(str_c(base_dir, '/', genotypes_file), exdir = out_dir)

    # move files one directory up
    fou_files = dir_ls(str_c(out_dir, '/', 'MMnGM')) %>% basename
    file_move(str_c(out_dir, '/MMnGM/', fou_files), 
	      str_c(out_dir, '/', fou_files))

    # clean old file
    file_delete(str_c(base_dir, '/', genotypes_file))
}

# read CC strain genotypes
# NOTE: using SEQgenotypes.csv (135K)
cc_geno = read_csv(str_c(base_dir, '/cc_genotypes/', 'SEQgenotypes.csv'), 
	 col_types = cols(marker = 'c', chromosome = 'c',
			  'position(b38)' = 'i', .default = 'c'))
cc_geno_long = cc_geno %>% 
    pivot_longer(cols = matches('^CC'), 
		 names_to = 'strain',
		 values_to = 'gt') %>%
    rename(pos = 'position(b38)')
cc_geno_long = cc_geno_long %>%
    mutate(gt = if_else(gt %in% c('N', 'H'), NA_character_, gt))
cc_geno_long = cc_geno_long %>%
    left_join(cc_geno_long %>%
		distinct(strain) %>%
		separate('strain', c('strain_name', 'desig'), '/', remove = FALSE), by = 'strain') %>%
    select(marker, chr = chromosome, pos, strain = strain_name, gt)

# load the allele codes
fga = read_csv(file.path(base_dir, 'fou_genotypes', 'MMnGM_allelecodes.csv'),
	       comment = '#', col_types = cols(.default = 'c'))
fga = fga %>% filter(chr %in% c(1:19,"X"))

# omit markers not in founder genotypes (reciprocally)
cc_geno_long = cc_geno_long %>%
    semi_join(fga %>% distinct(marker), by = 'marker')
fga = fga %>%
    semi_join(cc_geno_long %>% distinct(marker), by = 'marker')

# check
# setdiff(cc_geno_long %>% pull(marker) %>% unique, fga %>% pull(marker) %>% unique)
# setdiff(fga %>% pull(marker) %>% unique, cc_geno_long %>% pull(marker) %>% unique)

# reorder rows of genotype data to match fga
cc_geno = cc_geno_long %>% 
    pivot_wider(id_cols = 'marker', names_from = 'strain', values_from = 'gt')
fga_mat = fga %>% column_to_rownames(var = 'marker') %>% as.matrix
# head(fga_mat)
cc_geno = cc_geno[match(cc_geno$marker, rownames(fga_mat)),]

# check order
# all(cc_geno$marker == rownames(fga_mat))

# cut down to just the genotypes
cc_geno_mat = cc_geno %>% column_to_rownames(var = 'marker') %>% as.matrix
# head(cc_geno_df)

# encode genotypes
cc_geno_mat_enc = qtl2convert::encode_geno(cc_geno_mat, fga_mat[,c("A","B")], cores=0)
# head(cc_geno_mat_enc)

# omit markers with no data
cc_geno_mat_enc = cc_geno_mat_enc[rowSums(cc_geno_mat_enc!="-") > 0, , drop=FALSE]

# write genotypes to files, one chromosome at a time
# chr = 1
for(chr in c(1:19,"X")) {
    out_file = str_c(cc_dir, '/cc_geno', chr, ".csv")
    geno_sub = cc_geno_mat_enc[rownames(cc_geno_mat_enc) %in% rownames(fga_mat)[fga_mat[,'chr']==chr], , drop=FALSE]
    
    # check
    # cc_geno_mat_enc %>% dim
    # geno_sub %>% dim

    # set up comment line
    comment_line = paste0("# Chromosome ", chr, " genotypes ",
                        "for Collaborative Cross (CC) lines inferred from ",
                        "genotypes from Srivastava et al. (2017) ",
                        "doi:10.1534/genetics.116.198838, ",
                        "data at doi:10.5281/zenodo.377036",
			"\n# nrow ", nrow(geno_sub),
			"\n# ncol ", ncol(geno_sub) + 1)
    write_lines(comment_line, out_file)
    write_csv(geno_sub %>% as.data.frame %>% rownames_to_column(var = 'marker'),
	      out_file, append = TRUE, col_names = TRUE)
}

# create JSON file for R/qtl2
chroms = c(1:19,"X")
qtl2::write_control_file(str_c(cc_dir, '/', 'cc.json'), crosstype="risib8",
                   geno_file = paste0("cc_geno", chroms, ".csv"),
                   # founder_geno_file = paste0("MMnGM/MMnGM_foundergeno", chroms, ".csv"),
                   # gmap_file = paste0("MMnGM/MMnGM_gmap", chroms, ".csv"),
                   # pmap_file = paste0("MMnGM/MMnGM_pmap", chroms, ".csv"),
                   founder_geno_file = paste0("MMnGM_foundergeno", chroms, ".csv"),
                   gmap_file = paste0("MMnGM_gmap", chroms, ".csv"),
                   pmap_file = paste0("MMnGM_pmap", chroms, ".csv"),
                   covar_file = "cc_covar.csv",
                   crossinfo_file = "cc_crossinfo.csv",
                   geno_codes = c("A"=1, "B"=3),
                   alleles = LETTERS[1:8],
                   xchr = "X",
                   geno_transposed = TRUE,
                   founder_geno_transposed = TRUE,
                   description = paste("Data for Collaborative Cross (CC) lines",
                                     "from Srivastava et al. (2017)",
                                     "doi:10.1534/genetics.116.198838,",
                                     "data at doi:10.5281/zenodo.377036"),
                   overwrite = TRUE)

# create zip file
files_to_zip = c(
    dir_ls(cc_dir) %>% as.character, 
    dir_ls(str_c(base_dir, '/fou_genotypes')) %>% as.character)
zip_name = str_c(base_dir, '/cc.zip')
if (file_exists(zip_name)) file_delete(zip_name)
zip(zip_name, files_to_zip, extras = '-j', zip = 'zip')

# try to open zip file with R/qtl2
cc_cross = qtl2::read_cross2(str_c(base_dir, '/', 'cc.zip'), quiet = FALSE)
