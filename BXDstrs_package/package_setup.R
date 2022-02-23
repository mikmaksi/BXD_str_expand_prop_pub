#!/home/momaksimov/anaconda3/envs/r4/bin/Rscript

# about: set up an R package for BXD data with commonly used functions and data

# libraries
library(devtools)
library(tidyverse)
library(fs)
library(DBI)
library(dbplyr)

# function for adding data description header to R/data.R
write_data_descr = function(cfg) {
    # example
    # cfg = list(
    #     name = 'STR regions,
    #     data_descr = 'a bed file with STR coordiantes',
    #     data_format = 'a data frame with x, y:',
    #     data_fields = list(price = 'price, in US dollars',
    #     		   carat = 'weight of the diamond, in carats'),
    #     data_source = '\\url{http://www.google.com}',
    #     varname = 'str_regions'
    # )

    # add name
    descr = c(cfg$name, '', 
	      (str_wrap(cfg$data_descr, 80) %>% stringr::str_split('\n'))[[1]], '', 
	      str_c('@format ', cfg$data_format), '\\describe{') 

    # add fields
    for (field in names(cfg$data_fields)) {
	descr = c(descr, sprintf('\t\\item{%s}{%s}', field, cfg$data_fields[[field]]))
    }

    # add data source
    descr = c(descr, '}', str_c('@source ', cfg$data_source))

    # add comments
    descr = str_c("#' ", descr)

    # add variable name
    descr = c(descr, str_c('"', cfg$varname, '"'), '')

    # write the lines
    write_lines(descr, 'R/data.R', append = TRUE)
}

# config
package = 'BXDstrs'; dir_create(package)

# clean up
if (dir_exists(path(package, 'man'))) dir_delete(path(package, 'man'))
if (dir_exists(path(package, 'data'))) dir_delete(path(package, 'data'))
if (file_exists(path(package, 'R', 'data.R'))) file_delete(path(package, 'R', 'data.R'))

# set up DESCRIPTION and NAMESPACE files w/ custom fields
package_fields = list(
    Package = package,
    Version = "0.0.0.9000",
    Title = "Mikhail's BXD STR data",
    Description = "Holds data and helper functions",
    `Authors@R` = "person(\"Mikhail\", \"M\", , \"momaksimov@health.ucsd.edu\", c(\"aut\", \"cre\"), comment = c(ORCID = \"YOUR-ORCID-ID\"))",
    License = "`use_mit_license()`, `use_gpl3_license()` or friends to pick a license",
    Encoding = "UTF-8",
    LazyData = "true",
    Roxygen = "list(markdown = TRUE)",
    RoxygenNote = "7.1.1"
)
create_package(package, fields = package_fields)

# debug
# proj_activate(package)

### code ###

# create R file for documenting included datasets
r_file_name = 'data'; full_r_file_name = str_c('R/', r_file_name, '.R')
if (file_exists(full_r_file_name)) file_delete(full_r_file_name)
use_r(r_file_name, open = FALSE)

### data ###

## strain information
{
    # config
    data_dir = '../../info'

    # load strain info
    strain_info = read_tsv(path(data_dir, 'bxd_strain_names_plus.tsv'),
	col_types = cols(
	  long_name      = 'c',
	  short_name     = 'c',
	  bxd_id         = 'c',
	  batch          = 'c',
	  off_epoch      = 'c',
	  type           = 'c',
	  gen_inbreeding = 'i',
	  gen_breeding   = 'i',
	  note           = 'c',
	  cross_type     = 'c'
	)
    )
    
    # remove redundant columns
    strain_info = strain_info %>% select(!c(batch, note, type))

    # get list of strains for which sequencing data is available (either have strain in snp vcf or bams for str genotyping)
    seqed_strains = list(snp = 'snp_strain_list', 
			 str = 'str_strain_list') %>%
	map_df(~read_tsv(path(data_dir, .x), 
			 col_names = 'short_name', 
			 col_types = cols('c')), .id = 'seq_data') %>%
	mutate(is_seq = TRUE) %>%
	pivot_wider(id_cols = short_name, 
		    names_from = seq_data, names_prefix = 'is_seq_',
		    values_from = is_seq, values_fill = list(is_seq = FALSE))
    strain_info = strain_info %>% 
	left_join(seqed_strains, by = 'short_name') %>%
	mutate(across(matches('is_seq_'), ~replace_na(.x, FALSE)))

    # check
    # strain_info %>% count(is_seq_str, is_seq_snp)

    # add to the data cache
    use_data(strain_info, overwrite = TRUE)

    # document 
    cfg = list(
	name = 'BXD strain info',
	data_descr = 'Strain information about BXD strains',
	data_format = 'a data frame with:',
	data_fields = list(
	  long_name      = 'Long name for the strain; e.g. 4512-JFI-0361_BXD001_TyJ',
	  short_name     = 'Short name for the strain; e.g. BXD001_TyJ_0361',
	  bxd_id         = 'BXD id for the strain; e.g. BXD001',
	  off_epoch      = 'Name of official epoch',
	  gen_inbreeding = 'Number of generations of inbreeding',
	  gen_breeding   = 'Number of generations of breeding',
	  cross_type     = 'Cross type'),
	data_source = '',
	varname = 'strain_info'
    )
    write_data_descr(cfg)
}


## regions files
{
    # str_regions.bed: successfully genotypes ones
    regions_files = list(
	segregating_strs = path('../../data/ref', 'str_segreg.tsv'),
	lowcr            = path('../../data/ref', 'str_lowcr.tsv'),
	segdup           = path('../../data/ref', 'str_segdup.tsv'),
	monoallelic      = path('../../data/ref', 'str_zv_hetonly.tsv')
    )

    # load regions
    regions = map(regions_files, 
		  read_tsv, 
		  col_types = c(chr = 'c', pos = 'i', end = 'i'))

    # get list of genotyped strs
    cmd = sprintf("zcat %s | awk '{print $1,$2,$3}' OFS=\"\t\"", '../../data/str_gts/all_repcn_raw.tsv.gz')
    regions$genotyped_strs = read_tsv(pipe(cmd)) %>% set_names(c('chr', 'pos', 'end')) %>% mutate(across(c(pos, end), as.integer))

    # find skipped loci those not processed by GangSTR
    # NOTE: this file not generated in this analysis
    gt_attempt = read_tsv('../../data/ref/str_regions.for_gangstr.bed',
			  col_names = c('chr', 'pos', 'end'), 
			  col_types = c(chr = 'c', pos = 'i', end = 'i'))
    skipped = gt_attempt %>% anti_join(regions$genotyped, by = c('chr', 'pos', 'end'))
    regions$skipped = skipped

    # add to data cache
    use_data(regions, overwrite = TRUE)

    # document
    cfg = list(
	name = 'STR regions',
	data_descr = 'a names list of bed files with STR coordinates',
	data_format = 'a data frame with:',
	data_fields = list(chr = 'chromosome name',
			   pos = 'start position bp',
			   end = 'end position bp'),
	data_source = '',
	varname = 'regions'
    )
    write_data_descr(cfg)
}

## motif composition data
{
    # NOTE: this file not generated in this analysis
    # config
    base_dir = '../../../090520_unified_workflow'
    data_dir = path(base_dir, 'data/ref/motif_comp')

    # load motif info
    motif_info = read_tsv(path(data_dir, 'str_regions_mm10_filt_w_hom.tsv.gz'), 
	col_types = cols(
	    chr = col_character(),
	    pos = 'i',
	    end = 'i',
	    motif_len = 'i',
	    motif = 'c',
	    unq_motif = 'c',
	    unq_motif_len = 'i',
	    A = 'i',
	    C = 'i',
	    G = 'i',
	    T = 'i',
	    canon_motif = 'c',
	    canon_unq_motif = 'c'
	)
    ) 

    # drop not useful columns
    motif_info = motif_info %>% select(!c(unq_motif, unq_motif_len, canon_unq_motif))

    # add to the data cache
    use_data(motif_info, overwrite = TRUE)

    # document 
    cfg = list(
	name = 'STR motif information',
	data_descr = 'Properties of STR motifs',
	data_format = 'a data frame with:',
	data_fields = list(chr = 'chromosome name',
			   pos = 'start position bp',
			   end = 'end position bp',
			   motif_len = 'length of the motif bp',
			   motif = 'motif identity',
			   A = 'number of As',
			   C = 'number of Cs',
			   G = 'number of Gs',
			   T = 'number of Ts',
			   canon_motif = 'canonical motif'),
	data_source = '',
	varname = 'motif_info'
    )
    write_data_descr(cfg)
}

## denovo str lists 
{
    # config
    data_dir = '../../data/denovo_info'

    # load motif info
    denovo_strs = read_tsv(path(data_dir, 'denovo_ri_gts_hom.tsv'), 
	col_types = cols(
	  chr         = col_character(),
	  pos         = 'i',
	  end         = 'i',
	  RN_A        = 'i',
	  RN_B        = 'i',
	  strain      = col_character(),
	  founder     = col_character(),
	  fou_rn  = 'i',
	  delta_fou   = 'i',
	  expand_sign = 'i',
	  expand_type = col_character()
	),
	comment = '#'
    )

    # format
    denovo_strs = denovo_strs %>% 
	mutate(RN_T = RN_A + RN_B) %>%
	rename(founder_rn = fou_rn)

    # check number of loci
    # denovo_strs %>% distinct(chr, pos, end)
    
    # add to the data cache
    use_data(denovo_strs, overwrite = TRUE)

    # document 
    cfg = list(
	name = 'Denovo STRs',
	data_descr = 'RIL genotypes for denovo STRs. Homozygous denovos only',
	data_format = 'a data frame with:',
	data_fields = list(
	  chr         = 'chromosome name',
	  pos         = 'start position',
	  end         = 'end position',
	  RN_A        = 'A allele repeat count',
	  RN_B        = 'B allele repeat count',
	  RN_T        = 'Total repeat count',
	  strain      = 'Strain name',
	  founder     = 'Founder label B/D',
	  founder_rn  = 'Repeat count in founder (single allele)',
	  delta_fou   = 'Absolute difference in repeat count b/w RIL and founder (RU)',
	  expand_sign = 'Direction of difference +1/-1',
	  expand_type = 'Expansion/contraction'),
	data_source = '',
	varname = 'denovo_strs'
    )
    write_data_descr(cfg)
}

## qtl2 formatted data
{
    # config
    data_dir = '../../data/snp_qtl2/gw'

    # load qtl2 formatted objects
    # genotype probabilities, physical map and strain kinship
    snp_probs   = readRDS(path(data_dir, 'probs.rds'))
    snp_pmap    = readRDS(path(data_dir, 'pmap.rds'))
    snp_kinship = readRDS(path(data_dir, 'kinship.rds'))

    # collate into one list
    qtl_data = list(
	snp_probs   = snp_probs, 
	snp_pmap    = snp_pmap, 
	snp_kinship = snp_kinship)
    # names(qtl_data)

    # add to the data cache
    use_data(qtl_data, overwrite = TRUE)

    # document 
    cfg = list(
	name = 'R/qtl2 formatted data',
	data_descr = 'genotype probabilities, physical map and strain kinship using SNPs',
	data_format = 'a list of qtl2 objects:',
	data_fields = '',
	data_source = '',
	varname = 'qtl_data'
    )
    write_data_descr(cfg)
}

## file paths
{
    # sqlite files
    sqlite_files = list(
	expr_db    = '../../../070521_expanded_GN_dbase/data/db_files/bxd_expr_db.sqlite',
	gn206_db   = '../../../070521_expanded_GN_dbase/data/db_files/gn206_db.sqlite'
    )
    sqlite_files = sqlite_files %>% map(path_real)

    # get schemas
    fpath = sqlite_files$expr_db
    db_schemas = map(list(expr_schema = sqlite_files$expr_db),
	function(fpath) {
	    conn = dbConnect(RSQLite::SQLite(), fpath)
	    tbls = dbListTables(conn) %>% set_names(.)
	    tbls = map(tbls, ~dbListFields(conn, .x))
	    dbDisconnect(conn)
	    return(tbls)
	}
    )

    # vcf files
    snp_vcf_dir = '../../../../resources/datasets/BXD'
    str_vcf_dir = '../../../051820_unified_workflow/data/call_sets/strs'
    vcf_files = list(
	bxd_snp_indel = path(snp_vcf_dir, 'david_dropbox', 'Merged_gvcf_files_all_chr_screen_recalibrated_INDEL_variants_99.9_PASSED_variants.recode_in_at_least_20.vcf.gz'),
	bxd_snp_indel_unfilt = path(snp_vcf_dir, 'david_dropbox', 'Merged_gvcf_files_all_chr_screen_recalibrated_INDEL_variants_99.9_PASSED_variants.recode.vcf.gz'),
	bxd_strs = path(str_vcf_dir, 'BXD_variants_MM_strs_highQ_calls_strain_filt.vcf.gz')
    )
    vcf_files = vcf_files %>% map(path_real)

    # fasta files
    data_dir = '../../../../resources/dbase/mouse/mm10'
    fasta_files = list(mm10 = path(data_dir, 'mm10.fa'))
    fasta_files = fasta_files %>% map(path_real)

    # bam files
    base_dir = '../../../090520_unified_workflow'
    data_dir = path(base_dir, 'data')
    bam_files = read_tsv(path(data_dir, 'strain_list', 'bxd_strain_list'), col_names = 'bam_path', col_types = cols())
    bam_files = bam_files %>%
	separate('bam_path', c('sub_dir', 'bam_file'), sep = '/') %>%
	mutate(bam_dir = '../../../../resources/datasets/BXD') %>%
	mutate(long_name = str_replace(bam_file, '_phased_possorted_bam.bam', '')) %>%
	left_join(strain_info %>% select(long_name, bxd_id), by = 'long_name') %>%
	mutate(bam_path = path(bam_dir, sub_dir, bam_file)) %>%
	select(bxd_id, bam_path)
    bam_files = bam_files %>% mutate(bam_path = bam_path %>% path_real)

    # check
    # bam_files %>% 
    #     # slice_sample(n = 1) %>% 
    #     pull(bam_path) %>% file_exists
     
    # add to the data cache
    use_data(sqlite_files, vcf_files, fasta_files, bam_files, db_schemas, overwrite = TRUE)

    # document 
    cfg = list(
	name = 'sqlite databases',
	data_descr = 'Addresses of relevant sqlite databases',
	data_format = 'a list of file paths:',
	data_fields = '',
	data_source = '',
	varname = 'sqlite_files'
    )
    write_data_descr(cfg)
    cfg = list(
	name = 'vcf files',
	data_descr = 'Addresses of relevant vcf files',
	data_format = 'a list of file paths:',
	data_fields = '',
	data_source = '',
	varname = 'vcf_files'
    )
    write_data_descr(cfg)
    cfg = list(
	name = 'fasta files',
	data_descr = 'Addresses of relevant fasta files',
	data_format = 'a list of file paths:',
	data_fields = '',
	data_source = '',
	varname = 'fasta_files'
    )
    write_data_descr(cfg)
    cfg = list(
	name = 'bam files',
	data_descr = 'Addresses of relevant bam files',
	data_format = 'a list of file paths:',
	data_fields = '',
	data_source = '',
	varname = 'bam_files'
    )
    write_data_descr(cfg)
    cfg = list(
	name = 'Database schema',
	data_descr = 'Names of tables and fields within table for relevant sqlite databases',
	data_format = 'a list of lists. Top level item are schema names. Next level are table names and finally field names.',
	data_fields = '',
	data_source = '',
	varname = 'db_schemas'
    )
    write_data_descr(cfg)
}

## GN dataset information
{
    # load dataset info
    gn_info = read_tsv(path('../../info', 'expr_dsets_plus.tsv'),
	col_types = cols(
	  dset         = 'c',
	  tissue       = 'c',
	  tissue_group = 'c',
	  tissue_id    = 'c'
	)
    )
    
    # add to the data cache
    use_data(gn_info, overwrite = TRUE)

    # document 
    cfg = list(
	name = 'GeneNetwork dataset info',
	data_descr = 'Tissue types for each GeneNetwork expression dataset',
	data_format = 'a data frame with:',
	data_fields = list(
	  dset         = 'Dataset id; e.g. GN206',
	  tissue       = 'Tissue type; e.g. eye_retina',
	  tissue_group = 'Tissue group; e.g. eye',
	  tissue_id    = 'Tissue id; unique identifier if multiple datasets with same tissue; e.g. brain_hippocampus_1'),
	data_source = '',
	varname = 'gn_info'
    )
    write_data_descr(cfg)
}

## strains per expression dataset
{
    # all datasets and GN206 separately
    # TODO: consolidate into a single database in the future
    dset_strains = map_df(list(sqlite_files$expr_db, sqlite_files$gn206_db), function(db) {
	# connect to database
	conn = dbConnect(RSQLite::SQLite(), db)
	dset_strains = tbl(conn, sql("select distinct strain_id, dset_id from expr_vals")) %>% collect()
	strains = tbl(conn, 'strains') %>% collect()
	dsets = tbl(conn, 'dset_names') %>% collect()
	dbDisconnect(conn)

	# add strain and dataset names
	dset_strains %>%
	    left_join(strains, by = 'strain_id') %>%
	    left_join(dsets %>% select(dset_id = id, GN), by = 'dset_id') %>%
	    distinct %>%
	    select(GN, strain)
    })

    # add to the data cache
    use_data(dset_strains, overwrite = TRUE)
    
    # document 
    cfg = list(
	name = 'List of BXD strains featured in each expression dataset',
	data_descr = 'Each dataset has its own unique set of strains available',
	data_format = 'a data frame with:',
	data_fields = list(
	    GN     = 'Dataset id; e.g. GN105',
	    strain = 'BXD strain name; e.g. BXD001'
	),
	data_source = '',
	varname = 'dset_strains'
    )
    write_data_descr(cfg)
}

## genotyped loci per strain
{
    # NOTE: that only segregating are used here
    # crosschecked with previous results: avg. difference of -37 genotyped loci per strain due to 93 new lowcr loci
    gt_strs = readRDS('../../data/str_gts/all_repcn_proc_nosegdup_nolowcr_segreg.rds')
    gtloc_per_strain = gt_strs %>%
	summarise(across(.cols = c(matches('BXD'), C57BL, DBA), 
			 .fns = list(ngt = ~sum(!is.na(.x)),
				     nmiss = ~sum(is.na(.x)),
				     ntot = ~length(.x)
				     ),
			 .names = '{col}:{fn}'
			 )) %>%
	pivot_longer(everything(), names_to = c('strain', 'var'), names_sep = ':', values_to = 'val') %>%
	pivot_wider(id_cols = strain, names_from = var, values_from = val)

    # reformat
    gtloc_per_strain = gtloc_per_strain %>% select(strain, n_loci = ngt)
    
    # add to the data cache
    use_data(gtloc_per_strain, overwrite = TRUE)

    # document 
    cfg = list(
	name = 'Genotyped loci per BXD strain',
	data_descr = 'Number of successfully genotyped STR loci per BXD strain',
	data_format = 'a data frame with:',
	data_fields = list(
	    strain = 'Strain id; e.g. BXD001',
	    n_loci = 'Number of loci'
	),
	data_source = '',
	varname = 'gtloc_per_strain'
    )
    write_data_descr(cfg)
}

## inheritance block starts/stops
{
    # config
    # data_dir = path(base_dir, 'data')
    # founder_blocks = readRDS(path(data_dir, 'inherit_smry', 'r_data', 'major_blocks.rds'))
    data_dir = '../../data/inherit_smry'
    founder_blocks = read_tsv(path(data_dir, 'snp_haplo_blocks.tsv'),
	col_types = cols(
	  strain = col_character(),
	  chr    = col_character(),
	  pos    = col_integer(),
	  end    = col_integer(),
	  lab    = col_character(),
	  n_loci = col_integer()
	)
    )

    # add to the data cache
    use_data(founder_blocks, overwrite = TRUE)

    # document 
    cfg = list(
	name = 'Founder blocks from snps',
	data_descr = 'bedgraph stype dataframe with starts and stops of B/D inheritance blocks within BXD',
	data_format = 'a data frame with:',
	data_fields = list(
	    strain = 'Strain id; e.g. BXD001',
	    chr    = 'Chromosome name',
	    pos    = 'Start position bp',
	    end    = 'End position bp',
	    lab    = 'Founder label',
	    n_loci = 'Number of loci'
	),
	data_source = '',
	varname = 'founder_blocks'
    )
    write_data_descr(cfg)
}

### packages ###
use_package('purrr')
use_package('stringr')
use_package('magrittr')
use_package('readr')
use_package('tidyr')
use_package('dplyr')
use_pipe()

### build ###

# update documentation
document()

# eck that the package is valid
check()
build()

# install the package 
# NOTE: during development phase, can simply use load_all()
# withr::with_libpaths('~/anaconda3/envs/r4.2/lib/R/library', install())

# d_env = new.env()
# load('../tests/gemma_input.Rdata', envir = d_env)
