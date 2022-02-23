#' Query vcf
#' 
#' Wrapper around bcftools to query data from vcf. Can read standard vcfs when 
#' `loc_type` is snp or GangSTR vcfs when `loc_type` is str. 
#' Can optionally specify a bed file with regions to only return variants at 
#' at those coordinates or a region in the chr1:32038-34364 format. 
#' If samples are specified, then data for only those samples can be returned.
#' query_type control whether awk is used to parse output of bcftools in long
#' format or if genotypes will be pulled into R in the wide format where
#' each sample is its own column
#' 
#' @param vcf string
#' @param loc_type string
#' @param query_type string
#' @param regions_file string
#' @param samples string
#' @param region string
#' @param verbose logical
#' 
#' @return tibble
#' @export 
#' @examples 
query_vcf = function(vcf, 
		     loc_type = c('str', 'snp')[1], query_type = c('wide', 'long')[1], # wide should be faster in most cases
		     regions_file = NULL, samples = NULL, region = NULL, cstm_filt = NULL, verbose = FALSE) {

    # user out
    if (verbose) print(vcf)

    # form the query string and format it
    if (query_type == 'long') {
	if (loc_type == 'str') {
	    format_string = "'[%CHROM\t%POS\t%INFO/END\t%SAMPLE\t%REPCN\n]'"
	    awk_script='{
		CHROM=$1;POS=$2;EN=$3;SAMPLE=$4;REPCN=$5;
		if(REPCN=="."){REPCN=".,."};
		split(REPCN, RN, ",");
		print CHROM,POS,EN,SAMPLE,RN[1],RN[2];
		}'
	} else {
	    format_string = "'[%CHROM\t%POS\t%POS\t%SAMPLE\t%GT\n]'"
	    awk_script='{
		CHROM=$1;POS=$2;EN=$3;SAMPLE=$4;GTVALS=$5;
		if(GTVALS=="."){GTVALS="./."};
		split(GTVALS, GT, "/");
		print CHROM,POS,EN,SAMPLE,GT[1],GT[2];
		}'
	    awk_script = gsub('\\n', '', awk_script)
	}
	awk_script = gsub('\\n', '', awk_script)
    } else if (query_type == 'wide') {
	if (loc_type == 'str') {
	    format_string = "'%CHROM\t%POS\t%INFO/END\t[%SAMPLE,]\t[%REPCN;]\n'"
	} else {
	    format_string = "'%CHROM\t%POS\t%POS\t[%SAMPLE,]\t[%GT;]\n'"
	}
    }

    # add options
    if (!is.null(cstm_filt)) {
	cmd = str_c("bcftools view ", cstm_filt, " ", vcf, " | bcftools query -f ", format_string, " - ")
    } else {
	cmd = str_c("bcftools query -f ", format_string, " ", vcf)
    }
    if (!is.null(samples)) cmd = str_c(cmd, " -s ", str_c(samples, collapse = ','))
    if (!is.null(regions_file)) cmd = str_c(cmd, " -R ", regions_file)
    if (!is.null(region)) cmd = str_c(cmd, " -r ", region)

    # finish command string
    if (query_type == 'long') cmd = str_c(cmd, " | awk '", awk_script, "' OFS='\t'")
    
    # user out
    if (verbose) cat(cmd) 

    # read data
    if (query_type == 'long') {
	calls = suppressWarnings(read_tsv(pipe(cmd),
	    col_names = c('chr', 'pos', 'end', 'strain', 'GT_A', 'GT_B'),
	    col_types = cols(chr = 'c', strain = 'c', pos = 'i', end = 'i', GT_A = 'i', GT_B = 'i')
	))
	
	# check for zero length data.frame
	if (nrow(calls) == 0) return(NULL)
    } else if (query_type == 'wide') {
	calls = read_tsv(pipe(cmd), 
	    col_names = c('chr', 'pos', 'end', 'strain', 'GT'),
	    col_types = cols(chr = 'c', pos = 'i', end = 'i', strain = 'c', GT = 'c'))
	
	# check for zero length data.frame
	if (nrow(calls) == 0) return(NULL)
    
	# get rid of extra comma (not sure how to do this in bcftools)
	calls = calls %>%
	    mutate_at(c('strain'), function(x) str_replace(x, ',$', '')) %>%
	    mutate_at(c('GT'), function(x) str_replace(x, ';$', ''))
	
	# get sample_names
	strain_names = calls %>% pull(strain) %>% unique %>% str_split(',', simplify = TRUE)
	calls$strain = NULL

	# split genotype column and turn into long format
	calls = calls %>% separate('GT', strain_names, sep = ';')
    }

    # output
    return(calls)
}
