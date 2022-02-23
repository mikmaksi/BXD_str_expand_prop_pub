#' Run GEMMA
#' 
#' Wrapper around GEMMA from R
#' 
#' @param geno
#'     locus	al_A	al_B	BXD001	BXD002	BXD003
#'     locA	C	T		12	5	14
#'     locB	G	A		1	8	11
#' @param pheno
#'     strain	pheno1 pheno2
#'     BXD001	1.3	12.2
#'     BXD003	0.2	13.1
#'     ...
#'     BXD222	2	20
#' @param covar 
#'     strain	intercept	covar1	covar2
#'     BXD001	1	1
#'     BXD003	1	2
#'     ...
#'     BXD222	1	2
#' @param grm matrix
#'     strain BXD001 BXD002 ... BXD222
#'     BXD001 1      0.3	 0.89	 
#'     BXD002
#'     BXD222			 1
#' @param strains vector of strains corresponding to those in geno, pheno and covar to be used.
#'	Needs to be provided even without subselection to reorder the matrices.
#' @param gemma_exec string pointing to gemma executable
#' @param tmp_dir
#' @param tmp_file
#' @param permute logical
#' @param verbose logical
#' @return tibble
#' @export 
#' @examples 
run_gemma = function(geno, pheno, covar, grm, strains, covar_cols,
		     gemma_exec, gemma_options, tmp_dir, tmp_file, permute = FALSE, verbose = FALSE) {

    # check data
    if (!all(c('locus', 'al_A', 'al_B', strains) %in% names(geno)))
	stop('geno needs "locus, al_A, al_B, BXD1, BXD2 ..." columns')
    if (!all(c('strain', 'pheno') %in% names(pheno)))
	stop('pheno needs "strain, pheno" columns')
    if (!all(c('strain', 'intercept', covar_cols) %in% names(covar)))
	stop('covar needs "strain, pheno, $covar_cols" columns')
    if (!all(c('strain', strains) %in% names(grm)))
	stop('grm needs "strain, BXD1, BXD2 ..." columns')

    # enforce column order for geno
    geno = geno[,c('locus', 'al_A', 'al_B', strains)]

    # make strain the rowname
    covar = covar %>% column_to_rownames(var = 'strain')
    pheno = pheno %>% column_to_rownames(var = 'strain')
    grm = grm %>% column_to_rownames(var = 'strain')

    # enforce column order for expr values and phenotypes, unless permute is on, in which case shuffle
    if (!permute) {
	covar = covar[strains, c('intercept', covar_cols)]
	pheno = pheno[strains, c('pheno')]
	grm = grm[strains, strains]
    } else {
	reord_strains = sample(strains, size = length(strains), replace = FALSE)
	covar = covar[reord_strains, c('intercept', covar_cols)]
	pheno = pheno[reord_strains, c('pheno')]
	grm = grm[reord_strains, reord_strains]
    }

    # final check on strains
    if (!all(c(dim(covar)[1], 
	       dim(pheno)[1], 
	       (dim(geno)[2]-3)) == length(strains))) {
	stop('Unequal number of strains detected in GEMMA input')
    }

    # write the genotype file (single locus)
    write.table(geno, paste0(tmp_dir, '/', tmp_file, '.geno'), sep = ',', quote = FALSE, row.names = FALSE, col.names = FALSE)

    # write the phenotype file (many phenotypes)
    write.table(pheno, paste0(tmp_dir, '/', tmp_file, '.pheno'), 
		sep = '\t', quote = FALSE, row.names = FALSE, col.names = FALSE)

    # write the covariates file with intercept
    write.table(covar, paste0(tmp_dir, '/', tmp_file, '.covar'), 
		sep = '\t', quote = FALSE, row.names = FALSE, col.names = FALSE)

    # write the subset grm
    write.table(grm, paste0(tmp_dir, '/', tmp_file, '.grm'), sep = '\t', quote = FALSE, row.names = FALSE, col.names = FALSE)

    # send to GEMMA
    # `-n 1` means first and only phenotype
    # `-r2 1` threshold for filtering based on correlation with covariations; 1 means no filter
    # `-miss 1` X threshold for missingness; 1 means no filter
    if (verbose) print('Starting GEMMA')
    cmd = paste(gemma_exec, 
		'-g', paste0(tmp_dir, '/', tmp_file, '.geno'), 
		'-p', paste0(tmp_dir, '/', tmp_file, '.pheno'), 
		'-c', paste0(tmp_dir, '/', tmp_file, '.covar'), 
		'-k', paste0(tmp_dir, '/', tmp_file, '.grm'), 
		gemma_options, # see global config
		'-outdir', tmp_dir,
		'-o gemma_out')
    
    # debug
    if (verbose) print(cmd)

    ret = system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE, wait = TRUE)
    # ret = system(cmd, ignore.stdout = FALSE, ignore.stderr = FALSE, wait = TRUE)
    if (verbose) print('GEMMA done')
	
    # read data
    gemma_res = try(read.table(paste0(tmp_dir, '/', 'gemma_out.assoc.txt'), header = TRUE, sep = '\t'))

    # check for errors and open outputs
    if (class(gemma_res)[1] == 'try-error') {
	stop('GEMMA failed')
    }

    # clear out temporary files
    tmp_files = paste0(tmp_file, c('.geno', '.pheno', '.covar', '.grm'))
    for (f in paste0(tmp_dir, '/', c(tmp_files, 'gemma_out.assoc.txt', 'gemma_out.log.txt'))) {
        if (file.exists(f)) {
            file.remove(f)
        }
    }

    # output
    return(gemma_res %>% as_tibble %>% select(rs, n_miss, beta, se, p_wald, p_lrt, p_score))
}
