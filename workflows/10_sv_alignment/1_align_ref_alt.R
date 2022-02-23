#!/home/momaksimov/anaconda3/envs/r4.2/bin/Rscript
#PBS -q hotel
#PBS -N align_ref_alt
#PBS -l nodes=1:ppn=1
#PBS -l walltime=01:00:00
#PBS -o job.o
#PBS -e job.e 
#PBS -V
#PBS -M mikmaksi@gmail.com
#PBS -m abe
#PBS -t 1-10
#PBS -A gymreklab-group

# about: for each structural variant, align the REF and ALT alleles to each other to
# find where the actual insertion/deletion is
# this is necessary because REF was chosen in a somewhat arbitrary way and deletion may
# somewhere within the reference, but not corresponding to the start position in the vcf

# options
options(stringsAsFactors = FALSE, dplyr.summarise.inform = FALSE)

# libraries
library(tidyverse)
library(cowplot)
library(fs)
library(progress)

# set directory
setwd('/projects/ps-gymreklab/mikhail/072321_bxd_mutator_paper/workflows/10_sv_alignment')

# get job array id
batch = as.integer(Sys.getenv('PBS_ARRAYID'))
# batch = 5L

# config
sv_dir = '/projects/ps-gymreklab/resources/datasets/BXD/pangenome/symlinks'
# sv_file = 'chr13.pang.5mb.large.vcf.gz'
sv_file = 'chr13.pang.10mb.large.vcf.gz'

# read precomputed file
cmd = sprintf("bcftools query -H -f '%%CHROM\t%%POS\t%%REF\t%%ALT\n' %s", path(sv_dir, sv_file))
sv_locs = read_tsv(pipe(cmd), col_types = cols(.default = 'c'))

# debug
# bcftools view ../data/vep_annot/bxd_svs.annot.vcf.gz chr13:92211872-92355003 | less -S

# rename columns
sv_locs = sv_locs %>% 
    rename_all(~str_replace(.x, '\\[.*\\]', '') %>% str_replace('# ', '') %>% str_replace(':GT', '')) %>%
    mutate(across(c(POS), as.integer)) %>%
    rename(chr = CHROM, pos = POS)

# remove duplicates
sv_locs = sv_locs %>% distinct

# calculate END position and size of svs
sv_locs = sv_locs %>%
    mutate(across(c(REF, ALT), nchar, .names = '{.col}_len')) %>%
    mutate(diff_len = ALT_len - REF_len,
	   sv_type = case_when(sign(diff_len) == 1 ~ 'INS', 
			       sign(diff_len) == -1 ~ 'DEL',
			       TRUE ~ 'OTHER'),
	   sv_size = abs(diff_len),
	   end = pos + sv_size)

# check types
sv_locs %>% count(sv_type)

# remove "OTHER"
sv_locs = sv_locs %>% filter(sv_type %in% c('DEL', 'INS'))

# cut into batches
sv_locs = sv_locs %>%
    mutate(rid = 1:n(),
	   batch = cut_number(rid, n = 10, label = FALSE, boundary = 1)) %>%
    filter(batch == !!batch) %>%
    select(!c(rid, batch))

# functions for aligning ref and alt sequences
library(Biostrings)
align_seqs = function(ref, alt, sv_type) {
    # ref is the subject
    # alt is the pattern
    # we are searching for the pattern in the subject
    # INS: ALT is longer than the REF
    # DEL: ALT is shorter than the REF
    # we want to align the shorter sequence to the longer one
    #	- INS: pattern in REF
    #	- DEL: pattern is ALT

    # padd sequences on left and right so that alignment edges aren't a problem
    padd = rep('X', 10) %>% str_c(collapse = '')
    ref = str_c(padd, ref, padd)
    alt = str_c(padd, alt, padd)

    # determine which is subject
    if (sv_type == 'DEL') {
	subj = ref
	patt = alt
    } else if (sv_type == 'INS') {
	subj = alt
	patt = ref
    }

    # run alignment
    seq_alig = pairwiseAlignment(pattern = patt, subject = subj, type = 'global')

    # collated insertions and deletions
    dels = deletion(seq_alig) %>% 
	as.data.frame %>% 
	as_tibble %>%
	select(start, end, width)
    ins = insertion(seq_alig) %>%
	as.data.frame %>% 
	as_tibble %>%
	select(start, end, width)
    indels = bind_rows(dels %>% mutate(type = 'del'), ins %>% mutate(type = 'ins'))

    # adjust for padding
    indels = indels %>% mutate(across(c(start, end), ~.x - nchar(padd)))

    # count number of ins/del
    indel_cnts = indels %>% 
	count(type) %>% 
	complete(type = c('ins', 'del'), fill = list(n = 0L)) %>%
	split(.$type) %>% map(~.x %>% pull(n))

    # get edges (this shouldn't be relevant with padding
    edges = tibble(
	subj_st = subject(seq_alig) %>% start(),
	subj_en = subject(seq_alig) %>% end(),
	patt_st = pattern(seq_alig) %>% start(),
    	patt_en = pattern(seq_alig) %>% end()
    )

    # get aligned sequences and remove padding
    alig_seqs = list(
	subj = subject(seq_alig) %>% as.character,
	patt = pattern(seq_alig) %>% as.character
    )
    alig_seqs = alig_seqs %>% map(~.x %>% str_replace_all(padd, ''))
    
    # calculate the sv size
    largest_sv = indels %>% slice_max(n = 1, order_by = dplyr::desc(width), with_ties = FALSE) %>% pull(width)
    total_sv = indels %>% 
	mutate(width = if_else(type == 'del', -width, width)) %>%
	pull(width) %>% sum %>% abs

    # output
    return(list(indels = indels, edges = edges, 
		subj = alig_seqs$subj, 
		patt = alig_seqs$patt, 
		largest_sv = largest_sv, total_sv = total_sv,
		n_del = indel_cnts$del, n_ins = indel_cnts$ins))
}

# for testing
# eg = sv_locs %>% filter(REF_len != 1) %>% filter(sv_type == 'INS') %>% slice_sample(n = 1)
# align_seqs(eg$REF, eg$ALT, eg$sv_type)
# eg = (sv_locs %>% filter(pos == 87491945, end == 87491946))[1,]
# align_seqs(eg$REF, eg$ALT, eg$sv_type)
# eg = (sv_locs %>% filter(pos == 91358269, end == 91358321))[1,]
# align_seqs(eg$REF, eg$ALT, eg$sv_type)
# eg = (sv_locs %>% filter(pos == 87569457, end == 87576508))[1,]
# align_seqs(eg$REF, eg$ALT, eg$sv_type)

# pb <- progress_bar$new(total = sv_locs %>% nrow)
aligs = list()
# i = 2
for (i in 1:nrow(sv_locs)) {
# for (i in 1:10) {
    this_loc = sv_locs[i,]
    seq_alig = try(align_seqs(this_loc$REF, this_loc$ALT, this_loc$sv_type), silent = TRUE)
    if (class(seq_alig) != 'try-error') {
	aligs = c(aligs, list(this_loc %>% mutate(alig_data = list(seq_alig))))
    }
    # pb$tick()
    print(sprintf('Processed %d of %d', i, nrow(sv_locs)))
}
aligs = aligs %>% bind_rows()

# extract single row fields
aligs = aligs %>%
    mutate(
	subj       = map_chr(alig_data, ~.x$subj),
	patt       = map_chr(alig_data, ~.x$patt),
	largest_sv = map_int(alig_data, ~.x$largest_sv),
	total_sv   = map_int(alig_data, ~.x$total_sv),
	n_ins      = map_int(alig_data, ~.x$n_ins),
	n_del      = map_int(alig_data, ~.x$n_del)
    ) %>%
    mutate(alig_tbls = map(alig_data, ~.x[c('indels', 'edges')])) %>%
    select(!alig_data)

# check
# aligs %>% pull(alig_tbls)
# aligs %>% select(chr, pos, sv_type, sv_size, largest_sv, total_sv, n_ins, n_del)

# save data
out_dir = '../../data/sv_alig/refalt_alig'; dir_create(out_dir)
saveRDS(aligs, str_c(out_dir, '/batch_', batch, '.rds'))

# check
# readRDS(str_c(out_dir, '/ref_alt_alig.rds'))
