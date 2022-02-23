#!/bin/bash
#PBS -q hotel
#PBS -N vep
#PBS -l nodes=1:ppn=1
#PBS -l walltime=01:00:00
#PBS -o job.o
#PBS -e job.e 
#PBS -V
#PBS -M mikmaksi@gmail.com
#PBS -m abe
#PBS -t 1-10
#PBS -A gymreklab-group

# about:

# cd into main analysis directory
cd /projects/ps-gymreklab/mikhail/072321_bxd_mutator_paper/workflows/7b_variant_effect_analysis_new_wind

# activate conda environment
# source activate vep

# config
# vcf_dirs=(../../data/vep_annot_new_wind ../../data/vep_annot_new_wind ../../data/vep_annot_new_wind)
# vcf_files=(bxd_snp_indel.annot.vcf.gz bxd_strs.annot.vcf.gz bxd_svs.annot.vcf.gz)
vcf_dirs=(../../data/vep_annot_new_wind)
vcf_files=(bxd_svs.annot.vcf.gz)

# make output directory if necessary
out_dir='../../data/vep_annot_new_wind' && mkdir -p $out_dir

# function to run vep
run_vep() {
    # vep command line options
    # OPS='--ccds --uniprot --symbol --numbers --domains --regulatory --canonical --protein --biotype --uniprot --tsl --appris --gene_phenotype --af --af_1kg --af_esp --af_gnomad --max_af --pubmed --variant_class'
    OPS='--everything'
    cache_dir='/projects/ps-gymreklab/resources/dbase/mouse/vep/'
    # threads=14
    threads=1

    # params
    local file=$1
    local out=$2
    local stats_out=$3

    # form command
    cmd="vep -v --cache --offline \
	    --format vcf \
	    --dir $cache_dir \
	    --species mus_musculus \
	    --assembly GRCm38 \
	    --cache_version 102 \
	    -i ${file} \
	    -o ${out} \
	    --stats_file ${stats_out} \
	    --tab \
	    $OPS \
	    --fork $threads \
	    --force_overwrite"
    echo $cmd
}
# cmd=`run_vep test.annot.vcf.gz test.annot.tsv test.stats`
# eval $cmd

# loop over each file that needs to be processed
for ((i = 0; i < ${#vcf_files[@]}; i++))
do
    # run VEP
    vcf=${vcf_files[$i]}
    vcf_dir=${vcf_dirs[$i]}
    cmd=`run_vep ${vcf_dir}/${vcf} ${out_dir}/${vcf%.vcf.gz}.annot.tsv ${out_dir}/${vcf%.vcf.gz}.stats`
    echo $cmd
    eval $cmd
done
