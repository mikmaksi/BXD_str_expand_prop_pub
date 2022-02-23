#!/bin/bash
#PBS -q hotel
#PBS -N my_script
#PBS -l nodes=1:ppn=1
#PBS -l walltime=01:00:00
#PBS -o job.o
#PBS -e job.e 
#PBS -V
#PBS -M mikmaksi@gmail.com
#PBS -m abe
#PBS -t 1-19
#PBS -A gymreklab-group

# about: generate an ld-pruned snp list using PLINK

# cd into main analysis directory
# cd /projects/ps-gymreklab/mikhail/072321_bxd_mutator_paper/workflows/1_prep_snp_gts

# for running locally
# PBS_ARRAYID=$1

# debug
# PBS_ARRAYID=1

# config
vcf_dir='../../data/snp_gts/vcf'
out_dir='../../data/plink' && mkdir -p $out_dir
chrom=chr${PBS_ARRAYID}

# load PLINK module
module load plink
module load bcftools # needs to happend after plink loading

# NOTE: this is fairly fast so can do as a loop
for chrom in chr{1..19}
do
    echo $chrom

    # convert vcf into PLINK format
    plink --vcf ${vcf_dir}/${chrom}.vcf.gz \
	--double-id \
	--set-missing-var-ids @:# \
	--make-bed \
	--out ${out_dir}/${chrom}

    # calculate sets of LD-independent loci
    # want to prune moderately to make computations easier, but not too aggresively to get back to fully LD-indep set
    wind=10 # kb
    step=1 # doesn't make a whole lot of difference so set to 1
    ld_thresh=0.999 # no difference between 0.999 and 0.9999
    plink --bfile ${out_dir}/${chrom} \
	--indep-pairwise $wind kb $step $ld_thresh \
	--out ${out_dir}/${chrom}
    awk '{print $1,$2,$2}' FS=':' OFS='\t' ${out_dir}/${chrom}.prune.in > tmp && mv tmp ${out_dir}/${chrom}.prune.in
    wc -l ${out_dir}/*prune*
done
