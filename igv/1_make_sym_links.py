#!/home/momaksimov/anaconda3/envs/r4.2/bin/python3

# about: make symlinks for large bam files to serve up to IGV via http-server

# imports
import os
import csv
from pandas import read_csv

# read strain list
strain_info_file = '../info/bxd_strain_names_plus.tsv'
strain_info = read_csv(strain_info_file, sep = '\t')

# read list of bams
bxd_info_file = '/projects/ps-gymreklab/mikhail/090520_unified_workflow/data/strain_list/bxd_strain_list'
bxd_info = {}
with open(bxd_info_file) as f:
    # line = ['drive1/4512-JFI-0481_BXD19_TyJ_phased_possorted_bam.bam']
    for line in csv.reader(f, delimiter = '\t'):
        long_name = os.path.basename(line[0]).replace('_phased_possorted_bam.bam', '')
        bxd_id = strain_info[strain_info.long_name == long_name]['bxd_id']
        if (bxd_id.size == 1): bxd_id = bxd_id.iat[0]
        bxd_info[bxd_id] = line[0]

# set up directories
base_dir = '/projects/ps-gymreklab/resources/datasets/BXD'
out_dir = 'bxd_sym_links'
os.makedirs(out_dir, exist_ok = True)

# create symlinks
for strain, bam_path in bxd_info.items():
    # make links for '.bam' and '.bai'
    for ft in ['', '.bai']:
        src_file = os.path.join(base_dir, bam_path + ft)
        dest_file = os.path.join(out_dir, strain + '.bam' + ft)
        if not os.path.exists(dest_file):
            os.symlink(src_file, dest_file)

