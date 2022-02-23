#!/bin/bash

# download data for the Attie study on Diversity Outbred mouse crosses 
# doi: https://doi.org/10.1534/genetics.118.300864
# data: https://datadryad.org/stash/dataset/doi:10.5061/dryad.pj105
# ENA: https://www.ebi.ac.uk/ena/browser/view/PRJNA418100?show=reads (RNAseq data, not using for now
# visualizer: https://churchilllab.jax.org/qtlviewer/attie/islets

# download files
out_dir='../../data/diversity_outbred' && mkdir -p $out_dir
# wget https://datadryad.org/stash/downloads/file_stream/62694 -O ${out_dir}/Attie_DO378_eQTL_viewer_v1.Rdata
# wget https://datadryad.org/stash/downloads/file_stream/62695 -O ${out_dir}/README.txt

# just generate symlinks to save space
ln -s -t $out_dir /projects/ps-gymreklab/mikhail/2021/100421_cc_mutator_revisit/data/DO_attie_raw/*
