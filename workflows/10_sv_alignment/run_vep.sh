#!/bin/bash

vep -v --format ensembl --cache --offline --dir /projects/ps-gymreklab/resources/dbase/mouse/vep/ --species mus_musculus --assembly GRCm38 --cache_version 102 -i ../../data/sv_alig/vep/vep_in.vcf -o ../../data/sv_alig/vep/vep_out.tsv --stats_file ../../data/sv_alig/vep/stats_out.tsv --tab --everything --force_overwrite
