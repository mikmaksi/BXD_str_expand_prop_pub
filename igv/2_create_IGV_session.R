#!/home/momaksimov/anaconda3/envs/r4/bin/Rscript

# about: make an IGV session xml to stack CC strains in an order of % expanded

# options
options(stringsAsFactors = FALSE, dplyr.summarise.inform = FALSE)

# libraries
library(tidyverse)
library(cowplot)
library(fs)
library(xml2)

### FOR BAM ALIGNMENT

# read other strain information
strain_info = readRDS('../info/strain_info.rds')

# filter for all that are sequenced
strain_info = strain_info %>% filter(is_seq_str)

# loop over epochs
strains_per_epoch = strain_info %>% 
    select(bxd_id, off_epoch) %>%
    group_by(off_epoch) %>%
    summarise(strains = list(bxd_id))
pwalk(strains_per_epoch, function(off_epoch, strains) {
    # create a new doc
    mydoc    = read_xml('IGV_template.xml')
    resources = xml_child(mydoc, 'Resources')
    panel     = xml_child(mydoc, 'Panel')
    data_dir  = "http://localhost:8889/bxd_sym_links/"

    # add Resource nodes
    walk(str_c(strains, '.bam'), function(.x) {
	new_node = as_xml_document(list(Resource = structure(list(), path = str_c(data_dir, .x))))
	xml_add_child(resources, new_node)
    })

    # add Panel nodes
    walk(str_c(strains, '.bam'), function(.x) {
	new_node = as_xml_document(
	    list(Track = structure(
		list(),
		autoScale    = "true",
		clazz        = "org.broad.igv.sam.CoverageTrack",
		color        = "175,175,175",
		colorScale   = "ContinuousColorScale;0.0;10.0;255,255,255;175,175,175",
		fontSize     = "10",
		id           = str_c(data_dir, .x, "_coverage"),
		name         = str_c(.x, " Coverage"),
		snpThreshold = "0.2",
		visible      = "true"
	    ))
	)
	xml_add_child(panel, new_node)
    })

    # save document
    write_xml(mydoc, str_c(off_epoch, '_igv_session.xml'))
})
