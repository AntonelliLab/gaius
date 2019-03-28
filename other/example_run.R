# Lib
library(gaius)
# Set parameters
pset(max_name_dist = 0.15, column_cutoff = 0.5, tip_cutoff = 0.1,
     min_ntips = 5, max_ntips = 100, min_ngenes = 5, min_nbps = 250)
# Align
alignment_files <- align(flpths = list.files("phylotar_results", '.fasta'),
                         method = 'mafft')
# Identify groups
groups <- groups_get(alignment_files = alignment_files, tree_file = tree_file)
# Generate supermatrix
alignment_list <- alignment_read(flpths = alignment_files)
supermatrices <- supermatrices_get(groups = groups,
                                   alignment_list = alignment_list)
# Phylogeny
specieslevel_files <- phylo(supermatrices = filter(supermatrices, '-unmatched'),
                            method = 'astral')
backbone_files <- phylo(supermatrices = filter(supermatrices, '+unmatched'),
                        method = 'raxml')
# Supertree
supertree_file <- supertree(specieslevel = specieslevel_files,
                            backbone = backbone_files)