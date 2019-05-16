# Lib
library(gaius)
wd <- 'project'
# Set parameters
pset(max_name_dist = 0.15, column_cutoff = 0.5, tip_cutoff = 0.1,
     min_ntips = 5, max_ntips = 100, min_ngenes = 5, min_nbps = 250)
# Align
alignment_files <- mafft('--auto', sequence_files)
# Identify groups
groups <- groups_get(alignment_files = alignment_files, tree_file = tree_file)
# Generate supermatrix
alignment_list <- alignment_read(flpths = alignment_files)
supermatrices <- supermatrices_get(groups = groups,
                                   alignment_list = alignment_list)
# Phylogeny
specieslevel_files <- phylo(supermatrices = sift(supermatrices,
                                                 drop = 'unmatched'),
                            method = 'astral')
backbone_files <- phylo(supermatrices = sift(supermatrices, keep = 'unmatched'),
                        method = 'raxml')
# Supertree
supertree_file <- supertree(specieslevel = specieslevel_files,
                            backbone = backbone_files)
