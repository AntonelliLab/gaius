devtools::load_all()

wd <- file.path("/Users/djb208/Coding/supersmartR-workshop", 'pipelines',
                '3_supertree')
input_dir <- file.path(wd, '3_align')
tree_file <- file.path(wd, '6_supertree', 'supertree.tre')
alignment_files <- file.path(input_dir, list.files(path = input_dir,
                                                   pattern = '.fasta'))
alignment_names <- names_from_alignments(alignment_files)
tree_names <- names_from_tree(tree_file)
matched_names <- name_match(alignment_names = alignment_names,
                            tree_names = tree_names)
groups <- groups_get(tree_file = tree_file, matched_names = matched_names)
backup_groups <- groups
alignment_list <- alignment_read(flpths = alignment_files)



column_cutoff = .1
tip_cutoff = column_cutoff
min_ntips = 5
min_ngenes = 2
min_nbps = 200
res <- list()

groups <- backup_groups
all_tips <- character(0)
groups <- groups[names(groups) != 'unmatched']
grp_ids <- names(groups)
grp_ids <- c(grp_ids, 'backbone')
for (grp_id in grp_ids) {
  if (grp_id == 'backbone') {
    selected <- backbone_sequences_select(groups = groups,
                                          alignment_list = alignment_list)
  } else {
    nms <- groups[[grp_id]]
    # select sequences from alignments
    selected <- sequences_select(nms = nms, alignment_list = alignment_list)
  }
  # filter selected sequences
  filtered <- sequences_filter(alignment_list = selected,
                               cutoff = column_cutoff, min_nbps = min_nbps)
  if (length(filtered) == 0) {
    next
  }
  # merge into supermatrix
  supermatrix <- supermatrix_get(alignment_list = filtered)
  # drop tips
  supermatrix <- drop_tips(supermatrix = supermatrix, cutoff = tip_cutoff)
  if (length(supermatrix) >= min_ntips &
      length(attr(supermatrix, 'genes')) >= min_ngenes) {
    res[[grp_id]] <- supermatrix
    # record tips
    all_tips <- c(all_tips, names(supermatrix))
  }
}
attr(res, 'tips') <- all_tips
class(res) <- 'supermatrices'
res




