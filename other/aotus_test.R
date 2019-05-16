# Libs ----
library(phylotaR)

# Vars ----
outdir <- file.path('demos', 'aotus')

# Data ----
data('aotus')

# Parse ----
all_clusters <- aotus
cids <- all_clusters@cids
n_taxa <- get_ntaxa(phylota = all_clusters, cid = cids)
keep <- cids[n_taxa > 6]
selected <- drop_clstrs(phylota = all_clusters, cid = keep)
smmry <- summary(selected)
reduced <- drop_by_rank(phylota = selected, rnk = 'species', n = 1)

# Write out ----
for (i in seq_len(nrow(smmry))) {
  cid <- smmry[i, 'ID']
  txids <- get_txids(phylota = reduced, cid = cid, rnk = 'species')
  scientific_names <- get_tx_slot(phylota = reduced, txid = txids,
                                  slt_nm = 'scnm')
  scientific_names <- gsub('\\.', '', scientific_names)
  scientific_names <- gsub('\\s+', '_', scientific_names)
  sids <- reduced@clstrs[[cid]]@sids
  gene_nm <- sub(pattern = '\\s.*$', replacement = '', smmry[i, 'Feature'])
  write_sqs(phylota = reduced, sid = sids, sq_nm = scientific_names,
            outfile = file.path(outdir, paste0(gene_nm, '.fasta')))
}


# Align
library(gaius)
sequence_files = file.path(outdir, list.files(outdir, '.fasta'))
alignment_files <- align(sequence_files = sequence_files, method = 'mafft')
alignment_list <- alignment_read(flpths = alignment_files)

# Supermatrices
# Feedback:
# - annoying way of setting parameters
# - feedback from supermatrices would be helpful
pget()
tree <- ape::stree(n = 6)
tree$tip.label <- c('Aotus_nancymaaec', 'Aotus_azarai', 'Aotus_griseimembra',
                    'Aotus_trivirgatus', 'Aotus_nigriceps', 'Aotus_vociferans',
                    'Aotus_lemurinus')
ape::write.tree(tree, file = file.path(outdir, 'tree.tre'))
groups <- groups_get(alignment_files = alignment_files,
                     tree_file = file.path(outdir, 'tree.tre'))
alignment_list <- alignment_read(flpths = alignment_files)

pset(parameter = c('min_nbps', 'min_ntips', 'min_ngenes',
                   'column_cutoff'),
     val = c(50, 3, 1, 0.9))
supermatrices <- supermatrices_get(groups = groups, min_ntips = 3,
                                   min_ngenes = 1,
                                   alignment_list = alignment_list)

system(paste0('raxmlHPC -f a -m GTRGAMMA -T 2 -# 100 -p ',
              sample(0:10000000, 1), ' -x ', sample(0:10000000, 1),
              ' -n primates -s ', inpt, ' -q ', prttnfl))

raxml <- function(alignment_obj, partition = TRUE, bootstrap = 100L,
                  consensus = TRUE) {
  
}

raxml('-f', 'a', '-m', 'GTRGAMMA', -T 2 -# 100 -p ',
        sample(0:10000000, 1), ' -x ', sample(0:10000000, 1),
      ' -n primates -s ', inpt, ' -q ', prttnfl))

backbone_files <- phylo(supermatrices = sift(supermatrices, keep = 'unmatched'),
                        method = 'raxml')
