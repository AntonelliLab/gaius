# LIBS
library(gaius)
library(testthat)

# Vars ----
data_dir <- gaius:::datadir_get()
alignment_dir <- file.path(data_dir, 'alignment.fasta')
alignments_dir <- file.path(data_dir, 'alignments')
alignments_files <- file.path(alignments_dir, list.files(path = alignments_dir,
                                                         pattern = '.fasta'))
# sample for speed
alignments_files <- sample(alignments_files, 5)

# Running ----
context('Testing \'matrix\'')