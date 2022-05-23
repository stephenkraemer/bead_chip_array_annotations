# ** Imports

import os
# ** Project paths

project_dir = (
    "/omics/groups/OE0029/internal/kraemers/projects/epic-arrays-for-hematopoiesis"
)

# ** Gencode

gencode_download_url = "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.annotation.gtf.gz"

# unmodified gencode gtf (note: with chr prefix)

gencode_gtf = project_dir + "/gencode.vM25.annotation.gtf.gz"

# gencode filtered for principal transcripts of protein coding genes, note that that chromosome prefix ('chr') is removed in this file

gencode_coding_canonical_gtf = (
    project_dir + "/gencode.vM25.annotation_coding_canonical.gtf.gz"
)

# ** Probes

# path to original probe file

original_probes_bed = project_dir + "/2021-03-22_mmbc_probes.bed"

# !head {original_probes_bed}

# path to reformatted probe file

reformatted_probes_bed = project_dir + "/2021-03-22_mmbc_probes_reformatted.bed"

# full Illumina probe file

illumina_probes_url = "https://support.illumina.com/content/dam/illumina-support/documents/downloads/productfiles/mouse-methylation/Infinium%20Mouse%20Methylation%20v1.0%20A1%20GS%20Manifest%20File.csv"

illumina_probes_csv = (
    project_dir
    + "/Infinium_20Mouse_20Methylation_20v1.0_20A1_20GS_20Manifest_20File.csv"
)

# full Illumina probe file BED coordinates

illumina_coordinate_bed = project_dir + "/illumina-all-probes.bed"

# ** Gene annotation results

# *** gtfanno results

custom_intervals_results_dir = project_dir + "/custom-intervals_1500-500"
os.makedirs(custom_intervals_results_dir, exist_ok=True)

custom_intervals_trunk_path = (
    custom_intervals_results_dir + "/custom-intervals_1500-500"
)
custom_intervals_trunk_path

custom_intervals_results_paths_d = dict(
    primary_annos_bed=custom_intervals_trunk_path + "_primary-annotations.bed",
    primary_annos_p=custom_intervals_trunk_path + "_primary-annotations.p",
    all_annos_bed=custom_intervals_trunk_path + "_all-annotations.bed",
    all_annos_p=custom_intervals_trunk_path + "_all-annotations.p",
)

# *** Final tables

gene_annos_primary_one_row = project_dir + "/gene-annos_primary_one-row.bed"
gene_annos_primary_multi_row = project_dir + "/gene-annos_primary_multi-row.bed"
# ** Specificity
# *** Probes sequences

specificity_analysis_dir = project_dir + '/probe-specificity'
probe_sequences_dir = specificity_analysis_dir + '/probe-sequences'
os.makedirs(probe_sequences_dir, exist_ok = True)
prob_seqs_by_probe_set_name = probe_sequences_dir + '/probe-seqs_{probe_set_name}.fq'
probe_seqs_v1_fq = prob_seqs_by_probe_set_name.format(probe_set_name = 'probes_v1')

# *** Alignment

probe_alignments_dir = specificity_analysis_dir + '/probe-alignments'
os.makedirs(probe_alignments_dir, exist_ok = True)
probe_alignment_bam = probe_alignments_dir + '/probes.bam'
# ** Alignments

biscuit_bam_by_ref_genome_probe_set = project_dir + '/alignments/{ref_genome}/{probe_set_name}/{ref_genome}_{probe_set_name}.bam'
