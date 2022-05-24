from pathlib import Path

# # Paths


# # Imports

import os
# # Project paths

project_dir = (
    "/omics/groups/OE0029/internal/kraemers/projects/epic-arrays-for-hematopoiesis"
)

# # Reference genome

mm10_fa_url = 'http://hgdownload.cse.ucsc.edu/goldenpath/mm10/bigZips/mm10.fa.gz'
mm10_fa = project_dir + '/ucsc_mm10.fa.gz'
mm10_fa_bgz = project_dir + '/ucsc_mm10.fa.bgz'

ncbi_mm10_fa_url = 'https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Mus_musculus/all_assembly_versions/GCF_000001635.26_GRCm38.p6/GCF_000001635.26_GRCm38.p6_genomic.fna.gz'
ncbi_mm10_fa = project_dir + '/ncbi_mm10.fa.gz'
ncbi_mm10_fa_bgz = project_dir + '/ncbi_mm10.fa.bgz'

chrom_alias_url = 'https://hgdownload.cse.ucsc.edu/goldenPath/mm10/database/chromAlias.txt.gz'
chrom_alias_txt_gz = project_dir + '/chromAlias.txt.gz'

# # Gencode

gencode_download_url = "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.annotation.gtf.gz"

# unmodified gencode gtf (note: with chr prefix)

gencode_gtf = project_dir + "/gencode.vM25.annotation.gtf.gz"

# gencode filtered for principal transcripts of protein coding genes, note that that chromosome prefix ('chr') is removed in this file

gencode_coding_canonical_gtf = (
    project_dir + "/gencode.vM25.annotation_coding_canonical.gtf.gz"
)
# # CpG islands
cpg_islands_ucsc_unmasked_url = 'http://hgdownload.cse.ucsc.edu/goldenpath/mm10/database/cpgIslandExt.txt.gz'
cpg_islands_ucsc_unmasked_txt_gz = project_dir + "/cpg-islands_ucsc_unmasked.txt.gz"

# # Probes

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

# # Gene annotation results

# ## gtfanno results

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

# ## Final tables

# v1 annos do not annotate CHR in ['0', nan] probes

gene_annos_primary_one_row = project_dir + "/gene-annos_primary_one-row.bed"
gene_annos_primary_multi_row = project_dir + "/gene-annos_primary_multi-row.bed"

# v2 probes provide annos for all illumina probes in manifest order
# annos are NA for CHR in ['0', nan] probes

probe_annos_one_row_bed_csv = project_dir + "/gene-annos_primary_one-row_v2.bed"
probe_annos_one_row_bed_parquet = project_dir + "/gene-annos_primary_one-row_v2.parquet"

# # Specificity

# ## Probes sequences

specificity_analysis_dir = project_dir + '/probe-specificity'
probe_sequences_dir = specificity_analysis_dir + '/probe-sequences'
os.makedirs(probe_sequences_dir, exist_ok = True)
prob_seqs_by_probe_set_name = probe_sequences_dir + '/probe-seqs_{probe_set_name}.fq'
probe_seqs_v1_fq = prob_seqs_by_probe_set_name.format(probe_set_name = 'probes_v1')

# ## Alignment

probe_alignments_dir = specificity_analysis_dir + '/probe-alignments'
os.makedirs(probe_alignments_dir, exist_ok = True)
probe_alignment_bam = probe_alignments_dir + '/probes.bam'
# # Alignments

biscuit_bam_by_ref_genome_probe_set = project_dir + '/alignments/{ref_genome}/{probe_set_name}/{ref_genome}_{probe_set_name}.bam'

# # Auto-create dirs

for value in globals().copy().values():
    if (
        isinstance(value, str)
        and value.startswith(project_dir)
        and not "{" in value
    ):
        # print(value)
        Path(value).parent.mkdir(exist_ok=True, parents=True)
