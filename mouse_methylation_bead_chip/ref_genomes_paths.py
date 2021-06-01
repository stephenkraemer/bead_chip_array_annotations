# * General considerations

# - we use UCSC style identifiers, like the analysis sets for mm10 and hg38
# - do not use Alt contigs
#   - mapping quality estimates under consideration of Alt contigs for bisulfite data is not established
#   - and it is also not necessary / high priority
#     - for mouse
#       - variation between strands is better represented for our purposes by strand-specific analysis
#       - note that GRCm39 does not contain ALT contigs anymore
#     - for humans
#       - polymorphisms could be taken into account using a graph-based aligner, such as hisat-3n (we are still evaluating whether this is worth the effort)
# - no hard masking of small repeats
# - hard masking of PAR and centromere repeats is partially out of our control, this is indicated for each reference genome
# - use decoy sequences where appropriate since 'contaminations' could also affect probe specificity
#   - only available for human sequences at the moment


# * Sources
# ** Mouse

# *** C57BL/6J (GRCm*)

# **** GRCm38

# - analysis set without alt regions
# - no decoys

mm10_fa_gz_url = 'https://ftp.ncbi.nlm.nih.gov/genomes/archive/old_genbank/Eukaryotes/vertebrates_mammals/Mus_musculus/GRCm38.p3/seqs_for_alignment_pipelines/GCA_000001635.5_GRCm38.p3_no_alt_analysis_set.fna.gz'

# **** GRCm39

# - NCBI should provide analysis set ( https://www.ncbi.nlm.nih.gov/genome/doc/ftpfaq/#fastaseqid), but the https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/635/GCA_000001635.9_GRCm39/seqs_for_alignment_pipelines.ucsc_ids file does not exist
# - use ucsc soft masked, notes
#   - PAR regions are hard masked (checked by hand)
# - no decoy sequences

mm39_fa_gz_url = 'https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/mm39.fa.gz'

# *** Other common strains

# - these are assemblies from the MGP, retrieved via UCSC
# - the extent of hard masking is not documented
# - we attempt to cover some common strains, could be extended upon request
# - it seems that only 2-bit genomes are available via UCSC

# **** BALB_cJ

# - chrY missing

balb_fa_2bit_url = 'http://hgdownload.soe.ucsc.edu/hubs/mouseStrains/GCA_001632525.1_BALB_cJ_v1/GCA_001632525.1_BALB_cJ_v1.2bit'


# **** A_J

aj_fa_2bit_url = 'http://hgdownload.soe.ucsc.edu/hubs/mouseStrains/GCA_001624215.1_A_J_v1/GCA_001624215.1_A_J_v1.2bit'

# **** 129S1_SvImJ

s1_fa_2bit_url = 'http://hgdownload.soe.ucsc.edu/hubs/mouseStrains/GCA_001624185.1_129S1_SvImJ_v1/GCA_001624185.1_129S1_SvImJ_v1.2bit'

# **** C57BL_6NJ

bl6nj_fa_2bit_url = 'http://hgdownload.soe.ucsc.edu/hubs/mouseStrains/GCA_001632555.1_C57BL_6NJ_v1/GCA_001632555.1_C57BL_6NJ_v1.2bit'

# **** C3H_HeJ

c3h_fa_2bit_url = 'http://hgdownload.soe.ucsc.edu/hubs/mouseStrains/GCA_001632575.1_C3H_HeJ_v1/GCA_001632575.1_C3H_HeJ_v1.2bit'

# ** Human

# *** hg19

# - GRCh37.p13 analysis set as provided by NCBI
#  - masked PAR regions
#  - no alt sequences
#  - EBV decoy
#  - current chrM sequence

hg19_fa_gz_url = 'https://ftp.ncbi.nlm.nih.gov/genomes/archive/old_genbank/Eukaryotes/vertebrates_mammals/Homo_sapiens/GRCh37.p13/seqs_for_alignment_pipelines/GCA_000001405.14_GRCh37.p13_no_alt_analysis_set.fna.gz'


# *** hg38

# - analysis set
#   - masked PAR and centromeric and genomic repeat arrays
#   - no alt contigs
#   - with HS38d1 decoy, which contains viral sequences, eg EBV which may well be present in a human NGS sample, as well as human genomic sequences which are not included in the assembly

hg38_fa_gz_url = 'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz'

# * Setup

import mouse_methylation_bead_chip.beadchip_probe_annotation_paths as paths
import os

ref_genomes_dir = paths.project_dir + "/ref-genomes"
os.makedirs(ref_genomes_dir, exist_ok = True)

# # **** mm10

# mm10_genomes_dir = ref_genomes_dir + '/mm10'
# os.makedirs(mm10_genomes_dir, exist_ok = True)

# # full analysis set fasta download path
# mm10_full_fa_gz_url = 'https://ftp.ncbi.nlm.nih.gov/genomes/archive/old_genbank/Eukaryotes/vertebrates_mammals/Mus_musculus/GRCm38.p3/seqs_for_alignment_pipelines/GCA_000001635.5_GRCm38.p3_full_analysis_set.fna.gz'
# mm10_full_fa_gz = mm10_genomes_dir + '/GCA_000001635.5_GRCm38.p3_full_analysis_set.fna.gz'

# mm10_full_fa_re_gz = mm10_genomes_dir + '/GCA_000001635.5_GRCm38.p3_full_analysis_set_re-gz.fna.gz'

# mm10_full_fa_re_gz_chr19_only = mm10_genomes_dir + '/GCA_000001635.5_GRCm38.p3_full_analysis_set_re-gz_chrom-19.fna.gz'

# * Ref genome urls

ref_genome_urls = {
    'mm10': mm10_fa_gz_url,
    'mm39': mm39_fa_gz_url,
    'hg19': hg19_fa_gz_url,
    'hg38': hg38_fa_gz_url,
    'balb': balb_fa_2bit_url,
    'aj': aj_fa_2bit_url,
    's1': s1_fa_2bit_url,
    'c3h': c3h_fa_2bit_url,
    'bl6nj': bl6nj_fa_2bit_url,
}

# * File paths

# - for all files, indices etc. will be placed next to the file
# - directories are autocreated by snakemake during download

fa_gz_pattern = ref_genomes_dir + '/{ref_genome}/{ref_genome}.fa.gz'
twobit_pattern = ref_genomes_dir + '/{ref_genome}/{ref_genome}.2bit'
# ** Mouse reference genoms

mm10_fa_gz = ref_genomes_dir + '/mm10/mm10.fa.gz'
mm39_fa_gz = ref_genomes_dir + '/mm39/mm39.fa.gz'
balb_fa_2bit = ref_genomes_dir + '/balb/balb.2bit'
c3h_fa_2bit = ref_genomes_dir + '/c3h/c3h.2bit'
bl6nj_fa_2bit =  ref_genomes_dir + '/bl6nj/bl6nj.2bit'
s1_fa_2bit =  ref_genomes_dir + '/s1/s1.2bit'
aj_fa_2bit =  ref_genomes_dir + '/aj/aj.2bit'
balb_fa_gz =  ref_genomes_dir + '/balb/balb.fa.gz'
c3h_fa_gz =  ref_genomes_dir + '/c3h/c3h.fa.gz'
bl6nj_fa_gz =  ref_genomes_dir + '/bl6nj/bl6nj.fa.gz'
s1_fa_gz =  ref_genomes_dir + '/s1/s1.fa.gz'
aj_fa_gz =  ref_genomes_dir + '/aj/aj.fa.gz'

# ** Human reference genoms

hg19_fa_gz = ref_genomes_dir + '/hg19/hg19.fa.gz'
hg38_fa_gz = ref_genomes_dir + '/hg38/hg38.fa.gz'
