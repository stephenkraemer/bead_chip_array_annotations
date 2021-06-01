import mouse_methylation_bead_chip.ref_genomes_paths as ref_genome_paths
import mouse_methylation_bead_chip.beadchip_probe_annotation_paths as paths
from snakemake.io import temp, touch

log_dir = paths.project_dir + '/snakemake-logs'

# * All

rule all:
    input:
        ref_genome_paths.mm10_fa_gz
        ref_genome_paths.mm39_fa_gz
        ref_genome_paths.balb_fa_2bit
        ref_genome_paths.c3h_fa_2bit
        ref_genome_paths.bl6nj_fa_2bit
        ref_genome_paths.s1_fa_2bit
        ref_genome_paths.aj_fa_2bit
        # ref_genome_paths.balb_fa_gz
        # ref_genome_paths.c3h_fa_gz
        # ref_genome_paths.bl6nj_fa_gz
        # ref_genome_paths.s1_fa_gz
        # ref_genome_paths.aj_fa_gz
        ref_genome_paths.hg19_fa_gz
        ref_genome_paths.hg38_fa_gz

# * Download

rule download_and_bgzip_index_ref_genome_fa_gz:
    output:
        ref_fa_gz = ref_genome_paths.fa_gz_pattern,
        # ref_fai = ref_genome_paths.ref_genomes_dir + '/{ref_genome}/{ref_genome}.fai',
        # used either for
        # 1.save fa_gz to temp during download so that we can then bgzip it to the final destination
        # 2. 2bit file has to be saved to disk cannot be piped, then it still needs to be bgzipped and indexed
        ref_fa_gz_temp = temp(ref_genome_paths.fa_gz_pattern + '.tmp'),
        two_bit = ref_genome_paths.twobit_pattern,
    params:
        url = lambda wildcards: ref_genome_paths.ref_genome_urls[wildcards.ref_genome],
        avg_mem = 500,
        mem_mb = 500,
        walltime_min = lambda attempt: 10 * attempt,
    threads: 1
    log:
        out = log_dir + '/download_{ref_genome}.out',
        err = log_dir + '/download_{ref_genome}.err',
    shell:
        """
        exec >{log.out}
        exec 2>{log.err}
        # note that the temp file could also be a 2bit file, file name is not ideal
        if [[ {url} == *.2bit ]]; then
            wget -O {output.two_bit} {url}
            twoBitToFa {output.two_bit} {output.ref_fa_gz_temp}
        else
            wget -O {output.ref_fa_gz_temp} {url}
        fi
        zcat {output.ref_fa_gz_temp} | bgzip > {output.ref_fa_gz}
        samtools faidx {output.ref_fa_gz}
        """

# * Biscuit

biscuit_index_done_pattern = ref_genome_paths.ref_genomes_dir + '/{ref_genome}/{reg_genome}.fa.gz' + '_biscuit-index.done'
rule biscuit_index:
    input:
        ref_fa_gz = ref_genome_paths.ref_genomes_dir + '/{ref_genome}/{reg_genome}.fa.gz',
    output:
        touch(biscuit_index_done_pattern),
    params:
        avg_mem = 6000,
        mem_mb = 6000,
        walltime_min = lambda attempt: 7 * 60 * attempt,
    threads: 1
    log:
        out = log_dir + '/download_{ref_genome}.out',
        err = log_dir + '/download_{ref_genome}.err',
    shell:
        """
        exec >{log.out}
        exec 2>{log.err}
        biscuit index {input.ref_fa_gz}
        """

rule biscuit_align:
    input:
        biscuit_index = biscuit_index_done_pattern,
        fa_gz = fa_gz_pattern,
        probe_sequences_fq = paths.probe_sequences_fq,
    output:
        bam = paths.biscuit_bam_pattern,
    threads: 16
    params:
        avg_mem = 6000,
        mem_mb = 6000,
        walltime_min = lambda attempt: 20 * attempt,
    log:
        out = log_dir + '/download_{ref_genome}.out',
        err = log_dir + '/download_{ref_genome}.err',
    shell:
        """
        exec >{log.out}
        exec 2>{log.err}
        biscuit align \
                -t {threads} \
                {input.fa_gz} \
                {input.probe_sequences_fq} \
            > {bam}
        """
