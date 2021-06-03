import mouse_methylation_bead_chip.ref_genomes_paths as ref_genome_paths
import mouse_methylation_bead_chip.beadchip_probe_annotation_paths as paths
from snakemake.io import temp, touch, expand

log_dir = paths.project_dir + '/snakemake-logs'

# * All

rule all:
    input:
        expand(
            ref_genome_paths.biscuit_index_done_by_ref_genome,
            ref_genome = ref_genome_paths.ref_genomes,
            )

        # ref_genome_paths.balb_fa_gz,
        # ref_genome_paths.mm39_fa_gz,

        # ref_genome_paths.mm10_fa_gz,
        # ref_genome_paths.c3h_fa_gz,
        # ref_genome_paths.bl6nj_fa_gz,
        # ref_genome_paths.s1_fa_gz,
        # ref_genome_paths.aj_fa_gz,
        # ref_genome_paths.hg19_fa_gz,
        # ref_genome_paths.hg38_fa_gz,

        # ref_genome_paths.balb_fa_2bit,
        # ref_genome_paths.c3h_fa_2bit,
        # ref_genome_paths.bl6nj_fa_2bit,
        # ref_genome_paths.s1_fa_2bit,
        # ref_genome_paths.aj_fa_2bit,

# * Download

rule download_and_bgzip_index_ref_genome_fa_gz:
    output:
        ref_fa_gz = ref_genome_paths.fa_gz_by_ref_genome,
        # ref_fai = ref_genome_paths.ref_genomes_dir + '/{ref_genome}/{ref_genome}.fai',
        # used either for
        # 1.save fa_gz to temp during download so that we can then bgzip it to the final destination
        # 2. 2bit file has to be saved to disk cannot be piped, then it still needs to be bgzipped and indexed
        ref_fa_gz_temp = temp(ref_genome_paths.fa_gz_by_ref_genome + '.tmp'),
        two_bit = temp(ref_genome_paths.twobit_by_ref_genome),
    params:
        url = lambda wildcards: ref_genome_paths.ref_genome_urls[wildcards.ref_genome],
    resources:
        avg_mem = 500,
        mem_mb = 500,
        walltime_min = lambda wildcards, attempt: 10 * attempt,
    threads: 1
    log:
        out = log_dir + '/download_{ref_genome}.out',
        err = log_dir + '/download_{ref_genome}.err',
    shell:
        """
        exec >{log.out}
        exec 2>{log.err}

        echo "starting download"
        if [[ {params.url} == *.2bit ]]; then
            wget -O {output.two_bit} {params.url}
            # ~ 1 min
            twoBitToFa {output.two_bit} {output.ref_fa_gz_temp}
            echo "bgzip"
            # ~ 1 min
            bgzip --stdout {output.ref_fa_gz_temp} > {output.ref_fa_gz}
        else
            # is this necessary?
            touch {output.two_bit}
            wget -O {output.ref_fa_gz_temp} {params.url}
            echo "bgzip"
            zcat {output.ref_fa_gz_temp} | bgzip > {output.ref_fa_gz}
        fi

        echo "index"
        samtools faidx {output.ref_fa_gz}
        """

# * Biscuit

rule biscuit_index:
    input:
        ref_fa_gz = ref_genome_paths.fa_gz_by_ref_genome,
    output:
        touch(ref_genome_paths.biscuit_index_done_by_ref_genome),
    resources:
        avg_mem = 6000,
        mem_mb = 6000,
        walltime_min = lambda wildcards, attempt: 7 * 60 * attempt,
    threads: 1
    log:
        out = log_dir + '/biscuit-index_{ref_genome}.out',
        err = log_dir + '/biscuit-index_{ref_genome}.err',
    shell:
        """
        exec >{log.out}
        exec 2>{log.err}
        biscuit index {input.ref_fa_gz}
        """

rule biscuit_align:
    input:
        biscuit_index = ref_genome_paths.biscuit_index_done_by_ref_genome,
        fa_gz = ref_genome_paths.fa_gz_by_ref_genome,
        probe_sequences_fq = paths.prob_seqs_by_probe_set_name,
    output:
        bam = paths.biscuit_bam_by_ref_genome_probe_set,
    threads: 16
    resources:
        avg_mem = 6000,
        mem_mb = 6000,
        walltime_min = lambda wildcards, attempt: 20 * attempt,
    log:
        out = log_dir + '/align_{ref_genome}_{probe_set_name}.out',
        err = log_dir + '/align_{ref_genome}_{probe_set_name}.err',
    shell:
        """
        exec >{log.out}
        exec 2>{log.err}
        biscuit align \
                -t {threads} \
                {input.fa_gz} \
                {input.probe_sequences_fq} \
            > {output.bam}
        """
