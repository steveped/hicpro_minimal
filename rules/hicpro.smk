rule get_chrom_sizes:
    output: chr_sizes
    params:
        ftp = "ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/vertebrate_mammalian/Homo_sapiens/all_assembly_versions",
        genbank = config['ref']['genbank'],
        build = config['ref']['build']
    threads: 1
    shell:
        """

        # Download the assembly report
        TEMPDIR=$(mktemp -d -t chrXXXXXXXXXX)
        REPORT="{params.genbank}_{params.build}_assembly_report.txt"
        URL="{params.ftp}/{params.genbank}_{params.build}/{params.genbank}_{params.build}_assembly_report.txt"
        wget -O "$TEMPDIR/$REPORT" $URL

        # Extract the chrom_sizes
        egrep 'assembled-molecule' "$TEMPDIR/$REPORT" | \
          awk '{{print $11"\t"$10}}' > {output}

        rm -rf $TEMPDIR

        """

rule find_rs_fragments:
    input: '/hpcfs/users/a1018048/refs/gencode-release-33/GRCh37/dna/GRCh37.primary_assembly.genome.fa'
    output:
        script = "scripts/digest_genome.py",
        rs = rs_frags
    params:
        enzyme = config['hicpro']['enzyme']
    threads: 1
    conda: "../envs/python2.7.yml"
    shell:
        """
        # Get v2.9.0 from the HiC-Pro repo
        wget \
          -O {output.script} \
          "https://raw.githubusercontent.com/nservant/HiC-Pro/2d15209fbb75ce3278d68801bd98be4b2416e5b5/bin/utils/digest_genome.py"

        # Run the python script
        python scripts/digest_genome.py \
          -r {params.enzyme} \
          -o {output.rs} \
          {input}
        """

rule make_hicpro_config:
    input:
        idx = '/hpcfs/users/a1018048/refs/gencode-release-33/GRCh37/dna/bt2/',
        rs = rs_frags,
        chr_sizes = chr_sizes
    output:
        hicpro_config
    conda: "../envs/stringr.yml"
    threads: 1
    shell:
        """
        Rscript --vanilla \
          scripts/write_hicpro_config.R \
          {input.idx} \
          {input.chr_sizes} \
          {input.rs} \
          {output}
        """

rule hicpro_mapping:
    input:
        config = hicpro_config,
        files = expand(["data/test_data/{sample}{reads}{suffix}"],
                       sample = samples, suffix = suffix,
                       reads = [config['hicpro']['pair1_ext'], config['hicpro']['pair2_ext']])
    output:
        bam = expand(["data/hic/bowtie_results/bwt2/{sample}_{reads}_" + build + "." + assembly + ".bwt2merged.bam"],
              reads = ['R1', 'R2'], sample = samples)
    params:
        indir = "data/test_data",
        outdir = "data/hic"
    log: "logs/hicpro/hicpro_mapping.log"
    threads: config['hicpro']['ncpu']
    shell:
        """
        ######################################
        ## Specific to phoenix for now only ##
        ######################################
        ## Load modules
        module load HiC-Pro/2.9.0-foss-2016b

        ## Remove any existing data as leaving this here causes HicPro to
        ## make an interactive request. Piping `yes` into HicPro may be the
        ## source of some current problems
        if [[ -d {params.outdir} ]]; then
          rm -rf {params.outdir}
        fi

        ##Run HiC-pro responding to yes to any interactive requests
        HiC-Pro \
          -s mapping \
          -c {input.config} \
          -i {params.indir} \
          -o {params.outdir} &> {log}
        """

rule hicpro_proc:
    input:
        config = hicpro_config,
        files = rules.hicpro_mapping.output.bam
    output:
        bam = expand(["data/hic/bowtie_results/bwt2/{sample}_" + build + "." + assembly + ".bwt2pairs.bam"],
                     sample = samples),
        pairs = expand(["data/hic/hic_results/data/{sample}_" + build + "." + assembly + ".bwt2pairs.validPairs"],
                     sample = samples)
    params:
        indir = "data/hic/bowtie_results/bwt2",
        outdir = "data/hic"
    log: "logs/hicpro/hicpro_proc.log"
    threads: config['hicpro']['ncpu']
    shell:
        """
        ######################################
        ## Specific to phoenix for now only ##
        ######################################
        ## Load modules
        module load HiC-Pro/2.9.0-foss-2016b

        ##Run HiC-pro responding to yes to any interactive requests
        HiC-Pro \
          -s proc_hic \
          -c {input.config} \
          -i {params.indir} \
          -o {params.outdir} &> {log}
        """

rule hicpro_qc:
    input:
        config = hicpro_config,
        files = rules.hicpro_mapping.output.bam
    output:
        pic = directory("data/hic/hic_results/pic")
    params:
        indir = "data/hic/bowtie_results/bwt2",
        outdir = "data/hic"
    log: "logs/hicpro/hicpro_qc.log"
    threads: config['hicpro']['ncpu']
    shell:
        """
        ######################################
        ## Specific to phoenix for now only ##
        ######################################
        ## Load modules
        module load HiC-Pro/2.9.0-foss-2016b

        ##Run HiC-pro responding to yes to any interactive requests
        HiC-Pro \
          -s quality_checks \
          -c {input.config} \
          -i {params.indir} \
          -o {params.outdir} &> {log}
        """

rule hicpro_merge:
    input:
        config = hicpro_config,
        files = rules.hicpro_proc.output.pairs
    output:
        pairs = expand(["data/hic/hic_results/data/{rep}/{rep}_allValidPairs"],
                       rep = df['sample'])
    params:
        indir = "data/hic/hic_results/data",
        outdir = "data/hic"
    log: "logs/hicpro/hicpro_merge.log"
    threads: config['hicpro']['ncpu']
    shell:
        """
        ######################################
        ## Specific to phoenix for now only ##
        ######################################
        ## Load modules
        module load HiC-Pro/2.9.0-foss-2016b

        ##Run HiC-pro responding to yes to any interactive requests
        HiC-Pro \
          -s merge_persample \
          -c {input.config} \
          -i {params.indir} \
          -o {params.outdir} &> {log}
        """

rule build_contact_maps:
    input:
        config = hicpro_config,
        pairs = rules.hicpro_merge.otput.pairs
    output:
        bed = expand(["data/hic/hic_results/matrix/{rep}/{bin}/{rep}_{bin}_{type}.bed"],
                     rep = df['sample'], bin = bins, type = ['abs', 'ord']),
        mat = expand(["data/hic/hic_results/matrix/{rep}/{bin}/{rep}_{bin}.matrix"],
                     rep = df['sample'], bin = bins)
    params:
        indir = "data/hic/hic_results/data",
        outdir = "data/hic"
    log: "logs/hicpro/build_contact_maps.log"
    threads: config['hicpro']['ncpu']
    shell:
        """
        ######################################
        ## Specific to phoenix for now only ##
        ######################################
        ## Load modules
        module load HiC-Pro/2.9.0-foss-2016b

        ##Run HiC-pro responding to yes to any interactive requests
        HiC-Pro \
          -s build_contact_maps \
          -c {input.config} \
          -i {params.indir} \
          -o {params.outdir} &> {log}
        """
