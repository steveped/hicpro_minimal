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

rule run_hicpro_mapping:
    input:
        config = hicpro_config,
        files = expand(["data/test_data/{sample}{reads}{suffix}"],
                       sample = samples, suffix = suffix,
                       reads = [config['hicpro']['pair1_ext'], config['hicpro']['pair2_ext']])
    output:
        bam = expand(["data/hic/bowtie_results/bwt2/{sample}_" + build + "." + assembly + ".bwt2pairs.bam"],
                     sample = samples)
    param:
        indir = "data/test_data",
        outdir = "data/hic"
    log: "logs/hicpro/run_hicpro_mapping.log"
    threads: config['hicpro']['ncpu']
    shell:
        """
        ######################################
        ## Specific to phoenix for now only ##
        ######################################
        ## Load modules
        module load HiC-Pro/2.9.0-foss-2016b

        ##Run HiC-pro responding to yes to any interactive requests
        yes | HiC-Pro \
          -s mapping \
          -c {input.config} \
          -i {params.indir} \
          -o {params.outdir} &> {log}
        """



