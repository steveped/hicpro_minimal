import pandas as pd
import os
import re

configfile: "config/config.yml"

# Samples
df = pd.read_table(config["samples"])
samples = df['sample'] + "/" + df['file']
suffix = config['suffix']

#################################
## Variables for the reference ##
#################################
build = config['ref']['build']
ref_root = os.path.join(config['ref']['root'], "gencode-release-" + str(config['ref']['gencode']),
                        build, "dna")
# Key output files
assembly = config['ref']['assembly'] + ".genome"
ref_fa = build + "." + assembly + ".fa"
ref_fagz = ref_fa + ".gz"
chr_sizes = os.path.join(os.getcwd(), "config", build + ".chr_sizes.tsv")
rs_frags = os.path.join(os.getcwd(), "config", build + "_" + config['hicpro']['enzyme'] + "_fragment.bed")

#####################
## HiC-Pro outputs ##
#####################
bins = re.split(r" ", config['hicpro']['bin_size'])
hicpro_config = "config/hicpro-config.txt"
digest_script = "scripts/digest_genome.py"
MAPPING = expand(["data/hic/bowtie_results/bwt2/{sample}_{reads}_" + build + "." + assembly + ".bwt2merged.bam"],
              reads = ['R1', 'R2'], sample = samples)
PROC_BAM = expand(["data/hic/bowtie_results/bwt2/{sample}_" + build + "." + assembly + ".bwt2pairs.bam"],
                  sample = samples)
PROC_PAIRS = expand(["data/hic/hic_results/data/{sample}_" + build + "." + assembly + ".bwt2pairs.validPairs"],
                    sample = samples)
HIC_QC = ['data/hic/hic_results/pic']

#HIC_PAIRS = expand(["data/hic/hic_results/data/{sample}_allValidPairs"],
#                   sample = samples)
#HIC_MAT = expand(["data/hic/hic_results/matrix/{sample}/raw/{bin}/{sample}_{bin}.matrix"],
#                 bin = bins, sample = samples)
#HIC_BED = expand(["data/hic/hic_results/matrix/{sample}/raw/{bin}/{sample}_{bin}_abs.bed"],
#                 bin = bins, sample = samples)

## Define all the required outputs as a single object
REFS = [chr_sizes, rs_frags]
ALL_OUTPUTS = []
ALL_OUTPUTS.extend(REFS)
ALL_OUTPUTS.extend([hicpro_config, digest_script])
ALL_OUTPUTS.extend(MAPPING)
ALL_OUTPUTS.extend(PROC_BAM)
ALL_OUTPUTS.extend(PROC_PAIRS)
ALL_OUTPUTS.extend(HIC_QC)

rule all:
    input:
        ALL_OUTPUTS

include: "rules/hicpro.smk"


