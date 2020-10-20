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
ref_root = os.path.join(config['ref']['root'], "gencode-release-" + str(config['ref']['gencode']),
                        config['ref']['build'], "dna")
# Key output files
assembly = config['ref']['assembly'] + ".genome"
ref_fa = config['ref']['build'] + "." + assembly + ".fa"
ref_fagz = ref_fa + ".gz"
chr_sizes = os.path.join(os.getcwd(), "config", config['ref']['build'] + ".chr_sizes.tsv")
rs_frags = os.path.join(os.getcwd(), "config", config['ref']['build'] + "_" + config['hicpro']['enzyme'] + "_fragment.bed")

#####################
## HiC-Pro outputs ##
#####################
bins = re.split(r" ", config['hicpro']['bin_size'])
hicpro_config = "config/hicpro-config.txt"
digest_script = "scripts/digest_genome.py"
HIC_PAIRS = expand(["data/hic/hic_results/data/{sample}_allValidPairs"],
                   sample = samples)
HIC_MAT = expand(["data/hic/hic_results/matrix/{sample}/raw/{bin}/{sample}_{bin}.matrix"],
                 bin = bins, sample = samples)
HIC_BED = expand(["data/hic/hic_results/matrix/{sample}/raw/{bin}/{sample}_{bin}_abs.bed"],
                 bin = bins, sample = samples)

## Define all the required outputs as a single object
REFS = [chr_sizes, rs_frags]
FAGZ = [os.path.join(ref_root, ref_fagz)]
BOWTIEIDX = expand([ref_root + "/bt2/{prefix}.{sub}.bt2"],
               prefix = config['ref']['build'] + "." + assembly,
               sub = ['1', '2', '3', '4', 'rev.1', 'rev.2'] )
ALL_OUTPUTS = []
ALL_OUTPUTS.extend(REFS)
ALL_OUTPUTS.extend([BOWTIEIDX])
ALL_OUTPUTS.extend(FAGZ)
ALL_OUTPUTS.extend([hicpro_config, digest_script])
ALL_OUTPUTS.extend(HIC_PAIRS)
ALL_OUTPUTS.extend(HIC_MAT)
ALL_OUTPUTS.extend(HIC_BED)

rule all:
    input:
        ALL_OUTPUTS

include: "rules/reference.smk"
include: "rules/hicpro.smk"

