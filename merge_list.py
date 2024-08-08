import subprocess
import sys
import os
import yaml
import math

# Load configuration from config/config.yaml
config_path = 'config/config.yaml'
if not os.path.exists(config_path):
    sys.exit(f"Configuration file not found: {config_path}")

with open(config_path, 'r') as config_file:
    config = yaml.safe_load(config_file)

# Set parameters from the configuration file
SM = config["sample"]
W = config["window"]
F = config["mm_f"]
SLIDE = config.get("slide", 0)
N = config.get("nbatch", 1)

# Check if the input fasta file is indexed
FASTA = os.path.abspath(config["fasta"])
FAI = f"{FASTA}.fai"
if not os.path.exists(FAI):
    sys.exit(f"Input fasta must be indexed, try:\nsamtools faidx {FASTA}")

# Determine the IDs dynamically based on sequence splitting
lengths = [math.ceil(int(line.split()[1]) / W) for line in open(FAI).readlines()]
names = [line.split()[0] for line in open(FAI).readlines()]
REF_IDS = [f"ref_{i}" for i in range(len(names))] if SLIDE > 0 else ["ref_0"]
N_WINDOWS = int(sum(lengths))
N = min(N, N_WINDOWS)
IDS = list(range(N))
sys.stderr.write(f"[INFO] The sequence will be split into {N} batches.\n")

# Output and log file paths
output_list = f"temp/{SM}.{W}.{F}.list"
log_file = f"logs/merge_list.{SM}.{W}.{F}.log"

# Create the list of BAM files
bam_files = [f"temp/{SM}.{W}.{F}.{ID}.{REF_ID}.bam" for ID in IDS for REF_ID in REF_IDS]

# Write the list of BAM files to the output list file
try:
    os.makedirs(os.path.dirname(output_list), exist_ok=True)
    os.makedirs(os.path.dirname(log_file), exist_ok=True)
    with open(output_list, "w") as list_file:
        list_file.write("\n".join(bam_files) + "\n")
    print(f"List file created: {output_list}")
except Exception as e:
    print(f"Error creating list file: {e}", file=sys.stderr)
    sys.exit(1)
