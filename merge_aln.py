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
FASTA = os.path.abspath(config["fasta"])
SM = config["sample"]
W = config["window"]
F = config["mm_f"]
SLIDE = config.get("slide", 0)

# Calculate N based on the configuration and the sequence lengths
def calculate_n(fasta, window):
    fai = f"{fasta}.fai"
    if not os.path.exists(fai):
        sys.exit(f"Input fasta must be indexed, try:\nsamtools faidx {fasta}")

    lengths = [math.ceil(int(line.split()[1]) / window) for line in open(fai).readlines()]
    return int(sum(lengths))

N = calculate_n(FASTA, W)

# Output and log file paths
input_list = f"temp/{SM}.{W}.{F}.list"
output_bam = f"temp/{SM}.{W}.{F}.bam"
log_file = f"logs/merge_aln.{SM}.{W}.{F}.log"

# Verify that the list file exists
if not os.path.exists(input_list):
    sys.exit(f"List file not found: {input_list}")

# Verify that all BAM files referenced in the list file exist
missing_files = []
with open(input_list, 'r') as file:
    for line in file:
        bam_file = line.strip()
        if not os.path.exists(bam_file):
            missing_files.append(bam_file)

if missing_files:
    sys.exit(f"The following BAM files are missing:\n" + "\n".join(missing_files))

# Create the command
command = f"samtools merge -b {input_list} {output_bam}"

# Execute the command
try:
    os.makedirs(os.path.dirname(output_bam), exist_ok=True)
    os.makedirs(os.path.dirname(log_file), exist_ok=True)
    with open(log_file, "w") as log:
        result = subprocess.run(command, shell=True, check=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
        log.write(result.stderr.decode())
        log.write(result.stdout.decode())
        print(f"Merged BAM file created: {output_bam}")
except subprocess.CalledProcessError as e:
    print(f"Error running samtools merge: {e.stderr.decode()}", file=sys.stderr)
    sys.exit(1)
except Exception as e:
    print(f"Unexpected error: {str(e)}", file=sys.stderr)
    sys.exit(1)
