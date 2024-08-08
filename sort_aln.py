
import subprocess
import sys
import os
import yaml

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
SAMTOOLS_MEM = config.get("samtools_mem", 1)

# Output and log file paths
input_bam = f"temp/{SM}.{W}.{F}.bam"
output_bam = f"results/{SM}.{W}.{F}.sorted.bam"
log_file = f"logs/sort_aln.{SM}.{W}.{F}.log"

# Create the command
command = f"samtools sort -m {SAMTOOLS_MEM}G -o {output_bam} {input_bam}"

# Execute the command
try:
    os.makedirs(os.path.dirname(output_bam), exist_ok=True)
    os.makedirs(os.path.dirname(log_file), exist_ok=True)
    with open(log_file, "w") as log:
        subprocess.run(command, shell=True, check=True, stderr=log)
    print(f"Sorted BAM file created: {output_bam}")
except subprocess.CalledProcessError as e:
    print(f"Error running samtools sort: {e}", file=sys.stderr)
    sys.exit(1)
