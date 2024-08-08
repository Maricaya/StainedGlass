
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
FASTA = os.path.abspath(config["fasta"])
SM = config["sample"]
W = config["window"]

# Function to check if FAI file exists
def check_fai(fasta):
    fai = f"{fasta}.fai"
    if not os.path.exists(fai):
        sys.exit(f"Input fasta must be indexed, try:\nsamtools faidx {fasta}")
    return fai

# Check if FAI file exists
FAI = check_fai(FASTA)

# Output and log file paths
output_fasta = f"results/{SM}.{W}.fasta"
log_file = f"logs/window_fa.{SM}.{W}.log"
input_bed = f"temp/{SM}.{W}.bed"

# Create the command
command = f"bedtools getfasta -fi {FASTA} -bed {input_bed} > {output_fasta}"

# Execute the command
try:
    os.makedirs(os.path.dirname(output_fasta), exist_ok=True)
    os.makedirs(os.path.dirname(log_file), exist_ok=True)
    with open(log_file, "w") as log:
        subprocess.run(command, shell=True, check=True, stderr=log)
    print(f"FASTA file created: {output_fasta}")
except subprocess.CalledProcessError as e:
    print(f"Error running bedtools getfasta: {e}", file=sys.stderr)
    sys.exit(1)
