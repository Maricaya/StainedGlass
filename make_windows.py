
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
SLIDE = config.get("slide", 0)  # Default to 0 if not specified in the config

# Function to check if FAI file exists
def check_fai(fasta):
    fai = f"{fasta}.fai"
    if not os.path.exists(fai):
        sys.exit(f"Input fasta must be indexed, try:\nsamtools faidx {fasta}")
    return fai

# Check if FAI file exists
FAI = check_fai(FASTA)

# Output and log file paths
output_bed = f"temp/{SM}.{W}.bed"
log_file = f"logs/make_windows.{SM}.{W}.log"

# Create the command
slide_param = f"-s {SLIDE}" if SLIDE > 0 else ""
command = f"bedtools makewindows -g {FAI} -w {W} {slide_param} > {output_bed}"

# Execute the command
try:
    os.makedirs(os.path.dirname(output_bed), exist_ok=True)
    os.makedirs(os.path.dirname(log_file), exist_ok=True)
    with open(log_file, "w") as log:
        subprocess.run(command, shell=True, check=True, stderr=log)
    print(f"Windows file created: {output_bed}")
except subprocess.CalledProcessError as e:
    print(f"Error running bedtools makewindows: {e}", file=sys.stderr)
    sys.exit(1)
