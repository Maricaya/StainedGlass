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
input_bed = f"temp/{SM}.{W}.bed"
output_bed_pattern = f"temp/{SM}.{W}.{{ID}}.bed"
log_file = f"logs/split_windows.{SM}.{W}.log"
script = "./workflow/scripts/batch_bed_files.py"

# Generate output bed file paths
N = config.get("nbatch", 1)
output_beds = [output_bed_pattern.format(ID=i) for i in range(N)]

# Create the command
output_beds_str = " ".join(output_beds)
command = f"python {script} {input_bed} --outputs {output_beds_str}"

# Execute the command
try:
    os.makedirs(os.path.dirname(input_bed), exist_ok=True)
    os.makedirs(os.path.dirname(log_file), exist_ok=True)
    with open(log_file, "w") as log:
        result = subprocess.run(command, shell=True, check=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
    print(f"Split BED files created: {output_beds}")
except subprocess.CalledProcessError as e:
    print(f"Error running batch_bed_files.py: {e}", file=sys.stderr)
    print(e.stderr.decode(), file=sys.stderr)
    sys.exit(1)
