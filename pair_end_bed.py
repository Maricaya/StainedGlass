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
SLIDE = config.get("slide", 0)

# Paths
input_tbl = f"temp/{SM}.{W}.{F}.tbl.gz"
fai_path = f"{os.path.abspath(config['fasta'])}.fai"
output_bed = f"results/{SM}.{W}.{F}.bed.gz"
output_full = f"results/{SM}.{W}.{F}.full.tbl.gz"
log_file = f"logs/pair_end_bed.{SM}.{W}.{F}.log"
refmt_script = "workflow/scripts/refmt.py"  # Correct path

# Verify that the input table file exists and is not empty
if not os.path.exists(input_tbl):
    sys.exit(f"Input table file not found: {input_tbl}")
elif os.path.getsize(input_tbl) == 0:
    sys.exit(f"Input table file is empty: {input_tbl}")

# Verify that the FAI file exists
if not os.path.exists(fai_path):
    sys.exit(f"FAI file not found: {fai_path}")

# Additional parameters
one_param = "--one" if SLIDE > 0 else ""

# Create the command
command = f"python {refmt_script} --window {W} --fai {fai_path} --full {output_full} {one_param} {input_tbl} {output_bed}"

# Execute the command
try:
    os.makedirs(os.path.dirname(output_bed), exist_ok=True)
    os.makedirs(os.path.dirname(output_full), exist_ok=True)
    os.makedirs(os.path.dirname(log_file), exist_ok=True)
    with open(log_file, "w") as log:
        result = subprocess.run(command, shell=True, check=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
        log.write(result.stderr.decode())
        log.write(result.stdout.decode())
        print(f"Pair-end BED and full table created: {output_bed}, {output_full}")
except subprocess.CalledProcessError as e:
    print(f"Error running pair-end bed reformatting: {e.stderr.decode()}", file=sys.stderr)
    sys.exit(1)
except Exception as e:
    print(f"Unexpected error: {str(e)}", file=sys.stderr)
    sys.exit(1)
