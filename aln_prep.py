import subprocess
import sys
import os
import yaml
import math

# Function to check if a tool is installed
def check_tool(tool_name):
    try:
        result = subprocess.run([tool_name, "--version"], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        print(result.stdout.decode())
    except subprocess.CalledProcessError:
        sys.exit(f"Error: {tool_name} is not installed or not available in the PATH. Please install {tool_name} and try again.")
    except FileNotFoundError:
        sys.exit(f"Error: {tool_name} is not found in the PATH. Please install {tool_name} and try again.")

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
S = config.get("mm_s", int(W / 5))
MAP_PARAMS = "-ax ava-ont " if W < 500 else "-ax map-ont --secondary=no"
SLIDE = config.get("slide", 0)

# Determine input reference based on SLIDE value
input_ref = f"temp/{SM}.{W}.ref_{0}.split.fasta" if SLIDE > 0 else f"results/{SM}.{W}.fasta"

# Function to check if FAI file exists
def check_fai(fasta):
    fai = f"{fasta}.fai"
    if not os.path.exists(fai):
        sys.exit(f"Input fasta must be indexed, try:\nsamtools faidx {fasta}")
    return fai

# Check if minimap2 is installed
check_tool("minimap2")

# Check if FAI file exists
FAI = check_fai(FASTA)

# Determine the REF_IDs dynamically based on sequence splitting
lengths = [math.ceil(int(line.split()[1]) / W) for line in open(FAI).readlines()]
names = [line.split()[0] for line in open(FAI).readlines()]
REF_IDS = [f"ref_{i}" for i in range(len(names))] if SLIDE > 0 else ["ref_0"]

# Create and execute the command for each REF_ID
for REF_ID in REF_IDS:
    split_ref_index = f"temp/{SM}.{W}.{F}.{REF_ID}.fasta.mmi"
    log_file = f"logs/aln_prep.{SM}.{W}.{F}.{REF_ID}.log"

    # Create the command
    command = f"minimap2 -f {F} -s {S} {MAP_PARAMS} -d {split_ref_index} {input_ref}"

    # Execute the command
    try:
        os.makedirs(os.path.dirname(split_ref_index), exist_ok=True)
        os.makedirs(os.path.dirname(log_file), exist_ok=True)
        with open(log_file, "w") as log:
            result = subprocess.run(command, shell=True, check=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
            log.write(result.stderr.decode())
            log.write(result.stdout.decode())
            print(result.stderr.decode())
            print(result.stdout.decode())
            print(f"Index file created: {split_ref_index}")
    except subprocess.CalledProcessError as e:
        print(f"Error running minimap2 for REF_ID {REF_ID}: {e.stderr.decode()}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Unexpected error for REF_ID {REF_ID}: {str(e)}", file=sys.stderr)
        sys.exit(1)
