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
SLIDE = config.get("slide", 0)
N = config.pop("nbatch", 1)

# Determine the IDs dynamically based on sequence splitting
FAI = f"{FASTA}.fai"
if not os.path.exists(FAI):
    sys.exit(f"Input fasta must be indexed, try:\nsamtools faidx {FASTA}")

lengths = [math.ceil(int(line.split()[1]) / W) for line in open(FAI).readlines()]
N_WINDOWS = int(sum(lengths))
N = min(N, N_WINDOWS)
IDS = list(range(N))

# Function to check if bedtools is installed
check_tool("bedtools")

# Create bed files if they do not exist
for ID in IDS:
    bed_file = f"temp/{SM}.{W}.{ID}.bed"
    if not os.path.exists(bed_file):
        with open(bed_file, 'w') as bed:
            bed.write(f"{ID}\t0\t{W}\n")

# Create and execute the commands for each ID
for ID in IDS:
    bed_file = f"temp/{SM}.{W}.{ID}.bed"
    query_fasta = f"temp/{SM}.{W}.{ID}.query.fasta"
    log_file = f"logs/query_prep.{SM}.{W}.{ID}.log"

    # Check if bed file exists
    if not os.path.exists(bed_file):
        print(f"Bed file not found for ID {ID}: {bed_file}", file=sys.stderr)
        continue

    # Create the command
    command = f"bedtools getfasta -fi {FASTA} -bed {bed_file} > {query_fasta}"

    # Execute the command
    try:
        os.makedirs(os.path.dirname(query_fasta), exist_ok=True)
        os.makedirs(os.path.dirname(log_file), exist_ok=True)
        with open(log_file, "w") as log:
            result = subprocess.run(command, shell=True, check=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
            log.write(result.stderr.decode())
            log.write(result.stdout.decode())
            print(result.stderr.decode())
            print(result.stdout.decode())
            print(f"Query fasta file created for ID {ID}: {query_fasta}")
    except subprocess.CalledProcessError as e:
        print(f"Error running bedtools for ID {ID}: {e.stderr.decode()}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Unexpected error for ID {ID}: {str(e)}", file=sys.stderr)
        sys.exit(1)
