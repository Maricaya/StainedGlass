import subprocess
import sys
import os
import yaml
import math

# Function to check if a tool is installed
def check_tool(tool_name):
    try:
        result = subprocess.run([tool_name, "--version"], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        #print(result.stdout.decode())
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
ALN_T = config.get("alnthreads", 4)
MAP_PARAMS = "-ax ava-ont " if W < 500 else "-ax map-ont --secondary=no"
MEM = config.get("mem", 2)
SLIDE = config.get("slide", 0)

# Function to check if FAI file exists
def check_fai(fasta):
    fai = f"{fasta}.fai"
    if not os.path.exists(fai):
        sys.exit(f"Input fasta must be indexed, try:\nsamtools faidx {fasta}")
    return fai

# Check if minimap2 and samtools are installed
check_tool("minimap2")
check_tool("samtools")

# Check if FAI file exists
FAI = check_fai(FASTA)

# Determine the IDs and REF_IDs dynamically based on sequence splitting
lengths = [math.ceil(int(line.split()[1]) / W) for line in open(FAI).readlines()]
names = [line.split()[0] for line in open(FAI).readlines()]
REF_IDS = [f"ref_{i}" for i in range(len(names))] if SLIDE > 0 else ["ref_0"]
N_WINDOWS = int(sum(lengths))
N = config.pop("nbatch", 1)
N = min(N, N_WINDOWS)
print("===== N", N)
IDS = list(range(N))

# Create and execute the commands for each ID and REF_ID
for REF_ID in REF_IDS:
    for ID in IDS:
        query_fasta = f"temp/{SM}.{W}.{ID}.query.fasta"
        
        # Check if the query fasta file exists
        if not os.path.exists(query_fasta):
            #print(f"Query fasta file not found for ID {ID}: {query_fasta}", file=sys.stderr)
            continue
        
        split_ref = f"temp/{SM}.{W}.{F}.{REF_ID}.fasta.mmi"
        output_bam = f"temp/{SM}.{W}.{F}.{ID}.{REF_ID}.bam"
        intermediate_sam = f"temp/{SM}.{W}.{F}.{ID}.{REF_ID}.sam"
        log_file_minimap2 = f"logs/aln_minimap2.{SM}.{W}.{F}.{ID}.{REF_ID}.log"
        log_file_samtools = f"logs/aln_samtools.{SM}.{W}.{F}.{ID}.{REF_ID}.log"

        # Create the minimap2 command
        minimap2_command = f"minimap2 -t {ALN_T} -f {F} -s {S} {MAP_PARAMS} --dual=yes --eqx {split_ref} {query_fasta} 2> {log_file_minimap2} > /dev/null"

        # Create the samtools command
        samtools_command = f"samtools sort -m {MEM}G -@ {ALN_T} -o {output_bam} {intermediate_sam} 2> {log_file_samtools} > /dev/null"

        # Execute the minimap2 command
        try:
            os.makedirs(os.path.dirname(intermediate_sam), exist_ok=True)
            result = subprocess.run(minimap2_command, shell=True, check=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
            #print(result.stderr.decode())
        except subprocess.CalledProcessError as e:
            print(f"Error running minimap2 for ID {ID} and REF_ID {REF_ID}: {e.stderr.decode()}", file=sys.stderr)
            sys.exit(1)
        except Exception as e:
            print(f"Unexpected error for ID {ID} and REF_ID {REF_ID}: {str(e)}", file=sys.stderr)
            sys.exit(1)

        # Execute the samtools command
        try:
            os.makedirs(os.path.dirname(output_bam), exist_ok=True)
            result = subprocess.run(samtools_command, shell=True, check=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
            #print(result.stderr.decode())
            print(f"BAM file created for ID {ID} and REF_ID {REF_ID}: {output_bam}")
        except subprocess.CalledProcessError as e:
            print(f"Error running samtools for ID {ID} and REF_ID {REF_ID}: {e.stderr.decode()}", file=sys.stderr)
            sys.exit(1)
        except Exception as e:
            print(f"Unexpected error for ID {ID} and REF_ID {REF_ID}: {str(e)}", file=sys.stderr)
            sys.exit(1)
