import subprocess
import sys
import os
import yaml
import runpy

# 加载配置文件 config/config.yaml
config_path = 'config/config.yaml'
if not os.path.exists(config_path):
    sys.exit(f"Configuration file not found: {config_path}")

with open(config_path, 'r') as config_file:
    config = yaml.safe_load(config_file)

# 从配置文件中读取参数
SM = config["sample"]
W = config["window"]
F = config["mm_f"]
S = config.get("mm_s", int(W / 5))
THREADS = config.get("alnthreads", 8)
MEM = config.get("mem", 8)

# 定义路径
input_bam = f"results/{SM}.{W}.{F}.sorted.bam"
output_tbl = f"temp/{SM}.{W}.{F}.tbl.gz"
log_file = f"logs/identity.{SM}.{W}.{F}.log"
sam_identity_script = "workflow/scripts/samIdentity.py"  # 脚本路径

# 确保输入 BAM 文件存在
if not os.path.exists(input_bam):
    sys.exit(f"Input BAM file not found: {input_bam}")

# 设置参数
sys.argv = [
    sam_identity_script,  # 模拟命令行执行脚本名
    '--threads', str(THREADS),
    '--matches', str(S),
    '--header',
    input_bam
]

# 执行脚本
try:
    os.makedirs(os.path.dirname(output_tbl), exist_ok=True)
    os.makedirs(os.path.dirname(log_file), exist_ok=True)

    # 使用 runpy 执行脚本
    runpy.run_path(sam_identity_script, run_name="__main__")
    
    print(f"Identity table created: {output_tbl}")
    
except Exception as e:
    print(f"Unexpected error: {str(e)}", file=sys.stderr)
    sys.exit(1)
