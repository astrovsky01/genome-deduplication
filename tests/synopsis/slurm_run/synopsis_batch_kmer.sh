#!/bin/bash
#SBATCH --job-name="synopsis-test"
#SBATCH --account="account-test"
#SBATCH --partition="partition-test"
#SBATCH --time="00:05:00"
#SBATCH --mem="80"
#SBATCH --cpus-per-task="2"
#SBATCH --ntasks="1"
#SBATCH --nodes="1"
#SBATCH --gres="gpu:1"
#SBATCH --output="slurm-%A_%a.out"
#SBATCH --error="slurm-%A_%a.err"
#SBATCH --mail-user="test@example.com"
#SBATCH --mail-type="END,FAIL"
#SBATCH --array=0-2%2

echo setup-from-test-config

## Using SLURM_ARRAY_TASK_ID to select the input file from the list
FASTA_LIST="/var/folders/1_/sg134bt167vbp_2v993ch7x40000gn/T/synopsis-tests.g9yog9/gca_inputs.txt"
FASTA_PATH=$(awk -v index="$SLURM_ARRAY_TASK_ID" 'NF && $0 !~ /^[[:space:]]*#/ { if (count++ == index) { print; exit } }' "$FASTA_LIST")
BASENAME=$(basename "$FASTA_PATH")
BASENAME="${BASENAME%.gz}"
BASENAME="${BASENAME%.*}"

"/var/folders/1_/sg134bt167vbp_2v993ch7x40000gn/T//synopsis-tests.g9yog9/bin/fake-dedup" "-e" "per_kmer" -o "./${BASENAME}" "$FASTA_PATH"
