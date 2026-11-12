# Synopsis

`synopsis` is the command-line entry point for running FASTA deduplication. It supports running a single dataset, running many files independently in sequence, running many files independently in parallel on a Slurm cluster, or merging/shuffling previously deduplicated files together.

## Requirements

- Bash 4+ or Zsh
- bedtools
- C++17 compiler and zlib

## Building

```sh
make
```

This installs the `synopsis` command to a standard writable directory on your
`PATH` when one is available. On Apple Silicon macOS, this is normally
`/opt/homebrew/bin`; otherwise it falls back to `~/.local/bin`.

If the fallback directory is used and is not already on your `PATH`, add it to
your shell profile:

```sh
export PATH="$HOME/.local/bin:$PATH"
```

You can then run `synopsis` from any directory:

```sh
synopsis --help
```

For a system-wide installation, use a prefix on the system `PATH`, such as
`/usr/local`:

```sh
make PREFIX=/usr/local
```

On macOS, the launcher uses Zsh when available because the system Bash 3.2 is
older than the Bash 4+ required by `synopsis.sh`.

## Usage

```
synopsis [kmer|sample] [options] <input file|list>
synopsis individual [kmer|sample] [options] <input file|list>
synopsis individual [kmer|sample] --config <config_file> [-y] [options] <input file|list>
synopsis balanced [kmer|sample] [options] <file_list> <indir>
```

### Modes

| Mode | Description |
|---|---|
| **Direct** (default) | Runs deduplication once on the given input. If the input is a `.txt` list of fasta paths, it's treated as a single combined dataset in order. |
| **Individual** | Runs deduplication independently on each file in a list, one after another. Each file gets its own output subdirectory. With `--config`, generates a Slurm array script instead. When --config is set, additionally setting -y will automatically submit the generated slurm script using sbatch. |
| **Shuffle** | Merges adjacent samples from a directory of previously, independently deduplicated files, shuffles them together, re-deduplicates, and generates new samples. |

### Choosing a mode and evaluation method

The evaluation method is optional and defaults to `kmer`. Use `sample` for per-sample threshold evaluation.

Examples:

```sh
Direct mode (default):   synopsis [kmer|sample] [options] <input file|list>
Individual mode:         synopsis individual [kmer|sample] [options] <input file|list>
Batch mode:              synopsis individual [kmer|sample] --config <config_file> [-y] [options] <input file|list>
Balanced mode:           synopsis balanced [kmer|sample] [options] <file_list> <indir>
```

## Input detection

The final positional argument is inspected by extension:

| Extension | Treated as |
|---|---|
| `.txt` | A list of fasta file paths, one per line |
| `.fa`, `.fna`, `.fasta`, `.fa.gz`, `.fna.gz`, `.fasta.gz` | A single fasta file |
| anything else (including `.fq`/`.fastq`) | Not supported |

When a `.txt` list is given to `individual` mode:

- Blank lines are skipped.
- Lines beginning with `#` are treated as comments and skipped.
- Every other line is assumed to be a **full, absolute path** to a fasta file. 

## Output directories

`-o <path>` sets the base output directory and is shared across all modes.

- In **direct** mode, `-o` is used exactly as given.
- In **individual** mode, each fasta file gets its own subdirectory under `-o`, named after that file's basename with its extension(s) stripped (including both parts of a `.gz`-compressed double extension, e.g. `sample1.fasta.gz` → `sample1`):

  ```
  <outdir>/<basename_without_extension>/
  ```

## Dedup options (shared across all modes)

| Flag | Description |
|---|---|
| `-d <float>` | Deduplication retention rate. Default 0.0|
| `-k <int>` | K-mer length |
| `-l <int>` | Standard sample length |
| `-m <int>` | Minimum sample length (defaults to standard sample length) |
| `-o <path>` | Output directory |
| `-p <path>` | Seen kmers set from previous run to resume |
| `-s <int>` | Save seen kmer set frequency |
| `-v <int>` | Number of overlapping bases between samples. Default is k-1. |
| `-r` | Retain information as scanning each contig. Recommended for large sample lengths. |
| `-f` | Generate final sample files |
| `--save_kmers_at_end` | |
| `--write_ambiguous_beds` | |
| `--write_ignored_beds` | |
| `--write_masked_beds` | |
| `--seed <int>` | |
| `-h`, `--help` | Show usage |

## Shared Post-processing parameters
| Flag | Description |
|---|---|
| `-t <float>`  | Specify the amount of samples to reserve for test set relative to the total number of samples. If not specified, no separate test set will be generated | 
| `--dev_rate <float>` | Specify the fraction of samples to reserve for dev set. Default 0.1 |
| `--max_dev <int>` |  Specify the maximum number of samples to reserve for dev set. If this count is lower than the dev_rate fraction, the dev set will be capped at this number. |

## Shuffle mode options

```
synopsis shuffle [kmer|sample] [options] <file_list> <indir>
```

| Flag | Description |
|---|---|
| `-i` | Generate final ignored bed files |
| `-z` | Generate final masked bed files |
| `-x` | Add originally ignored regions before merging |

`<file_list>` and `<indir>` specify the previously, independently deduplicated files to merge, shuffle, and re-deduplicate. The evaluation method chosen (default `kmer`) governs the *inner* re-deduplication step. Nested modes are not allowed — you cannot run `synopsis shuffle individual`.

## Batch generation via individual mode

```
synopsis individual [kmer|sample] --config <config_file> [-y] [options] <input file|list>
```

When `individual` receives `--config`, it resolves the input into a list of fasta files and generates a Slurm script with one array task per file. Each array task deduplicates a single assigned fasta file, writing to its own output subdirectory as described above. The script is written locally; it is not executed by `synopsis` unless `-y` is supplied in an environment with Slurm.

### `--config <file>` (required)

Points to a Slurm configuration file (see format below) describing `#SBATCH` settings and any environment setup needed on the cluster (module loads, conda/venv activation, etc).

### `-y`

If given, after generating the Slurm script, `synopsis` attempts to submit it immediately via `sbatch`. Without `-y`, the script is only written to disk.

### Generated script

The output script is written to the **caller's current working directory** and named:

```
synopsis_batch_<eval_method>.sh
```

e.g. `synopsis_batch_kmer.sh` or `synopsis_batch_sample.sh`.

The `--array` range is computed automatically based on the number of resolved input files. Set `max_concurrent_jobs` to a positive integer to limit how many array tasks run concurrently; it is emitted as the Slurm `%<number>` array limit.

### Slurm config file format

The config file has a Slurm settings section and an optional setup section:

```
###SLURM VARIABLES###
#_____________________________
job-name=""
account=""
partition=""
time=""
mem=""
cpus-per-task=""
ntasks=""
nodes=""
gres=""
output=""
error=""
mail-user=""
mail-type=""
max_concurrent_jobs=""

###LOCAL SETUP (ENVIRONMENTS, MODULES, ETC)###
#_____________________________
module load anaconda3
source activate my_env
```

- Put Slurm options in `key="value"` form under `###SLURM VARIABLES###`.
- Empty values are ignored.
- `max_concurrent_jobs` limits concurrent array tasks and must be a positive integer.
- Put optional module loads or environment setup below `###LOCAL SETUP (ENVIRONMENTS, MODULES, ETC)###`.
- Setup text is copied into the generated script before the deduplication command.

# Post-processing

To generate the final sample files, the following commands will generate the same output files


```sh
<synopsis_root_dir>/code/get_samples.sh -i <deduplicated_data_dir> [-p <float>] [-m <int>] [-s <int>] [-t]
```

#### Options

- `-i <deduplicated_data_dir>` — Path to the deduplicated data directory.
- `-p <float>` — Specify the fraction of samples to reserve for dev set. Default 0.1
- `-m <int>` —  Specify the maximum number of samples to reserve for dev set. If this count is lower than the dev_rate fraction, the dev set will be capped at this number.
- `-s <int>` — Seed for shuffling.
- `-t` — Disable `.txt` output.