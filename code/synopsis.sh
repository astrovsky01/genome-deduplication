#!/usr/bin/env bash

# --- Shell compatibility: allow running under bash or zsh ---
if [ -n "${ZSH_VERSION:-}" ]; then
    emulate -L bash   # 0-indexed arrays, BASH_SOURCE, word-splitting, etc.
fi

set -euo pipefail

# Require Bash >= 4 for associative arrays when actually running under bash
# (zsh's bash emulation supports associative arrays regardless of version)
if [ -z "${ZSH_VERSION:-}" ] && [ -n "${BASH_VERSION:-}" ]; then
    if ((BASH_VERSINFO[0] < 4)); then
        echo "Error: this script requires Bash >= 4 or Zsh (found Bash ${BASH_VERSION})." >&2
        echo "On macOS, install a newer bash via: brew install bash" >&2
        exit 1
    fi
fi

primary_binary="dedup6"
ROOTDIR="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"
DEDUP_BIN="${DEDUP_BIN:-$ROOTDIR/$primary_binary}"
SHUFFLE_SCRIPT="${SHUFFLE_SCRIPT:-$ROOTDIR/shuffle_method.sh}"
FINAL_SAMPLES_SCRIPT="${FINAL_SAMPLES_SCRIPT:-$ROOTDIR/get_samples.sh}"

readonly MODE_WORDS=(individual balanced)
readonly VALID_EVAL_METHODS=(kmer sample)

declare -A EVAL_METHOD_MAP=(
    [kmer]="per_kmer"
    [sample]="per_sample_threshold"
)
readonly EVAL_METHOD_MAP

readonly SLURM_KEYS=(job-name account partition time mem cpus-per-task ntasks nodes gres output error mail-user mail-type max_concurrent_jobs)

usage() {
    echo "Usage:"
    echo "  Direct mode (default):   $0 [kmer|sample] [options] <input file|list>"
    echo "  Individual mode:         $0 individual [kmer|sample] [options] <input file|list>"
    echo "  Batch mode:              $0 individual [kmer|sample] --config <config_file> [-y] [options] <input file|list>"
    echo "  Balanced mode:           $0 balanced [kmer|sample] [options] <file_list> <indir>"
    echo ""
    echo "Dedup options (shared across all modes):"
    echo "  -d <float>  Duplicate retention rate. Default 0.0 (strict)"
    echo "  -k <int>    K-mer length"
    echo "  -l <int>    Sample length"
    echo "  -m <int>    Minimum sample length. Allows smaller samples than -l to be retained. Defaults to sample length"
    echo "  -o <path>   Output directory path"
    echo "  -p <path>   Seen kmers file path"
    echo "  -s <int>    Save kmer set every <int> sequences"
    echo "  -v <int>    Overlap samples by specified base count"
    echo "  -r          Retain information as scanning each contig. Recommended for large sample lengths."
    echo "  --save_kmers_at_end"
    echo "  --write_ambiguous_beds (not used in balanced)"
    echo "  --write_ignored_beds (not used in balanced)"
    echo "  --write_masked_beds (not used in balanced)"
    echo "  --seed <int>"
    echo "  -h, --help"
    echo ""
    echo "Mode-specific options:"
    echo " Balanced"
    echo "  -i           generate final ignored bed files for each original input file"
    echo "  -z           generate final masked bed files for each original input file"
    echo "  -x           add originally ignored files to secondary deuplication step after shuffling. Must have run previous step with --write-ignored-beds."
    echo " Individual/batch:"
    echo "  --config <file>  generate a batch script with this config"
    echo "  -y         submit batch script immediately"
    echo ""
    echo "Post-processing options:"
    echo "  -f          Generate final training bed and fasta files from all samples. Creates <output dir>/all_train.bed, <output dir>/all_train.fasta, and <output dir>/all_dev.bed, all_dev.fasta"
    echo "  -t <int>    Specify the number of samples to reserve for test set. If not specified, no separate test set will be generated."
    echo "  --dev_rate <float>  Specify the fraction of samples to reserve for dev set. Default 0.1"
    echo "  --max_dev <int>  Specify the maximum number of samples to reserve for dev set. If this count is lower than the dev_rate fraction, the dev set will be capped at this number."
    echo ""

}

die() {
    echo "Error: $*" >&2
    exit 1
}

is_mode_word() {
    local tok="$1" w
    for w in "${MODE_WORDS[@]}"; do [[ "$tok" == "$w" ]] && return 0; done
    return 1
}

is_eval_method() {
    local tok="$1" w
    for w in "${VALID_EVAL_METHODS[@]}"; do [[ "$tok" == "$w" ]] && return 0; done
    return 1
}

has_valid_input_extension() {
    local tok="$1"
    case "$tok" in
        *.txt) return 0 ;;
        *.fa|*.fna|*.fasta) return 0 ;;
        *.fa.gz|*.fna.gz|*.fasta.gz) return 0 ;;
        *) return 1 ;;
    esac
}

is_valid_input_token() {
    local tok="$1"
    [[ -f "$tok" ]] && has_valid_input_extension "$tok"
}

basename_noext() {
    local path="$1" base
    base="${path##*/}"
    if [[ "$base" == *.gz ]]; then
        base="${base%.gz}"
    fi
    base="${base%.*}"
    echo "$base"
}

detect_input_type() {
    local path="$1"
    case "$path" in
        *.txt) echo "list" ;;
        *) echo "fasta" ;;
    esac
}

resolve_file_list() {
    local path="$1"
    local out_ref="$2"
    local -a result=()
    if [[ "$(detect_input_type "$path")" == "list" ]]; then
        [[ -f "$path" ]] || die "input list not found: $path"
        local line
        while IFS= read -r line || [[ -n "$line" ]]; do
            [[ -z "$line" ]] && continue
            [[ "$line" == \#* ]] && continue
            result+=("$line")
        done < "$path"
    else
        result=("$path")
    fi
    eval "$out_ref=(\"\${result[@]}\")"
}

mode="direct"
eval_method=""

if [[ $# -gt 0 ]] && is_mode_word "$1"; then
    mode="$1"
    shift
fi

if [[ $# -gt 0 ]] && is_eval_method "$1"; then
    eval_method="$1"
    shift
fi

if [[ -z "$eval_method" ]]; then
    eval_method="kmer"
fi

if [[ "$mode" == "shuffle" ]] && [[ "$eval_method" != "kmer" && "$eval_method" != "sample" ]]; then
    die "shuffle mode only accepts kmer or sample as its eval method"
fi

if [[ $# -gt 0 ]]; then
    next_tok="$1"
    if [[ "$next_tok" != -* ]] && ! is_valid_input_token "$next_tok"; then
        die "unrecognized argument: $next_tok (expected a flag, or an existing input file with a valid extension)"
    fi
fi

mapped_eval="${EVAL_METHOD_MAP[$eval_method]}"

dedup_param=""
kmer_size=""
sample_len=""
min_sample_len=""
output_dir=""
seen_kmers=""
save_every=""
overlap=""
retain_info=0
save_kmers_at_end=0
write_ambiguous_beds=0
write_ignored_beds=0
write_masked_beds=0
seed_val=""

shuffle_i=0
shuffle_z=0
shuffle_x=0

config_file=""
submit_now=0
max_concurrent_jobs=""

final_samples_out=0
test_set_count=""
dev_rate="0.1"
max_dev=""

args=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        --save_kmers_at_end)
            save_kmers_at_end=1
            shift
            ;;
        --write_ambiguous_beds)
            write_ambiguous_beds=1
            shift
            ;;
        --write_ignored_beds)
            write_ignored_beds=1
            shift
            ;;
        --write_masked_beds)
            write_masked_beds=1
            shift
            ;;
        --seed)
            seed_val="$2"
            shift 2
            ;;
        --overlap)
            overlap="$2"
            shift 2
            ;;
        --config)
            config_file="$2"
            shift 2
            ;;
        --dev_rate)
            dev_rate="$2"
            shift 2
            ;;
        --max_dev)
            max_dev="$2"
            shift 2
            ;;
        --help)
            usage
            exit 0
            ;;
        --*)
            die "unknown long option: $1"
            ;;
        *)
            args+=("$1")
            shift
            ;;
    esac
done
set -- "${args[@]}"

while getopts ":d:k:l:m:o:p:s:v:t:rizxyhf" opt; do
    case "$opt" in
        d) dedup_param="$OPTARG" ;;
        k) kmer_size="$OPTARG" ;;
        l) sample_len="$OPTARG" ;;
        m) min_sample_len="$OPTARG" ;;
        o) output_dir="$OPTARG" ;;
        p) seen_kmers="$OPTARG" ;;
        s) save_every="$OPTARG" ;;
        t) test_set_count="$OPTARG" ;;
        v) overlap="$OPTARG" ;;
        f) final_samples_out=1;;
        r) retain_info=1 ;;
        i) shuffle_i=1 ;;
        z) shuffle_z=1 ;;
        x) shuffle_x=1 ;;
        y) submit_now=1 ;;
        h) usage; exit 0 ;;
        \?) die "unknown option: -$OPTARG" ;;
        :) die "option -$OPTARG requires an argument" ;;
    esac
done
shift $((OPTIND - 1))

if [[ "$mode" == "individual" && -n "$config_file" ]]; then
    mode="batch"
fi

if [[ "$submit_now" -eq 1 && "$mode" != "batch" ]]; then
    die "-y requires --config in individual mode"
fi

if [[ "$mode" == "batch" ]]; then
    [[ -n "$config_file" ]] || die "batch mode requires --config <file>"
    [[ -f "$config_file" ]] || die "config file not found: $config_file"
fi

if [[ "$mode" == "balanced" ]]; then
    [[ $# -eq 2 ]] || die "balanced mode requires exactly two positional arguments: <file_list> <indir>"
else
    [[ $# -eq 1 ]] || die "expected exactly one positional input argument (file or list), got $#"
fi

build_dedup_args() {
    local out_ref="$1"
    local -a result=()
    result+=(-e "$mapped_eval")
    [[ -n "$dedup_param" ]] && result+=(-d "$dedup_param")
    [[ -n "$kmer_size" ]] && result+=(-k "$kmer_size")
    [[ -n "$sample_len" ]] && result+=(-l "$sample_len")
    [[ -n "$min_sample_len" ]] && result+=(-m "$min_sample_len")
    [[ -n "$seen_kmers" ]] && result+=(-p "$seen_kmers")
    [[ -n "$save_every" ]] && result+=(-s "$save_every")
    [[ -n "$overlap" ]] && result+=(-v "$overlap")
    [[ "$retain_info" -eq 1 ]] && result+=(-r)
    [[ "$save_kmers_at_end" -eq 1 ]] && result+=(--save_kmers_at_end)
    [[ "$write_ambiguous_beds" -eq 1 ]] && result+=(--write_ambiguous_beds)
    [[ "$write_ignored_beds" -eq 1 ]] && result+=(--write_ignored_beds)
    [[ "$write_masked_beds" -eq 1 ]] && result+=(--write_masked_beds)
    [[ -n "$seed_val" ]] && result+=(--seed "$seed_val")
    eval "$out_ref=(\"\${result[@]}\")"
}

build_final_file_args() {
    local out_ref="$1"
    local -a result=()
    [[ -n "$test_set_count" ]] && result+=(-t "$test_set_count")
    [[ -n "$dev_rate" ]] && result+=(-p "$dev_rate")
    [[ -n "$max_dev" ]] && result+=(-m "$max_dev")
    [[ -n "$seed_val" ]] && result+=(-s "$seed_val")
    eval "$out_ref=(\"\${result[@]}\")"
}

parse_slurm_config() {
    local cfg="$1"
    local sbatch_ref="$2"
    local freetext_ref="$3"
    local -a sbatch_result=()
    local -a freetext_result=()

    local section="" line key val
    while IFS= read -r line || [[ -n "$line" ]]; do
        if [[ "$line" == "###SLURM VARIABLES###" ]]; then
            section="slurm"
            continue
        fi
        if [[ "$line" == "###LOCAL SETUP (ENVIRONMENTS, MODULES, ETC)###" ]]; then
            section="freetext"
            continue
        fi
        if [[ "$line" == \#* ]]; then
            continue
        fi

        if [[ "$section" == "slurm" ]]; then
            [[ -z "$line" ]] && continue
            [[ "$line" != *=* ]] && continue
            key="${line%%=*}"
            val="${line#*=}"
            val="${val%\"}"
            val="${val#\"}"
            [[ -z "$val" ]] && continue

            local known=0 k
            for k in "${SLURM_KEYS[@]}"; do
                [[ "$k" == "$key" ]] && known=1 && break
            done
            [[ "$known" -eq 1 ]] || die "unknown slurm config key: $key"

            if [[ "$key" == "max_concurrent_jobs" ]]; then
                [[ "$val" =~ ^[1-9][0-9]*$ ]] || die "max_concurrent_jobs must be a positive integer"
                max_concurrent_jobs="$val"
                continue
            fi

            sbatch_result+=("#SBATCH --${key}=\"${val}\"")
        elif [[ "$section" == "freetext" ]]; then
            freetext_result+=("$line")
fi
    done < "$cfg"
    eval "$sbatch_ref=(\"\${sbatch_result[@]}\")"
    eval "$freetext_ref=(\"\${freetext_result[@]}\")"
}

## BALANCED MODE
if [[ "$mode" == "balanced" ]]; then
    file_list_arg="$1"
    indir_arg="$2"

    shuffle_args=()
    [[ -n "$dedup_param" ]] && shuffle_args+=(-d "$dedup_param")
    [[ -n "$kmer_size" ]] && shuffle_args+=(-k "$kmer_size")
    [[ -n "$sample_len" ]] && shuffle_args+=(-l "$sample_len")
    [[ -n "$min_sample_len" ]] && shuffle_args+=(-m "$min_sample_len")
    [[ -n "$output_dir" ]] && shuffle_args+=(-o "$output_dir")
    [[ -n "$seen_kmers" ]] && shuffle_args+=(-p "$seen_kmers")
    [[ -n "$save_every" ]] && shuffle_args+=(-s "$save_every")
    [[ -n "$seed_val" ]] && shuffle_args+=(--seed "$seed_val")
    [[ -n "$overlap" ]] && shuffle_args+=(-v "$overlap")
    [[ "$retain_info" -eq 1 ]] && shuffle_args+=(-r)
    [[ "$shuffle_i" -eq 1 ]] && shuffle_args+=(-i)
    [[ "$shuffle_z" -eq 1 ]] && shuffle_args+=(-z)
    [[ "$shuffle_x" -eq 1 ]] && shuffle_args+=(-x)
    shuffle_args+=(-e "$mapped_eval")


    # exec "$SHUFFLE_SCRIPT" "${shuffle_args[@]}" "$file_list_arg" "$indir_arg"
    "$SHUFFLE_SCRIPT" "${shuffle_args[@]}" "$file_list_arg" "$indir_arg"
    # printf '# exec'
    # printf ' %q' "$SHUFFLE_SCRIPT" "${shuffle_args[@]}" "$file_list_arg" "$indir_arg"
    # printf '\n'


## INDIVIDUAL RUN
elif [[ "$mode" == "individual" ]]; then
    input_arg="$1"
    resolve_file_list "$input_arg" fasta_list

    [[ ${#fasta_list[@]} -gt 0 ]] || die "no input files found"

    build_dedup_args dedup_args

    for fasta in "${fasta_list[@]}"; do
        [[ -f "$fasta" ]] || die "input file not found: $fasta"
        sub="$(basename_noext "$fasta")"
        this_out="${output_dir:-.}/${sub}"
        "$DEDUP_BIN" "${dedup_args[@]}" -o "$this_out" "$fasta"
    done


## BATCH MODE
elif [[ "$mode" == "batch" ]]; then
    input_arg="$1"
    input_path="$(cd "$(dirname "$input_arg")" && pwd)/${input_arg##*/}"
    resolve_file_list "$input_arg" fasta_list

    [[ ${#fasta_list[@]} -gt 0 ]] || die "no input files found"

    build_dedup_args dedup_args
    if [[ "$final_samples_out" -eq 1 ]]; then
        build_final_file_args final_args
    fi

    parse_slurm_config "$config_file" sbatch_lines freetext_lines

    out_script="synopsis_batch_${eval_method}.sh"

    max_idx=$(( ${#fasta_list[@]} - 1 ))
    array_range="0-${max_idx}"
    [[ -n "$max_concurrent_jobs" ]] && array_range+="%${max_concurrent_jobs}"

    {
        echo "#!/bin/bash"
        for line in "${sbatch_lines[@]}"; do
            echo "$line"
        done
        echo "#SBATCH --array=${array_range}"
        echo ""

        for line in "${freetext_lines[@]}"; do
            echo "$line"
        done
        echo ""

        if [[ "$input_arg" == *.txt ]]; then
            echo "FASTA_LIST=\"$input_path\""
            echo 'FASTA_PATH=$(awk -v idx="$SLURM_ARRAY_TASK_ID" '\''NF && $0 !~ /^[[:space:]]*#/ { if (count++ == idx) { print; exit } }'\'' "$FASTA_LIST")'
        else
            echo "FASTA_PATH=\"$input_path\""
        fi

        echo 'BASENAME=$(basename "$FASTA_PATH")'
        echo 'BASENAME="${BASENAME%.gz}"'
        echo 'BASENAME="${BASENAME%.*}"'
        echo ""

        echo 'OUTPUT_DIR="'"${output_dir:-.}"'/${BASENAME}"'

        echo ""
        printf '"%s"' "$DEDUP_BIN"
        for a in "${dedup_args[@]}"; do
            printf ' "%s"' "$a"
        done
        echo ' -o "$OUTPUT_DIR" "$FASTA_PATH"'

        if [[ "$final_samples_out" -eq 1 ]]; then
            echo "FINAL_SAMPLES_SCRIPT=\"$FINAL_SAMPLES_SCRIPT\""
            echo -n '"$FINAL_SAMPLES_SCRIPT" "$OUTPUT_DIR"'
            for a in "${final_args[@]}"; do
                echo -n " \"$a\""
            done
            echo
        fi
    } > "$out_script"

    chmod +x "$out_script"
    echo "Generated: $out_script"

    if [[ "$submit_now" -eq 1 ]]; then
        sbatch "$out_script"
    fi

## STANDARD MODE
else
    input_arg="$1"
    build_dedup_args dedup_args
    dedup_opts=()
    [[ -n "$output_dir" ]] && dedup_opts+=(-o "$output_dir")
    # exec "$DEDUP_BIN" "${dedup_args[@]}" "${dedup_opts[@]}" "$input_arg"
    "$DEDUP_BIN" "${dedup_args[@]}" "${dedup_opts[@]}" "$input_arg"

fi


# FINAL SAMPLE GENERATION
if [[ "$final_samples_out" -eq 1 ]]; then
    if [[ "$mode" == "individual" ]]; then
    # Make final files for each subdirectory in the output directory
        for subdir in "$output_dir"/*/; do
            if [[ -d "$subdir" ]]; then
                build_final_file_args final_args
                # exec "$FINAL_SAMPLES_SCRIPT" $subdir "${final_args[@]}"
                "$FINAL_SAMPLES_SCRIPT" $subdir "${final_args[@]}"
            fi
        done
    elif [[ "$mode" == "balanced" ]]; then
        build_final_file_args final_args

        mv "${indir_arg}/final_files/basename_fasta_match.txt" \
        "${indir_arg}/final_files/shuffled_basename_fasta_match.txt"

        # Rebuild the mapping file from each individual directory
        for fasta in "${indir_arg}"/*.merged.bed; do
            individual_dir="${fasta%.merged.bed}"
            cat "${individual_dir}/basename_fasta_match.txt" \
                >> "${indir_arg}/final_files/basename_fasta_match.txt"
        done
        "$FINAL_SAMPLES_SCRIPT" ${indir_arg}/final_files "${final_args[@]}"
    elif [[ "$mode" == "direct" ]]; then
        build_final_file_args final_args

        input_dir="${output_dir:-dedup_out}"
        echo "$FINAL_SAMPLES_SCRIPT" -i "$input_dir" "${final_args[@]}"
        "$FINAL_SAMPLES_SCRIPT" -i "$input_dir" "${final_args[@]}"
    fi
fi
