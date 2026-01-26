import os
import random
import argparse

parser = argparse.ArgumentParser(description="Get deduplicated samples based on control line counts.")
parser.add_argument('-f', '--deduped_path', type=str, required=True)
parser.add_argument('-t', '--control_train_lines', type=int, required=True)
parser.add_argument('-d', '--control_dev_lines', type=int, required=True)
parser.add_argument('-o', '--output_dir', type=str, default="control_samples")
parser.add_argument('-s', '--seed', type=int, default=42)
args = parser.parse_args()
# deduped_path="/scratch4/mschatz1/aostrov4/retraining_project/genome-deduplication/deduplicated_files"
# control_lines=os.path.join(deduped_path, "full_hprc_control/line_counts.txt")
# samples_list=["full_hprc", "full_hprc_99", "full_hprc_97"]

def get_sample_counts(sample_file_path):
    train_file=os.path.join(sample_file_path, "all_train.bed")
    dev_file=os.path.join(sample_file_path, "all_dev.bed")
    with open(train_file, 'r') as f:
        train_count = sum(1 for _ in f)
    with open(dev_file, 'r') as f:
        dev_count = sum(1 for _ in f)
    total_count = train_count + dev_count
    return train_count, dev_count, total_count


def get_deduped_lines(sample_dir_path, sample_output_dir, control_train_counts, control_dev_counts, seed=42):
    random.seed(seed)
    sample_name=os.path.basename(sample_dir_path)
    sample_train_counts, sample_dev_counts, sample_total_counts = get_sample_counts(sample_dir_path)
    # Sed is 1-based, so change range to accomodate
    train_lines = random.sample(range(1, control_train_counts + 1), sample_train_counts)
    dev_lines = random.sample(range(1, control_dev_counts + 1), sample_dev_counts)
    if not os.path.exists(sample_output_dir):
        os.makedirs(sample_output_dir)
    with open(os.path.join(sample_output_dir, f"{sample_name}_train_sampled.lines"), 'w') as train_out:
        for i in train_lines:
            train_out.write(str(i)+"\n")
    with open(os.path.join(sample_output_dir, f"{sample_name}_dev_sampled.lines"), 'w') as dev_out:
        for i in dev_lines:
            dev_out.write(str(i)+"\n")

if __name__ == "__main__":
    get_deduped_lines(
        os.path.join(args.deduped_path),
        os.path.join(args.deduped_path, args.output_dir),
        args.control_train_lines,
        args.control_dev_lines,
        seed=args.seed
    )