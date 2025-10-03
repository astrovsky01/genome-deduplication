### For creating the multi-species datasets, we need to compile training and dev samples for each species and put them into one dataset
### Read in each dataset, shuffle them so that the model isn't biased by sample order, and write to final files
### Usage: python compile_pretraining_data.py [in_folder] [category] [dedup status] [out_folder]
### Category: {human, primate, mammal, vertebrate, animal, eukaryote, organism}
### Dedup status: {raw, deduplicated/dedup}

###============================================================================
### Imports

import os
import random
import sys

###============================================================================
### Process input args

if len(sys.argv) < 3:
    print("Error: must supply input args for in_folder, category, dedup_status, and out_folder")
    print("Usage: python aggregate_samples.py [in_folder] [seed] [...text file of sample basenames to include]")
    exit()

in_folder = sys.argv[1]

try:
    random_seed = int(sys.argv[2])
except Exception as e:
    print(f"Error: {sys.argv[2]} is not a valid seed")
    exit()

if len(sys.argv) > 3:
    samples_file = sys.argv[3]
    if not os.path.isfile(samples_file):
        print(f"Error: could not find samples file at {samples_file}")
        exit()

if not os.path.isdir(in_folder):
    print(f"Error: could not find input folder at {in_folder}")
    exit()

###============================================================================
### Make output folder at in_folder/final

out_folder = os.path.join(in_folder, "final")
if not os.path.isdir(out_folder):
    os.makedirs(out_folder, exist_ok=True)

###============================================================================
### Global args
### Can change these to match filetype naming convention

filetype_suffix_table = {
        "train_sample": ".samples.train.txt",
        "train_coords": ".samples.train.bed",
        "dev_sample": ".samples.dev.txt",
        "dev_coords": ".samples.dev.bed"
}

###============================================================================
### Get all needed files for each sample

in_folder_files = os.listdir(in_folder)
samples = [f[:-12] for f in in_folder_files if f.endswith(".samples.bed")]
file_dict = {}
for sample in samples:
    sample_file_dict = {}
    for file_name, filetype_suffix in filetype_suffix_table.items():
        sample_file = f"{sample}{filetype_suffix}"
        if sample_file not in in_folder_files:
            print(f"Could not find {file_name} file for {sample}; Excluding this sample")
            break
        sample_file_dict[file_name] = os.path.join(in_folder, sample_file)
    if len(sample_file_dict.keys()) == len(filetype_suffix_table.keys()):
        file_dict[sample] = sample_file_dict

print(f"Found all needed files for {len(file_dict.keys())} samples")

###============================================================================
### Read in all files

sample_data = {sample: {} for sample in samples}
for sample in samples:
    print(f"Reading data for {sample}")
    for file_type, file_path in file_dict[sample].items():
        with open(file_path, 'r') as f:
            lines = f.readlines()
        sample_data[sample][file_type] = lines

###============================================================================
### Assemble data

# Append sample info to all coords data for better clarity when combining everything
print("Adding sample name to coordinates data")
for sample in samples:
    for file_type in ["train_coords", "dev_coords"]:
        sample_data[sample][file_type] = [f"{line[:-1]}\t{sample}\n" for line in sample_data[sample][file_type]]

# Compile all data of each type
print("Compiling datasets")
file_type_data = {}
for file_type in filetype_suffix_table.keys():
    file_type_lines = [line for sample in samples for line in sample_data[sample][file_type]]
    file_type_data[file_type] = file_type_lines

# Assert data lengths match
assert(len(file_type_data["train_sample"]) == len(file_type_data["train_coords"]))
assert(len(file_type_data["dev_sample"]) == len(file_type_data["dev_coords"]))

###============================================================================
### Shuffle data

print(f"Shuffling data w/ seed={random_seed}")

# Shuffle training samples
n_train_samples = len(file_type_data["train_sample"])
train_sample_order = [i for i in range(n_train_samples)]
random.Random(random_seed).shuffle(train_sample_order)
file_type_data["train_sample"] = [file_type_data["train_sample"][i] for i in train_sample_order]
file_type_data["train_coords"] = [file_type_data["train_coords"][i] for i in train_sample_order]

# Shuffle dev samples - is this necessary? We do it anyway just to be safe
n_dev_samples = len(file_type_data["dev_sample"])
dev_sample_order = [i for i in range(n_dev_samples)]
random.Random(random_seed).shuffle(dev_sample_order)
file_type_data["dev_sample"] = [file_type_data["dev_sample"][i] for i in dev_sample_order]
file_type_data["dev_coords"] = [file_type_data["dev_coords"][i] for i in dev_sample_order]

###============================================================================
### Write shuffled data to output folder

print("Writing files")

with open(os.path.join(out_folder, "train.txt"), 'w') as f:
    for line in file_type_data["train_sample"]:
        f.write(line)
with open(os.path.join(out_folder, "train.bed"), 'w') as f:
    for line in file_type_data["train_coords"]:
        f.write(line)
with open(os.path.join(out_folder, "dev.txt"), 'w') as f:
    for line in file_type_data["dev_sample"]:
        f.write(line)
with open(os.path.join(out_folder, "dev.bed"), 'w') as f:
    for line in file_type_data["dev_coords"]:
        f.write(line)

