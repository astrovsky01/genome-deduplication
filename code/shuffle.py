### General script for shuffling an input file
### Usage: python shuffle.py [in file] [seed] [out file]

# Imports
import random
import sys

# Collect args
in_file = sys.argv[1]
random_seed = int(sys.argv[2])
out_file = sys.argv[3]

# Read input file
with open(in_file, 'r') as f:
    lines = f.readlines()

# Produce shuffled line order
rng = random.Random(random_seed)
line_order = [i for i in range(len(lines))]
rng.shuffle(line_order)

# Write lines in shuffled order
with open(out_file, 'w') as f:
    for i in line_order:
        f.write(lines[i])
