from matplotlib import pyplot as plt
import numpy as np
from matplotlib_venn import venn3, venn3_circles
import argparse

# input file format:
"""
Savvy   18423   0.12
CNVKit  3286 0.12
GATK    15196   0.12
Savvy-GATK  757  0.12
Savvy-CNVKit    133    0.12
CNVKit-GATK 22  0.12
Triple  18   0.12
"""

def parse_args():
    # input file
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', type=str, help='input file')
    parser.add_argument('-o', type=str, help='output file')
    args = parser.parse_args()
    return args

def load_data(input_file):
    # load the data
    data = {}
    with open(input_file, 'r') as f:
        for line in f:
            line = line.strip().split('\t')
            data[line[0]] = float(line[1])
    assert 'Savvy' in data.keys(), 'Savvy not in input file'
    assert 'CNVKit' in data.keys(), 'CNVKit not in input file'
    assert 'GATK' in data.keys(), 'GATK not in input file'
    assert 'Savvy-GATK' in data.keys(), 'Savvy-GATK not in input file'
    assert 'Savvy-CNVKit' in data.keys(), 'Savvy-CNVKit not in input file'
    assert 'CNVKit-GATK' in data.keys(), 'CNVKit-GATK not in input file'
    assert 'Triple' in data.keys(), 'Triple not in input file'
    return data

def main():
    args = parse_args()
    data = load_data(args.i)
    # Make a Basic Venn
    subsets = (data['GATK'], data['Savvy'], data['Savvy-GATK'], data['CNVKit'], data['CNVKit-GATK'], data['Savvy-CNVKit'], data['Triple'])
    v = venn3(subsets=subsets, set_labels = ('gCNV', 'Savvy', 'CNVKit'))    
    # Add title and annotation
    plt.title("Sample Venn diagram")
    # Show it
    plt.savefig(args.o,dpi=300)

if __name__ == "__main__":
    main()