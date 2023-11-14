import argparse
import matplotlib.pyplot as plt
from upsetplot import from_memberships, plot

# input file format:
"""
Savvy: 18423
CNVKit: 3286
GATK: 15196
Savvy-GATK: 757
Savvy-CNVKit: 133
CNVKit-GATK: 22
Triple: 18
"""

# function for getting arguments
def get_args():
    parser = argparse.ArgumentParser(description='Plot the sizes of calls from bed file(s)')
    parser.add_argument('-i', '--inputs', action='append', required=True, help='Input file')
    parser.add_argument('-l', '--labels', action='append', required=True, help='Labels for input file(s)')
    parser.add_argument('-o', '--outputs', action='append', required=True, help='Output file(s)')
    parser.add_argument('--logscale', default=1, help='Should log scale be used 0 or 1', type=int)
    parser.add_argument('--percent', default=False, help='Use floats instead of ints', action='store_true')
    args = parser.parse_args()
    return args

def load_data(input_file):
    # load the data
    data = {}
    with open(input_file, 'r') as f:
        for line in f:
            line = line.strip().split(': ')
            data[line[0]] = float(line[1])
    assert 'Savvy' in data.keys(), 'Savvy not in input file'
    assert 'CNVKit' in data.keys(), 'CNVKit not in input file'
    assert 'GATK' in data.keys(), 'GATK not in input file'
    assert 'Savvy-GATK' in data.keys() or 'GATK-Savvy' in data.keys() , 'Savvy-GATK not in input file'
    assert 'Savvy-CNVKit' in data.keys() or 'CNVKit-Savvy' in data.keys(), 'Savvy-CNVKit not in input file'
    assert 'CNVKit-GATK' in data.keys() or 'GATK-CNVKit' in data.keys(), 'CNVKit-GATK not in input file'
    assert 'Triple' in data.keys(), 'Triple not in input file'
    return data

def main():
    args = get_args()
    for file,label,output in zip(args.inputs, args.labels, args.outputs):
        data = load_data(file)

        # if 'Savvy-GATK' is in keys sg = 'Savvy-GATK' else sg = 'GATK-Savvy' as ternal operator
        sg = 'Savvy-GATK' if 'Savvy-GATK' in data.keys() else 'GATK-Savvy'
        sc = 'Savvy-CNVKit' if 'Savvy-CNVKit' in data.keys() else 'CNVKit-Savvy'
        cg = 'CNVKit-GATK' if 'CNVKit-GATK' in data.keys() else 'GATK-CNVKit'
        
        key_order = ['Savvy', 'CNVKit', 'GATK', sg, sc, cg, 'Triple']
        groupings = [['Savvy'],['CNVkit'],['gCNV'],['Savvy','gCNV'],['Savvy','CNVkit'],['CNVkit','gCNV'],['Savvy','CNVkit','gCNV']]
        upset_data = [ data[k] for k in key_order]
        upset = from_memberships(groupings, upset_data)
        plt.rcParams.update({'font.size': 6})
        counts = not args.percent
        plot(upset, show_counts=counts)
        # change the font size of the text annotations
        plt.suptitle(label)
        print('saving to {}'.format(output))
        
        # force y axis to be log scale
        if args.logscale:
            plt.yscale('log')
        plt.savefig(output,dpi=300)
        plt.clf()

if __name__ == '__main__':
    main()
