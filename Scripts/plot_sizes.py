import argparse
import matplotlib.pyplot as plt

# get the arguments function

def get_args():
    # inputs -i allow multiple
    parser = argparse.ArgumentParser(description='Plot the sizes of calls from bed file(s)')
    parser.add_argument('-i', '--input', action='append', required=True, help='Input bed file(s)')
    # labels from the input groups
    parser.add_argument('-l', '--labels', action='append', required=True, help='Labels for input file(s)')
    parser.add_argument('-o', '--output', required=True, help='Output file')
    args = parser.parse_args()
    return args

# function to get sizes of calls from bed file
def get_sizes(bed_file):
    sizes_dels = []
    sizes_dups = []
    with open(bed_file, 'r') as bed:
        for line in bed:
            line = line.strip().split()
            if 'DEL' in line:
                sizes_dels.append(int(line[2]) - int(line[1]))
            elif 'DUP' in line:
                sizes_dups.append(int(line[2]) - int(line[1]))
            else:
                print('ERROR: line not a DEL or DUP')
                print(line)
    return sizes_dels, sizes_dups

def main():
    args = get_args()
    assert len(args.input) == len(args.labels), 'Number of labels must match number of input files'
    num_inputs = len(args.input)
    # give me a list of 5 color blind friendly colors as hex codes
    colors = ['#0173b2', '#de8f05', '#029e73', '#d55e00', '#cc78bc']
    fig, ax = plt.subplots(2,num_inputs, figsize= (5 * num_inputs,5))
    print(args.input)
    print(args.labels)
    for i,bed_file,label in zip(list(range(num_inputs)),args.input, args.labels):
        del_sizes, dup_sizes = get_sizes(bed_file)
        ax[0,i].hist(del_sizes, bins=100, color=colors[i])
        ax[1,i].hist(dup_sizes, bins=100, color=colors[i])
        ax[0,i].set_title(chr(i + 65), loc='left')
        ax[0,i].set_title(label + ' DEL\nn={}'.format(str(len(del_sizes))), loc='center')
        ax[1,i].set_title(chr(i + 65 + num_inputs), loc='left')
        ax[1,i].set_title( label +  ' DUP\nn={}'.format(str(len(dup_sizes))), loc='center')
        ax[0,i].set_xlabel('Size (bp)')
        ax[1,i].set_xlabel('Size (bp)')
        ax[0,i].set_ylabel('Number of calls')
        ax[1,i].set_ylabel('Number of calls')
        # remove top and right borders
        ax[0,i].spines['top'].set_visible(False)
        ax[0,i].spines['right'].set_visible(False)
        ax[1,i].spines['top'].set_visible(False)
        ax[1,i].spines['right'].set_visible(False)
        # log x and y
        ax[0,i].set_yscale('log')
        ax[1,i].set_yscale('log')
    plt.tight_layout()
    plt.savefig(args.output, dpi=300)

        




if __name__ == '__main__':
    main()
    


