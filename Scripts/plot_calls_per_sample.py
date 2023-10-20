import argparse
import pandas as pd
import matplotlib.pyplot as plt

def get_args():
    # --dels --dups --output
    parser = argparse.ArgumentParser()
    parser.add_argument('--dels', type=str, help='dels bed file')
    parser.add_argument('--dups', type=str, help='dups bed file')
    parser.add_argument('--output', type=str, help='output file')
    args = parser.parse_args()
    return args

def load_bed_file(file):
    df = pd.read_csv(file, sep='\t', header=None)
    df.columns = ['chr', 'start', 'end', 'call_type', 'sample', 'caller']
    return df

def main():
    args = get_args()
    fig, axes = plt.subplots(1,2, sharey=True)
    # grouping by caller and call_type plot the number of callers per sample as a violin plot
    dels = load_bed_file(args.dels)
    dups = load_bed_file(args.dups)
    x_labels = []
    for i,caller in enumerate(dels['caller'].unique()):
        x_labels.append(caller.replace('GATK Cohort Mode','gCNV'))
        sub = dels[dels['caller'] == caller]
        calls_per_samp = [sum(sub['sample'] == sample) for sample in sub['sample'].unique()]
        axes[0].violinplot(calls_per_samp,[i], showmeans=True)
    axes[0].set_xticks(range(len(x_labels)))
    axes[0].set_xticklabels(x_labels, rotation=90)
    axes[0].set_title('Deletions',loc='center')
    axes[0].set_title('A',loc='left')
    
    x_labels = []
    for i,caller in enumerate(dups['caller'].unique()):
        x_labels.append(caller.replace('GATK Cohort Mode','gCNV'))
        sub = dups[dups['caller'] == caller]
        calls_per_samp = [sum(sub['sample'] == sample) for sample in sub['sample'].unique()]
        axes[1].violinplot(calls_per_samp,[i], showmeans=True)
    axes[1].set_xticks(range(len(x_labels)))
    axes[1].set_xticklabels(x_labels, rotation=90)
    axes[1].set_title('Duplications',loc='center')
    axes[1].set_title('B',loc='left')

    # log scale the y axis
    axes[0].set_yscale('log')
    axes[1].set_yscale('log')

    # remove top and right borders
    axes[0].spines['top'].set_visible(False)
    axes[0].spines['right'].set_visible(False)
    axes[1].spines['top'].set_visible(False)
    axes[1].spines['right'].set_visible(False)

    plt.tight_layout()
    plt.savefig(args.output, dpi=300)


if __name__ == "__main__":
    main()
