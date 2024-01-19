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
    fig, axes = plt.subplots(2,4)
    fig.set_size_inches(16,8)
    """
    plot format
    DELS: per sample, gCNV sizes, Savvy sizes, CNVkit sizes
    DUPS: per sample, gCNV sizes, Savvy sizes, CNVkit sizes
    per sample (s) should share y, sizes should share x and y
    """
    callers = ['gCNV', 'Savvy', 'CNVkit']
    # grouping by caller and call_type plot the number of callers per sample as a violin plot
    dels = load_bed_file(args.dels)
    dups = load_bed_file(args.dups)

    # calc the average number of calls per sample
    print('Average number of calls per sample')
    for caller in callers:
        sub = dels[dels['caller'] == caller]
        print('Deletion', caller, 'average=',sub.shape[0]/len(sub['sample'].unique()))
    for caller in callers:
        sub = dups[dups['caller'] == caller]
        print('Duplication', caller, 'average=',sub.shape[0]/len(sub['sample'].unique()))

    caller2color = {'gCNV':'#2873ae', 'CNVkit':'#2c982c', 'Savvy':'#ff7614'} # blue, green, orange

    x_labels = []
    for i,caller in enumerate(callers):
        x_labels.append(caller.replace('GATK Cohort Mode','gCNV'))
        sub = dels[dels['caller'] == caller]
        calls_per_samp = [sum(sub['sample'] == sample) for sample in sub['sample'].unique()]
        axes[0,0].violinplot(calls_per_samp,[i], showmeans=True)
    axes[0,0].set_xticks(range(len(x_labels)))
    axes[0,0].set_xticklabels(x_labels, rotation=90)
    axes[0,0].set_title('Deletions',loc='center')
    
    x_labels = []
    for i,caller in enumerate(callers):
        x_labels.append(caller.replace('GATK Cohort Mode','gCNV'))
        sub = dups[dups['caller'] == caller]
        calls_per_samp = [sum(sub['sample'] == sample) for sample in sub['sample'].unique()]
        axes[1,0].violinplot(calls_per_samp,[i], showmeans=True)
    
    axes[1,0].set_xticks(range(len(x_labels)))
    axes[1,0].set_xticklabels(x_labels, rotation=90)
    axes[1,0].set_title('Duplications',loc='center')
    
    axes[0,0].set_ylabel('Number of calls')
    axes[1,0].set_ylabel('Number of calls')

    # log scale the y axis
    axes[0,0].set_yscale('log')
    axes[1,0].set_yscale('log')

    # match y scale of 0,0 and 1,0
    # find which is bigger
    max_y = max(axes[0,0].get_ylim()[1], axes[1,0].get_ylim()[1])
    axes[0,0].set_ylim([1,max_y])

    # remove top and right borders
    axes[0,0].spines['top'].set_visible(False)
    axes[0,0].spines['right'].set_visible(False)
    axes[1,0].spines['top'].set_visible(False)
    axes[1,0].spines['right'].set_visible(False)

    
    # add a size column to the dataframes
    dels['size'] = dels['end'] - dels['start']
    dups['size'] = dups['end'] - dups['start']
    # plot the size distribution of each caller
    # dels
    max_y = 0
    max_x = 0
    for i,caller in enumerate(callers):
        sub = dels[dels['caller'] == caller]
        axes[0,i+1].hist(sub['size'], bins=100, color=caller2color[caller])
        axes[0,i+1].set_title(caller + ' DEL', loc='center')
        sub = dups[dups['caller'] == caller]
        axes[1,i+1].hist(sub['size'], bins=100, color=caller2color[caller])
        axes[1,i+1].set_title(caller + ' DUP', loc='center')
        # turn of top and right spines
        axes[0,i+1].spines['top'].set_visible(False)
        axes[0,i+1].spines['right'].set_visible(False)
        axes[1,i+1].spines['top'].set_visible(False)
        axes[1,i+1].spines['right'].set_visible(False)
        # track max y and x
        max_y = max(max_y, axes[0,i+1].get_ylim()[1], axes[1,i+1].get_ylim()[1])
        max_x = max(max_x, axes[0,i+1].get_xlim()[1], axes[1,i+1].get_xlim()[1])
        # log y
        axes[0,i+1].set_yscale('log')
        axes[1,i+1].set_yscale('log')
        # label x and y
        axes[0,i+1].set_xlabel('Size (bp)')
        axes[0,i+1].set_ylabel('Number of calls')
        axes[1,i+1].set_xlabel('Size (bp)')
        axes[1,i+1].set_ylabel('Number of calls')

    plt.tight_layout()
    plt.savefig(args.output, dpi=300)


if __name__ == "__main__":
    main()
