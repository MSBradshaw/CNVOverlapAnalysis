import argparse
import pandas as pd
import matplotlib.pyplot as plt

# inputs 1. slop df 2. recip df 3 output name
def get_args():
    parser = argparse.ArgumentParser(description='Plot overlap results')
    parser.add_argument('-s', '--slop', action='store', dest='slop', help='Slop overlap results')
    parser.add_argument('-r', '--recip', action='store', dest='recip', help='Reciprocal overlap results')
    parser.add_argument('-o', '--output', action='store', dest='output_name', help='Output name')
    args = parser.parse_args()
    return args

def main():
    args = get_args()
    # load slop df
    slop_df = pd.read_csv(args.slop, sep='\t', header=None)
    slop_df.columns = ['count', 'percent', 'cnv_type','slop_size']
    slop_df['percent'] = [ float('0.' + str(x).split('.')[-1]) for x in slop_df['percent']]
    # load recip df
    recip_df = pd.read_csv(args.recip, sep='\t', header=None)
    recip_df.columns = ['count', 'percent', 'cnv_type']

    fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(15, 10))
    
    # 0,0 Recip DEL
    call_level_del = recip_df[recip_df['cnv_type'] == 'DEL']
    # create a list of lists, each inner list is the values for a given percent
    del_values = []
    percents = []
    for percent in call_level_del['percent'].unique():
        percents.append(percent)
        del_values.append(call_level_del[call_level_del['percent'] == percent]['count'].values)
    axes[0,0].boxplot(del_values, labels=percents)
    axes[0,0].set_xlabel('Percent overlap')
    axes[0,0].set_ylabel('Number of overlaps')
    axes[0,0].spines['top'].set_visible(False)
    axes[0,0].spines['right'].set_visible(False)
    axes[0,0].set_title('A.                                  Reciprocal Overlap Deletions',loc='left')

    # 0,1 Recip DUP
    call_level_dup = recip_df[recip_df['cnv_type'] == 'DUP']
    # create a list of lists, each inner list is the values for a given percent
    del_values = []
    percents = []
    for percent in call_level_dup['percent'].unique():
        percents.append(percent)
        del_values.append(call_level_dup[call_level_dup['percent'] == percent]['count'].values)
    axes[0,1].boxplot(del_values, labels=percents)
    axes[0,1].set_xlabel('Percent overlap')
    axes[0,1].set_ylabel('Number of overlaps')
    axes[0,1].spines['top'].set_visible(False)
    axes[0,1].spines['right'].set_visible(False)
    axes[0,1].set_title('B.                                  Reciprocal Overlap Duplications',loc='left')

    # slop plots
    # top - percent overlap vs count
    del_sub = slop_df[slop_df['cnv_type'] == 'DEL']
    del_values = []
    percents = []
    for percent in del_sub['percent'].unique():
        percents.append(percent)
        del_values.append(del_sub[del_sub['percent'] == percent]['count'].values)
    axes[1,0].boxplot(del_values, labels=percents)
    axes[1,0].set_xlabel('Slop as a percent of size')
    axes[1,0].set_ylabel('# Overlaps')
    axes[1,0].spines['top'].set_visible(False)
    axes[1,0].spines['right'].set_visible(False)
    axes[1,0].set_title('C.                                  Breakpoint Overlap Deletions',loc='left')

    dup_sub = slop_df[slop_df['cnv_type'] == 'DUP']
    dup_values = []
    percents = []
    for percent in del_sub['percent'].unique():
        percents.append(percent)
        dup_values.append(dup_sub[dup_sub['percent'] == percent]['count'].values)
    axes[1,1].boxplot(dup_values, labels=percents)
    axes[1,1].set_xlabel('Slop as a percent of size')
    axes[1,1].set_ylabel('# Overlaps')
    axes[1,1].spines['top'].set_visible(False)
    axes[1,1].spines['right'].set_visible(False)
    axes[1,1].set_title('D.                                  Breakpoint Overlap Duplications',loc='left')

    # save figure
    # use tight layout to make sure all text is visible
    plt.tight_layout()
    plt.savefig(args.output_name,dpi=300)



if __name__ == '__main__':
    main()