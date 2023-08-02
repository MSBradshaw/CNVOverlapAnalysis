import pandas as pd
import matplotlib.pyplot as plt
import sys

"""
Inputs:
1. call_level_df - output from Scripts/compare_call_level_overlap.py
2. output_name - name for output file
3. overlap_type - 'reciprocal' or 'breakpoint'
"""

call_level_df = pd.read_csv(sys.argv[1], sep='\t', header=None)

ol_type = sys.argv[3]

call_level_df.columns = ['count', 'percent', 'cnv_type','slop_size']

call_level_df['percent'] = [ float('0.' + str(x).split('.')[-1]) for x in call_level_df['percent']]

fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(15, 10))

# 4 panel plot
# top # overlaps vs percent overlap
# bottom # overlaps vs slop size
# left DEL right DUP

# top - percent overlap vs count
del_sub = call_level_df[call_level_df['cnv_type'] == 'DEL']
del_values = []
percents = []
for percent in del_sub['percent'].unique():
    percents.append(percent)
    del_values.append(del_sub[del_sub['percent'] == percent]['count'].values)
axes[0,0].boxplot(del_values, labels=percents)
axes[0,0].set_title('DEL')
axes[0,0].set_xlabel('Slop as a percent of size')
axes[0,0].set_ylabel('# Overlaps')
# remove top and right spines
axes[0,0].spines['top'].set_visible(False)
axes[0,0].spines['right'].set_visible(False)

dup_sub = call_level_df[call_level_df['cnv_type'] == 'DUP']
dup_values = []
percents = []
for percent in del_sub['percent'].unique():
    percents.append(percent)
    dup_values.append(dup_sub[dup_sub['percent'] == percent]['count'].values)
axes[0,1].boxplot(dup_values, labels=percents)
axes[0,1].set_title('DUP')
axes[0,1].set_xlabel('Slop as a percent of size')
axes[0,1].set_ylabel('# Overlaps')
axes[0,1].spines['top'].set_visible(False)
axes[0,1].spines['right'].set_visible(False)

# bottom - actual slop size vs count
axes[1,0].scatter(del_sub['slop_size'], del_sub['count'],alpha=0.5, s=3)
axes[1,0].set_title('DEL')
axes[1,0].set_xlabel('Slop Size')
# log scale x axis
axes[1,0].set_xscale('log')
axes[1,0].set_ylabel('# Overlaps')
axes[1,0].spines['top'].set_visible(False)
axes[1,0].spines['right'].set_visible(False)

axes[1,1].scatter(dup_sub['slop_size'], dup_sub['count'],alpha=0.5, s=3)
axes[1,1].set_title('DUP')
axes[1,1].set_xlabel('Slop Size')
# log scale x axis
axes[1,1].set_xscale('log')
axes[1,1].set_ylabel('# Overlaps')
axes[1,1].spines['top'].set_visible(False)
axes[1,1].spines['right'].set_visible(False)

plt.tight_layout()
plt.savefig(sys.argv[2],dpi=300)