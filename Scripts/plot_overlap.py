import pandas as pd
import matplotlib.pyplot as plt
import sys

"""
Inputs:
1. call_level_df - output from Scripts/compare_call_level_overlap.py
2. all_v_all_df - output from Scripts/compare_all_v_all_overlap.py
3. output_prefix - prefix for output file
4. overlap_type - 'reciprocal' or 'breakpoint'
"""

call_level_df = pd.read_csv(sys.argv[1], sep='\t', header=None)
all_v_all_df = pd.read_csv(sys.argv[2], sep='\t', header=None)
ol_type = sys.argv[4]

call_level_df.columns = ['count', 'percent', 'cnv_type']
all_v_all_df.columns = ['count', 'percent', 'cnv_type']

# 4 panel plot
# 0,0 - all vs all DEL
# 1,0 - call level DEL - box plot of values at each % on y axis
# 0,1 - all vs all DUP - line plot of percent on x axis, total number of overlaps on y axis
# 1,1 - call level DUP - box plot of values at each % on y axis

fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(15, 10))

# all vs all DEL
all_v_all_del = all_v_all_df[all_v_all_df['cnv_type'] == 'DEL']
axes[0,0].plot(all_v_all_del['percent'], all_v_all_del['count'], color='black')
axes[0,0].set_xlabel('Percent {} overlap'.format(ol_type))
axes[0,0].set_ylabel('Number of overlaps')
# remove top and right spines
axes[0,0].spines['top'].set_visible(False)
axes[0,0].spines['right'].set_visible(False)
axes[0,0].set_title('A.                                          All Deletions Overlaps',loc='left')
'        '

# call level DEL
call_level_del = call_level_df[call_level_df['cnv_type'] == 'DEL']
# create a list of lists, each inner list is the values for a given percent
del_values = []
percents = []
for percent in call_level_del['percent'].unique():
    percents.append(percent)
    del_values.append(call_level_del[call_level_del['percent'] == percent]['count'].values)
axes[1,0].boxplot(del_values, labels=percents)
axes[1,0].set_xlabel('Percent {} overlap'.format(ol_type))
axes[1,0].set_ylabel('Number of overlaps')
axes[1,0].spines['top'].set_visible(False)
axes[1,0].spines['right'].set_visible(False)
axes[1,0].set_title('C.                                  Call Level Deletion Overlaps',loc='left')

# all vs all DUP
all_v_all_dup = all_v_all_df[all_v_all_df['cnv_type'] == 'DUP']
axes[0,1].plot(all_v_all_dup['percent'], all_v_all_dup['count'], color='black')
axes[0,1].set_xlabel('Percent {} overlap'.format(ol_type))
axes[0,1].set_ylabel('Number of overlaps')
axes[0,1].spines['top'].set_visible(False)
axes[0,1].spines['right'].set_visible(False)
axes[0,1].set_title('B.                                          All Duplication Overlaps',loc='left')

# call level DUP
call_level_dup = call_level_df[call_level_df['cnv_type'] == 'DUP']
# create a list of lists, each inner list is the values for a given percent
dup_values = []
percents = []
for percent in call_level_dup['percent'].unique():
    percents.append(percent)
    dup_values.append(call_level_dup[call_level_dup['percent'] == percent]['count'].values)
axes[1,1].boxplot(dup_values, labels=percents)
axes[1,1].set_xlabel('Percent {} overlap'.format(ol_type))
axes[1,1].set_ylabel('Number of overlaps')
axes[1,1].spines['top'].set_visible(False)
axes[1,1].spines['right'].set_visible(False)
axes[1,1].set_title('D.                                  Call Level Duplication Overlaps',loc='left')

# make 0,0 and 0,1 share the same y scale and 1,0 and 1,1 share the same y scale
# find the largest y value in 0,0 and 0,1
y_max_0 = max(axes[0,0].get_ylim()[1], axes[0,1].get_ylim()[1])
# find the smallest y value in 0,0 and 0,1
y_min_0 = min(axes[0,0].get_ylim()[0], axes[0,1].get_ylim()[0])
# set the y limits for 0,0 and 0,1 to the same values
axes[0,0].set_ylim(y_min_0, y_max_0)
axes[0,1].set_ylim(y_min_0, y_max_0)

# find the largest y value in 1,0 and 1,1
y_max_1 = max(axes[1,0].get_ylim()[1], axes[1,1].get_ylim()[1])
# find the smallest y value in 1,0 and 1,1
y_min_1 = min(axes[1,0].get_ylim()[0], axes[1,1].get_ylim()[0])
# set the y limits for 1,0 and 1,1 to the same values
axes[1,0].set_ylim(y_min_1, y_max_1)
axes[1,1].set_ylim(y_min_1, y_max_1)

plt.tight_layout()
plt.savefig(sys.argv[3],dpi=300)
