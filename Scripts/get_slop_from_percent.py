import sys

bedpe = sys.argv[1]
portion = float(sys.argv[2])

# 1	35255999	35256001	hs37d5	35256999	35257001	.	.	DUP	HG01378	Savvy	
row = None
with open(bedpe, 'r') as f:
    row = f.readline().strip().split('\t')

chr1, strart1, end1, chr2, start2, end2, skip1, skip2, skip3, skip4, call_type, sample_id, caller = row

size_of_call = int(start2) - int(end1)

print(int(size_of_call * portion))