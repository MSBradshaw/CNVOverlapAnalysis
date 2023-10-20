import sys

"""
Take as standard input a the first 12 columns of a bedtools -wao intersection, print out all lines that are not of the exact same call overlapping with itself.
"""

for line in sys.stdin:
    row = line.strip().split('\t')
    call1 = '\t'.join(row[:6])
    call2 = '\t'.join(row[6:])
    if call1 == call2:
        continue
    print(line.strip())