import sys

for line in sys.stdin:
    row = line.strip().split('\t')
    print('\t'.join(row[:6]))
    print('\t'.join(row[6:]))