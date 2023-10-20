import sys

# for line in stdin split on tabs and check 5 and 11 are the same string
seen = set()
for line in sys.stdin:
    line = line.strip()
    row = line.split('\t')
    if row[3] == row[9]:
        print(line)