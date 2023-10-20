import sys

# for line in stdin split on tabs and check 5 and 11 are the same string

for line in sys.stdin:
    line = line.strip()
    row = line.split('\t')
    # ensure the caller is not the same
    # callers should not be the same
    if row[5] == row[11]:
        continue
    # ensure the call is the same
    if row[4] == row[10]:
        print(line)