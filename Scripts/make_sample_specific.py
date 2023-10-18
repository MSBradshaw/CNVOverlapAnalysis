import sys

# for line in stdin split on tabs and check 5 and 11 are the same string
seen = set()
for line in sys.stdin:
    line = line.strip()
    row = line.split('\t')
    call_1 = str(row[:6])
    if call_1 in seen:
        continue
    seen.add(call_1)
    if row[5] == row[11]:
        print(line)