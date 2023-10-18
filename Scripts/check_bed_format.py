import sys

for line in open(sys.argv[1],'r'):
    row = line.strip().split('\t')
    if int(row[2]) < int(row[1]):
        print(line)
        sys.exit(1)