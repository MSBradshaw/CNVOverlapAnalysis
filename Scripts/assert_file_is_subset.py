import sys

# load the files in sys.argv[1] and sys.argv[2] into memory as sets

a = set([ line.strip() for line in open(sys.argv[1],'r')])
b = set([ line.strip() for line in open(sys.argv[2],'r')])

for x in a:
    assert x in b, x + ' not in ' + sys.argv[2]