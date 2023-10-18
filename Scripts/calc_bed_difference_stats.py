import argparse

def parse_args():
    # bed1, bed2
    parser = argparse.ArgumentParser()
    parser.add_argument('-a', type=str, help='bed 1')
    parser.add_argument('-b', type=str, help='bed 2')
    parser.add_argument('-l', type=str, help='label to prefix output with')
    
    args = parser.parse_args()
    return args

def count_lines(file):
    lines = 0
    for line in open(file,'r'):
        lines += 1
    return lines

def main():
    args = parse_args()
    b1 = count_lines(args.a)
    b2 = count_lines(args.b)
    a = max([b1,b2])
    b = min([b1,b2])
    print('{l}\t{b1}\t{b2}\t{diff}\t{percent}'.format(l=args.l, b1=a, b2=b, diff=b-a, percent=b/a))

if __name__ == "__main__":
    main()