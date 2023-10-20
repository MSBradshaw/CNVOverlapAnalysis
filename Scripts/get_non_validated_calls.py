import argparse

def parse_args():
    # bed1, bed2
    parser = argparse.ArgumentParser()
    parser.add_argument('-a', type=str, help='bed 1')
    parser.add_argument('-b', type=str, help='bed 2')
    parser.add_argument('-o', type=str, help='output file')
    
    args = parser.parse_args()
    return args

def get_lines(file):
    lines = set()
    for line in open(file,'r'):
        lines.add(line.strip())
    return lines

def main():
    args = parse_args()
    b1 = get_lines(args.a)
    b2 = get_lines(args.b)
    # get things from b1 not in b2
    b1 = b1.difference(b2)
    with open(args.o,'w') as out:
        for line in b1:
            out.write(line + '\n')
    

if __name__ == "__main__":
    main()