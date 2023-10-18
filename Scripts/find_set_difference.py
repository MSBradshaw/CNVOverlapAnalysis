import argparse

# params: -a bed file for series a, allows multiple files, -b bed file for series b, allows multiple files, -o output file name

def get_args():
    parser = argparse.ArgumentParser(description="finds the set difference between two sets of bed files")
    parser.add_argument("-a", "--a_bed", help="bed file for series a", action='append', required=True)
    parser.add_argument("-b", "--b_bed", help="bed file for series b", action='append', required=True)
    parser.add_argument("-o", "--output", help="output file name", required=True)
    return parser.parse_args()

def main():
    args = get_args()
    a_set = set()
    for file in args.a_bed:
        for line in open(file, 'r'):
            a_set.add(line.strip())
    b_set = set()
    for file in args.b_bed:
        for line in open(file, 'r'):
            b_set.add(line.strip())
    a_minus_b = a_set.difference(b_set)
    with open(args.output, 'w') as f:
        for line in a_minus_b:
            f.write(line + '\n')

if __name__ == "__main__":
    main()