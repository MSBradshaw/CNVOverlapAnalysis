import argparse

def get_args():
    """
    -a bed
    -b bed
    -c bed
    -o output
    """

    parser = argparse.ArgumentParser(description='Get triples from three bed files')
    parser.add_argument('-a', '--a_bed', type=str, required=True, help='First bed.txt file')
    parser.add_argument('-b', '--b_bed', type=str, required=True, help='Second bed.txt file')
    parser.add_argument('-c', '--c_bed', type=str, required=True, help='Third bed.txt file')
    parser.add_argument('-o', '--output', type=str, required=True, help='Output file')
    args = parser.parse_args()

    return args

def load_file_as_dict(file):
    calls_dict = {}
    for line in open(file,'r'):
        row = line.strip().split('\t')
        call1 = '\t'.join(row[:6])
        call2 = '\t'.join(row[6:])
        if call1 not in calls_dict:
            calls_dict[call1] = set()
        if call2 not in calls_dict:
            calls_dict[call2] = set()
        calls_dict[call1].add(call2)
        calls_dict[call2].add(call1)
    return calls_dict


def main():
    args = get_args()
    a_dict = load_file_as_dict(args.a_bed)
    b_dict = load_file_as_dict(args.b_bed)
    c_dict = load_file_as_dict(args.c_bed)
    # combine the three dicts into one
    combined_dict = {}
    for key in a_dict:
        combined_dict[key] = a_dict[key]
    for key in b_dict:
        if key in combined_dict:
            combined_dict[key] = combined_dict[key].union(b_dict[key])
        else:
            combined_dict[key] = b_dict[key]
    for key in c_dict:
        if key in combined_dict:
            combined_dict[key] = combined_dict[key].union(c_dict[key])
        else:
            combined_dict[key] = c_dict[key]
    
    triples = []
    for key_a in a_dict.keys():
        for key_b in b_dict.keys():
            for key_c in c_dict.keys():
                if key_b in combined_dict[key_a] and key_c in combined_dict[key_a] and key_c in combined_dict[key_b]:
                    tmp_triple = [key_b, key_c, key_a]
                    tmp_triple.sort()
                    triples.append(tmp_triple)
    # output triples
    with open(args.output, 'w') as f:
        for triple in triples:
            f.write('\n'.join(triple) + '\n')


                    




if __name__ == '__main__':
    main()