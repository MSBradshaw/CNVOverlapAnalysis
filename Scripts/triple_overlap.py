import argparse

# params s, c, g - all tsv files
# function to parse args

def parse_args():
    parser = argparse.ArgumentParser(description='Triple overlap')
    parser.add_argument('-s', type=str, help='savvy file')
    parser.add_argument('-c', type=str, help='cnvkit file')
    parser.add_argument('-g', type=str, help='gatk file')
    return parser.parse_args()

def get_overlap_lines_as_set(file, call):
    lines = set()
    for line in open(file,'r'):
        if call in line:
            lines.add(line.replace(call,'').strip())
    return lines

def get_all_call_pairs(file):
    pairs = {}
    for line in open(file,'r'):
        call1 = '-'.join(line.strip().split('\t')[:6])
        call2 = '-'.join(line.strip().split('\t')[6:-1])
        if call1 not in pairs:
            pairs[call1] = set()
        pairs[call1].add(call2)
        if call2 not in pairs:
            pairs[call2] = set()
        pairs[call2].add(call1)
    return pairs

def main():
    args = parse_args()
    all_triples = set()
    pairs = get_all_call_pairs(args.c)
    # print(len(pairs))
    # print(list(pairs.keys())[:10])
    i = 0
    savvy_hits = 0
    gatk_hits = 0
    for line in open(args.s,'r'):
        if 'GATK' not in line:
            continue
        if i > 100:
            break
        i += 1
        savvy_call = '-'.join(line.strip().split('\t')[:6])
        gatk_call = '-'.join(line.strip().split('\t')[6:-1])
        # get all the calls from the cnvkit file that overlap with this call
        # print('-', savvy_call, '-')
        # print('-', gatk_call, '-')
        if savvy_call in pairs:
            savvy_hits += 1
        else:
            continue
        if gatk_call in pairs:
            gatk_hits += 1
        else:
            continue
        savvy_cnvkit_ol = pairs[savvy_call]
        gatk_cnvkit_ol = pairs[gatk_call]
        # intersection of the two sets
        triple_ol_cnvkit = savvy_cnvkit_ol.intersection(gatk_cnvkit_ol)
        for call in triple_ol_cnvkit:
            three_calls = '\t'.join([savvy_call, gatk_call, call])
            all_triples.add(three_calls)
    # print('savvy hits', savvy_hits)
    # print('gatk hits', gatk_hits)
    # print so the output can be redirected to a file
    for t in all_triples:
        print(t.replace('-','\t'))

        
if __name__ == '__main__':
    main()
    