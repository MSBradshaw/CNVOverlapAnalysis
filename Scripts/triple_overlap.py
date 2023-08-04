import argparse

# params s, c, g - all tsv files
# function to parse args

def parse_args():
    parser = argparse.ArgumentParser(description='Triple overlap')
    parser.add_argument('-s', type=str, help='savvy file')
    parser.add_argument('-c', type=str, help='cnvkit file')
    parser.add_argument('-g', type=str, help='gatk file')
    parser.add_argument('-p', default=False, help='use bedpe files', action='store_true')
    return parser.parse_args()

def get_overlap_lines_as_set(file, call):
    lines = set()
    for line in open(file,'r'):
        if call in line:
            lines.add(line.replace(call,'').strip())
    return lines

def get_all_call_pairs(file,use_bedpe=False):
    if use_bedpe:
        idx = (None,13,None)
    else:
        idx = (None,6,-1)
    pairs = {}
    for line in open(file,'r'):
        call1 = '-'.join(line.strip().split('\t')[idx[0]:idx[1]])
        call2 = '-'.join(line.strip().split('\t')[idx[1]:idx[2]])
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
    pairs = get_all_call_pairs(args.c,args.p)
    # print('pairs',len(pairs))
    # print(list(pairs.keys())[len(pairs)-10:])
    if args.p:
        idx = (None,13,None)
    else:
        idx = (None,6,-1)
    i = 0
    savvy_hits = 0
    gatk_hits = 0
    for line in open(args.s,'r'):
        if 'GATK' not in line:
            continue
        # if i > 100:
        #     break
        # i += 1
        savvy_call = '-'.join(line.strip().split('\t')[idx[0]:idx[1]])
        gatk_call = '-'.join(line.strip().split('\t')[idx[1]:idx[2]])
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
    # print so the output can be redirected to a file
    for t in all_triples:
        print(t.replace('-','\t'))


if __name__ == '__main__':
    main()
    