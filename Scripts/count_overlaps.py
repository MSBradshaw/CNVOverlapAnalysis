import argparse

"""
Params:
    -s - savvy file
    -c - cnvkit file
    -g - gatk file
    --sx - savvy-X overlap file
    --cx - cnvkit-X overlap file
    --gx - gatk-X overlap file
    -t - triple overlap file
    -o - output file
"""

def parse_args():
    parser = argparse.ArgumentParser(description='Triple overlap')
    parser.add_argument('-s', type=str, help='savvy file')
    parser.add_argument('-c', type=str, help='cnvkit file')
    parser.add_argument('-g', type=str, help='gatk file')
    parser.add_argument('--sx', type=str, help='savvy-x overlap file')
    parser.add_argument('--cx', type=str, help='cnvkit-x overlap file')
    parser.add_argument('--gx', type=str, help='gatk-x overlap file')
    parser.add_argument('-t', type=str, help='triple overlap file')
    parser.add_argument('-o', type=str, help='output file')
    # binary flag for p, if using bedpe files or not
    parser.add_argument('-p', default=False, help='use bedpe files', action='store_true')
    parser.add_argument('--fs', default=False, help='force calls to only count if they are from the same sample', action='store_true')
    return parser.parse_args()

def load_single_file(file):
    calls = set()
    for line in open(file,'r'):
        line = line.strip()
        line.replace('\t','-')
        calls.add(line)
    return calls

def def_load_doubles_as_set(file, alt_caller,use_bedpe=False, match_samples=False):
    if use_bedpe:
        idx = (None,13,None)
    else:
        idx = (None,6,-1)
    calls = set()
    for line in open(file,'r'):
        if alt_caller not in line:
            continue

        row1 = line.strip().split('\t')[idx[0]:idx[1]]
        row2 = line.strip().split('\t')[idx[1]:idx[2]]

        if match_samples and row1[-2] != row2[-2]: 
            continue
        
        call1 = '-'.join(row1)
        call2 = '-'.join(row2)

        calls.add(call1)
        calls.add(call2)
    return calls

def load_triple_file(file,use_bedpe=False, match_samples=False):
    if use_bedpe:
        idx = (None,13,26,None)
    else:
        idx = (None,6,12,None)
    calls = set()
    lines = 0
    for line in open(file,'r'):
        row1 = line.strip().split('\t')[idx[0]:idx[1]]
        row2 = line.strip().split('\t')[idx[1]:idx[2]]
        row3 = line.strip().split('\t')[idx[2]:idx[3]]
        if match_samples and (row1[-2] != row2[-2] or row1[-2] != row3[-2] or row2[-2] != row3[-2]):
            continue
        call1 = '-'.join(row1)
        call2 = '-'.join(row2)
        call3 = '-'.join(row3)
        calls.add(call1)
        calls.add(call2)
        calls.add(call3)
        lines += 1
    return calls, lines

def main():
    args = parse_args()
    savvy_calls = load_single_file(args.s)
    cnvkit_calls = load_single_file(args.c)
    gatk_calls = load_single_file(args.g)
    
    savvy_gatk = def_load_doubles_as_set(args.sx, 'GATK', args.p, args.fs)
    savvy_cnvkit = def_load_doubles_as_set(args.sx, 'CNVkit', args.p, args.fs)
    cnvkit_gatk = def_load_doubles_as_set(args.cx, 'GATK', args.p, args.fs)

    triple_calls, triple_lines = load_triple_file(args.t, args.p, args.fs)
    print('Triple calls: ', len(triple_calls))
    print('Triple lines: ', triple_lines)

    just_single_savvy = savvy_calls.difference(savvy_gatk).difference(savvy_cnvkit)
    just_single_cnvkit = cnvkit_calls.difference(cnvkit_gatk).difference(savvy_cnvkit)
    just_single_gatk = gatk_calls.difference(savvy_gatk).difference(cnvkit_gatk)

    just_single_savvy = just_single_savvy.difference(triple_calls)
    just_single_cnvkit = just_single_cnvkit.difference(triple_calls)
    just_single_gatk = just_single_gatk.difference(triple_calls)

    just_savvy_gatk = savvy_gatk.difference(savvy_cnvkit).difference(triple_calls)
    just_savvy_cnvkit = savvy_cnvkit.difference(savvy_gatk).difference(triple_calls)
    just_cnvkit_gatk = cnvkit_gatk.difference(savvy_gatk).difference(triple_calls)

    print('Savvy: ', len(just_single_savvy))
    print('CNVKit: ', len(just_single_cnvkit))
    print('GATK: ', len(just_single_gatk))
    print('Savvy-GATK: ', len(just_savvy_gatk))
    print('Savvy-CNVKit: ', len(just_savvy_cnvkit))
    print('CNVKit-GATK: ', len(just_cnvkit_gatk))
    print('Triple: ', len(triple_calls))

    # write those print statements to a file
    with open(args.o,'w') as f:
        f.write('Savvy: ' + str(len(just_single_savvy)) + '\n')
        f.write('CNVKit: ' + str(len(just_single_cnvkit)) + '\n')
        f.write('GATK: ' + str(len(just_single_gatk)) + '\n')
        f.write('Savvy-GATK: ' + str(len(just_savvy_gatk)) + '\n')
        f.write('Savvy-CNVKit: ' + str(len(just_savvy_cnvkit)) + '\n')
        f.write('CNVKit-GATK: ' + str(len(just_cnvkit_gatk)) + '\n')
        f.write('Triple: ' + str(triple_lines) + '\n')


if __name__ == '__main__':
    main()