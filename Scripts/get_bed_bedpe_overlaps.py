import argparse

"""
Example command:
        python Scripts/get_bed_bedpe_overlaps.py \
            --sb work/CallerSpecificOverlapsBed/venn.{wildcards.cnv_type}.Savvy.bed.txt
            --sp work/CallerSpecificOverlapsBedPE/venn.{wildcards.cnv_type}.Savvy.bedpe.txt \
            --cb work/CallerSpecificOverlapsBed/venn.{wildcards.cnv_type}.CNVkit.bed.txt \
            --cp  work/CallerSpecificOverlapsBedPE/venn.{wildcards.cnv_type}.CNVkit.bedpe.txt \
            --gb work/CallerSpecificOverlapsBed/venn.{wildcards.cnv_type}.GATK.bed.txt \
            --gp  work/CallerSpecificOverlapsBedPE/venn.{wildcards.cnv_type}.GATK.bedpe.txt \
            --tb {input.triple} \
            --tp {input.bedpe_triple} > {output}
"""

# get args function
def parse_args():
    parser = argparse.ArgumentParser(description='Get overlaps between bed and bedpe files')
    parser.add_argument('--sb', help='Savvy bed file')
    parser.add_argument('--sp', help='Savvy bedpe file')
    parser.add_argument('--cb', help='CNVkit bed file')
    parser.add_argument('--cp', help='CNVkit bedpe file')
    parser.add_argument('--gb', help='GATK bed file')
    parser.add_argument('--gp', help='GATK bedpe file')
    parser.add_argument('--tb', help='Triple bed file')
    parser.add_argument('--tp', help='Triple bedpe file')
    parser.add_argument('--sx', help='temporary file for Savvy bed and bedpe overlaps')
    parser.add_argument('--cx', help='temporary file for CNVkit bed and bedpe overlaps')
    parser.add_argument('--gx', help='temporary file for GATK bed and bedpe overlaps')
    parser.add_argument('--tx', help='temporary file for triple bed and bedpe overlaps')

    return parser.parse_args()

def load_lines_as_set(f):
    calls = set()
    for line in open(f,'r'):
        line = line.replace('GATK Cohort Mode','GATK')
        row = line.strip().split('\t')
        line = '-'.join(row)
        calls.add(line)
    return calls

def check_double_overlaps(bedfile,bedpefile):
    """
    Check which calls overlap based on recirpocal overlap and breakpoint overlap criteria
    Return a set of calls in bed format that meet both criteria
    """
    bedpe_calls = load_lines_as_set(bedpefile)
    template = '{chrom}-{sm1}-{sp1}-{chrom}-{em1}-{ep1}-.-.-.-.-{cnv_type}-{samp}-{caller}'
    matches = set()
    for i,line in enumerate(open(bedfile,'r')):
        line = line.replace('GATK Cohort Mode','GATK')
        row = line.strip().split('\t')
        chrom1 = row[0]
        s1 = int(row[1])
        e1 = int(row[2])
        cnv_t1 = row[3]
        samp1 = row[4]
        caller1 = row[5]

        chrom2 = row[6]
        s2 = int(row[7])
        e2 = int(row[8])
        cnv_t2 = row[9]
        samp2 = row[10]
        caller2 = row[11]

        call1 = template.format(chrom=chrom1,sm1=str(s1-1),sp1=str(s1+1),
                                em1=str(e1-1),ep1=str(e1+1),cnv_type=cnv_t1, samp=samp1, caller=caller1)
        call2 = template.format(chrom=chrom2,sm1=str(s2-1),sp1=str(s2+1),
                                em1=str(e2-1),ep1=str(e2+1),cnv_type=cnv_t2, samp=samp2, caller=caller2)
        call = call1 + '-' + call2
        if call in bedpe_calls:
            matches.add(line.strip())
    return matches

def check_triple_overlaps(bedfile,bedpefile):
    bedpe_calls = load_lines_as_set(bedpefile)
    template = '{chrom}-{sm1}-{sp1}-{chrom}-{em1}-{ep1}-.-.-.-.-{cnv_type}-{samp}-{caller}'
    matches = set()
    for line in open(bedfile,'r'):
        line = line.replace('GATK Cohort Mode','GATK')
        row = line.strip().split('\t')
        chrom1 = row[0]
        s1 = int(row[1])
        e1 = int(row[2])
        cnv_t1 = row[3]
        samp1 = row[4]
        caller1 = row[5]

        chrom2 = row[6]
        s2 = int(row[7])
        e2 = int(row[8])
        cnv_t2 = row[9]
        samp2 = row[10]
        caller2 = row[11]

        chrom3 = row[12]
        s3 = int(row[13])
        e3 = int(row[14])
        cnv_t3 = row[15]
        samp3 = row[16]
        caller3 = row[17]

        call1 = template.format(chrom=chrom1,sm1=str(s1-1),sp1=str(s1+1),
                                em1=str(e1-1),ep1=str(e1+1),cnv_type=cnv_t1, samp=samp1, caller=caller1)
        call2 = template.format(chrom=chrom2,sm1=str(s2-1),sp1=str(s2+1),
                                em1=str(e2-1),ep1=str(e2+1),cnv_type=cnv_t2, samp=samp2, caller=caller2)
        call3 = template.format(chrom=chrom3,sm1=str(s3-1),sp1=str(s3+1),
                                em1=str(e3-1),ep1=str(e3+1),cnv_type=cnv_t3, samp=samp3, caller=caller3)
        combos = set()
        for c1 in [call1,call2,call3]:
            for c2 in [call1,call2,call3]:
                if c1 == c2:
                    continue
                for c3 in [call1,call2,call3]:
                    if c1 == c3 or c2 == c3:
                        continue
                    combos.add(c1 + '-' + c2 + '-' + c3)
        for call in combos:
            if call in bedpe_calls:
                matches.add(line.strip())
    return matches

def write_to_file(calls,file):
    with open(file,'w') as out:
        for call in calls:
            call = call.replace('-','\t')
            out.write(call + '\n')    

def main():
    args = parse_args()
    # check doubles for savvy bed and bedpe
    savvy_doubles = check_double_overlaps(args.sb,args.sp)
    # check doubles for cnvkit bed and bedpe
    cnvkit_doubles = check_double_overlaps(args.cb,args.cp)
    # check doubles for gatk bed and bedpe
    gatk_doubles = check_double_overlaps(args.gb,args.gp)
    # check triples bed and bedpe
    triples = check_triple_overlaps(args.tb,args.tp)
    # write to file
    write_to_file(savvy_doubles, args.sx)
    write_to_file(cnvkit_doubles, args.cx)
    write_to_file(gatk_doubles, args.gx)
    write_to_file(triples, args.tx)

if __name__ == '__main__':
    main()