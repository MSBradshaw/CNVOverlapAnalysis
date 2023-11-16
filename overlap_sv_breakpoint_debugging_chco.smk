cnv_calls = 'Data/all_calls.bed'

def convert_to_bedpe(infile,outfile):
    t = '{chrom}\t{s1}\t{e1}\t{chrom}\t{s1}\t{e1}\t.\t.\t.\t.\t{cnv_type}\t{sample}\t{caller}\n'
    with open(outfile,'w') as out:
        for line in open(infile,'r'):
            row = line.strip().split('\t')
            chrom = row[0]
            start = int(row[1])
            end = int(row[2])
            cnv_type = row[3]
            sample = row[4]
            caller = row[5]
            s1 = start - 1
            e1 = start + 1
            s2 = end - 1
            e2 = end + 1
            # **locals returns a dictionary of all local variables
            out.write(t.format(**locals()))

DIRS = ['CNVs']
CNV_TYPES = ['DEL','DUP']

rule all:
    input:
        'breakpointworkchco/triple_overlaps.DEL.bed',
        'breakpointworkchco/triple_overlaps.DUP.bed'
    
    

rule sort:
    input:
        cnv_calls = cnv_calls
    output:
        cnv_calls_sorted = 'breakpointworkchco/CNVs/all_calls.sorted.bed'
    shell:
        """
        mkdir -p breakpointworkchco
        mkdir -p breakpointworkchco/CNVs
        bedtools sort -i {input.cnv_calls} | sed -e 's/GATK Cohort Mode/GATK/g' > {output.cnv_calls_sorted}
        """

rule convert_to_bedpe:
    input:
        'breakpointworkchco/CNVs/all_calls.sorted.bed'
    output:
        'breakpointworkchco/CNVs/all_calls.bedpe.tmp'
    run:
        convert_to_bedpe(input[0],output[0])

rule clean_bedpe:
    input:
        'breakpointworkchco/CNVs/all_calls.bedpe.tmp'
    output:
        'breakpointworkchco/CNVs/all_calls.bedpe'
    shell:
        """
        # remove all lines that have a -1 in them as they just throw thousands of warnings
        grep -v "\-1" {input} > {output}
        """

rule intersect_bedpe:
    input:
        bedpe='breakpointworkchco/CNVs/all_calls.bedpe'
    output:
        bedpe_intersect_DEL='breakpointworkchco/CNVs/all_calls.intersect.DEL.bedpe',
        bedpe_intersect_DUP='breakpointworkchco/CNVs/all_calls.intersect.DUP.bedpe'
    shell:
        """
        bash Scripts/pairtopair_slop.sh {input.bedpe} {input.bedpe} 0.50 breakpointworkchco/CNVs/all_calls.intersect.bedpe
        cat breakpointworkchco/CNVs/all_calls.intersect.bedpe | grep DEL | grep -v DUP > {output.bedpe_intersect_DEL}.tmp
        cat breakpointworkchco/CNVs/all_calls.intersect.bedpe | grep DUP | grep -v DEL > {output.bedpe_intersect_DUP}.tmp
        
        # remove calls where the cnv_type is not the same

        # remove intersections where the sample is not the same
        # remove rows where 12 == 25 in awk
        cat breakpointworkchco/CNVs/all_calls.intersect.DEL.bedpe.tmp | awk '$12 != $25' > breakpointworkchco/CNVs/all_calls.intersect.DEL.filtered.bedpe.tmp2
        cat breakpointworkchco/CNVs/all_calls.intersect.DUP.bedpe.tmp | awk '$12 != $25' > breakpointworkchco/CNVs/all_calls.intersect.DUP.filtered.bedpe.tmp2
        
        # remove intersections where the caller is the same
        cat breakpointworkchco/CNVs/all_calls.intersect.DUP.bedpe.tmp | awk '$13 != $26' > breakpointworkchco/CNVs/all_calls.intersect.DUP.filtered.bedpe.tmp2
        cat breakpointworkchco/CNVs/all_calls.intersect.DEL.bedpe.tmp | awk '$13 != $26' > breakpointworkchco/CNVs/all_calls.intersect.DEL.filtered.bedpe.tmp2

        # split the lines in the bedpe so there is one call in each line
        # awk columns 1-13] and columns [14-26] with a new line between them
        cat breakpointworkchco/CNVs/all_calls.intersect.DEL.filtered.bedpe.tmp2 | awk -v OFS='\t' '{{print $1"\\t"$2"\\t"$3"\\t"$4"\\t"$5"\\t"$6"\\t"$7"\\t"$8"\\t"$9"\\t"$10"\\t"$11"\\t"$12"\\t"$13"\\n"$14"\\t"$15"\\t"$16"\\t"$17"\\t"$18"\\t"$19"\\t"$20"\\t"$21"\\t"$22"\\t"$23"\\t"$24"\\t"$25"\\t"$26}}' > {output.bedpe_intersect_DEL}
        cat breakpointworkchco/CNVs/all_calls.intersect.DUP.filtered.bedpe.tmp2 | awk -v OFS='\t' '{{print $1"\\t"$2"\\t"$3"\\t"$4"\\t"$5"\\t"$6"\\t"$7"\\t"$8"\\t"$9"\\t"$10"\\t"$11"\\t"$12"\\t"$13"\\n"$14"\\t"$15"\\t"$16"\\t"$17"\\t"$18"\\t"$19"\\t"$20"\\t"$21"\\t"$22"\\t"$23"\\t"$24"\\t"$25"\\t"$26}}' > {output.bedpe_intersect_DUP}
        """

def convert_bedpe_to_bed(infile,outfile):
    t = '{chrom}\t{start}\t{end}\t{cnv_type}\t{sample}\t{caller}\n'
    with open(outfile,'w') as out:
        for line in open(infile,'r'):
            row = line.strip().split('\t')
            chrom = row[0]
            start = int(row[1]) + 1
            end = int(row[4]) + 1
            cnv_type = row[10]
            sample = row[11]
            caller = row[12]
            # **locals returns a dictionary of all local variables
            out.write(t.format(**locals()))

rule convert_bedpe_to_bed:
    input:
        bedpe_intersect_DEL='breakpointworkchco/CNVs/all_calls.intersect.DEL.bedpe',
        bedpe_intersect_DUP='breakpointworkchco/CNVs/all_calls.intersect.DUP.bedpe'
    output:
        bed_intersect_DEL='breakpointworkchco/CNVs/all_calls.intersect.DEL.bed',
        bed_intersect_DUP='breakpointworkchco/CNVs/all_calls.intersect.DUP.bed'
    run:
        convert_bedpe_to_bed(input[0],output[0])
        convert_bedpe_to_bed(input[1],output[1])


def get_calls_from_bed(filename):
    cnvs = set()
    for line in open(filename):
        cnvs.add(line.strip())
    return cnvs

rule triple_overlap:
    input:
        cnv_calls_sorted = 'breakpointworkchco/CNVs/all_calls.bedpe'

    output:
        'breakpointworkchco/triple_overlaps.DEL.bed',
        'breakpointworkchco/triple_overlaps.DUP.bed'
    shell:
        """
            grep GATK breakpointworkchco/CNVs/all_calls.bedpe > breakpointworkchco/all_calls.GATK.sorted.bedpe
            grep CNVkit breakpointworkchco/CNVs/all_calls.bedpe > breakpointworkchco/all_calls.CNVkit.sorted.bedpe
            grep Savvy breakpointworkchco/CNVs/all_calls.bedpe > breakpointworkchco/all_calls.Savvy.sorted.bedpe

            grep DEL breakpointworkchco/all_calls.GATK.sorted.bedpe > breakpointworkchco/all_calls.GATK.DEL.sorted.bedpe
            grep DUP breakpointworkchco/all_calls.GATK.sorted.bedpe > breakpointworkchco/all_calls.GATK.DUP.sorted.bedpe
            grep DEL breakpointworkchco/all_calls.CNVkit.sorted.bedpe > breakpointworkchco/all_calls.CNVkit.DEL.sorted.bedpe
            grep DUP breakpointworkchco/all_calls.CNVkit.sorted.bedpe > breakpointworkchco/all_calls.CNVkit.DUP.sorted.bedpe
            grep DEL breakpointworkchco/all_calls.Savvy.sorted.bedpe > breakpointworkchco/all_calls.Savvy.DEL.sorted.bedpe
            grep DUP breakpointworkchco/all_calls.Savvy.sorted.bedpe > breakpointworkchco/all_calls.Savvy.DUP.sorted.bedpe


            bash Scripts/pairtopair_slop.sh breakpointworkchco/all_calls.GATK.DEL.sorted.bedpe breakpointworkchco/all_calls.Savvy.DEL.sorted.bedpe 0.50 breakpointworkchco/all_calls.GATK.DEL.x_Savvy.DEL.bedpe.tmp 
            cat breakpointworkchco/all_calls.GATK.DEL.x_Savvy.DEL.bedpe.tmp | awk '$12 != $25' > breakpointworkchco/all_calls.GATK.DEL.x_Savvy.DEL.bedpe

            bash Scripts/pairtopair_slop.sh breakpointworkchco/all_calls.GATK.DEL.sorted.bedpe breakpointworkchco/all_calls.CNVkit.DEL.sorted.bedpe  0.50 breakpointworkchco/all_calls.GATK.DEL.x_CNVkit.DEL.bedpe.tmp 
            cat breakpointworkchco/all_calls.GATK.DEL.x_CNVkit.DEL.bedpe.tmp | awk '$12 != $25' > breakpointworkchco/all_calls.GATK.DEL.x_CNVkit.DEL.bedpe

            bash Scripts/pairtopair_slop.sh breakpointworkchco/all_calls.Savvy.DEL.sorted.bedpe breakpointworkchco/all_calls.CNVkit.DEL.sorted.bedpe 0.50 breakpointworkchco/all_calls.Savvy.DEL.x_CNVkit.DEL.bedpe.tmp 
            cat breakpointworkchco/all_calls.Savvy.DEL.x_CNVkit.DEL.bedpe.tmp | awk '$12 != $25' > breakpointworkchco/all_calls.Savvy.DEL.x_CNVkit.DEL.bedpe

            
            bash Scripts/pairtopair_slop.sh breakpointworkchco/all_calls.GATK.DUP.sorted.bedpe breakpointworkchco/all_calls.Savvy.DUP.sorted.bedpe  0.50 breakpointworkchco/all_calls.GATK.DUP.x_Savvy.DUP.bedpe.tmp 
            cat breakpointworkchco/all_calls.GATK.DUP.x_Savvy.DUP.bedpe.tmp  | awk '$12 != $25' > breakpointworkchco/all_calls.GATK.DUP.x_Savvy.DUP.bedpe

            bash Scripts/pairtopair_slop.sh breakpointworkchco/all_calls.GATK.DUP.sorted.bedpe breakpointworkchco/all_calls.CNVkit.DUP.sorted.bedpe  0.50 breakpointworkchco/all_calls.GATK.DUP.x_CNVkit.DUP.bedpe.tmp 
            cat reakpointwork/all_calls.GATK.DUP.x_CNVkit.DUP.bedpe.tmp | awk '$12 != $25' > breakpointworkchco/all_calls.GATK.DUP.x_CNVkit.DUP.bedpe

            bash Scripts/pairtopair_slop.sh breakpointworkchco/all_calls.Savvy.DUP.sorted.bedpe breakpointworkchco/all_calls.CNVkit.DUP.sorted.bedpe  0.50 breakpointworkchco/all_calls.Savvy.DUP.x_CNVkit.DUP.bedpe.tmp 
            cat breakpointworkchco/all_calls.Savvy.DUP.x_CNVkit.DUP.bedpe.tmp | awk '$12 != $25' > breakpointworkchco/all_calls.Savvy.DUP.x_CNVkit.DUP.bedpe

            cat breakpointworkchco/all_calls.GATK.DEL.x_Savvy.DEL.bedpe | python Scripts/convert_bedpe_to_bed.py > breakpointworkchco/all_calls.GATK.DEL.x_Savvy.DEL.bed
            cat breakpointworkchco/all_calls.GATK.DEL.x_CNVkit.DEL.bedpe | python Scripts/convert_bedpe_to_bed.py > breakpointworkchco/all_calls.GATK.DEL.x_CNVkit.DEL.bed
            cat breakpointworkchco/all_calls.Savvy.DEL.x_CNVkit.DEL.bedpe | python Scripts/convert_bedpe_to_bed.py > breakpointworkchco/all_calls.Savvy.DEL.x_CNVkit.DEL.bed

            cat breakpointworkchco/all_calls.GATK.DUP.x_Savvy.DUP.bedpe | python Scripts/convert_bedpe_to_bed.py > breakpointworkchco/all_calls.GATK.DUP.x_Savvy.DUP.bed
            cat breakpointworkchco/all_calls.GATK.DUP.x_CNVkit.DUP.bedpe | python Scripts/convert_bedpe_to_bed.py > breakpointworkchco/all_calls.GATK.DUP.x_CNVkit.DUP.bed
            cat breakpointworkchco/all_calls.Savvy.DUP.x_CNVkit.DUP.bedpe | python Scripts/convert_bedpe_to_bed.py > breakpointworkchco/all_calls.Savvy.DUP.x_CNVkit.DUP.bed

            python Scripts/get_triples.py -a breakpointworkchco/all_calls.GATK.DEL.x_Savvy.DEL.bed -b breakpointworkchco/all_calls.GATK.DEL.x_CNVkit.DEL.bed -c breakpointworkchco/all_calls.Savvy.DEL.x_CNVkit.DEL.bed -o breakpointworkchco/triple_overlaps.DEL.bed
            python Scripts/get_triples.py -a breakpointworkchco/all_calls.GATK.DUP.x_Savvy.DUP.bed -b breakpointworkchco/all_calls.GATK.DUP.x_CNVkit.DUP.bed -c breakpointworkchco/all_calls.Savvy.DUP.x_CNVkit.DUP.bed -o breakpointworkchco/triple_overlaps.DUP.bed
        """