import numpy as np
import pandas as pd

CNV_TYPES = ['DEL', 'DUP']
ENDS=['savvy','cnvkit','gatk','triple','savvy_gatk','savvy_cnvkit','cnvkit_gatk']
rule all:
    input:
        expand('work/ValidatedCallsByCategory/{cnv_type}-{end}.bed.txt',cnv_type=CNV_TYPES,end=ENDS),
        expand('work/ValidatedCallsByCategoryBedPE/{cnv_type}-{end}.bed.txt',cnv_type=CNV_TYPES,end=ENDS),
        'work/ValidatedCallsByCategoryStats/agg_stats.txt',
        expand('work/ValidatedCallsByCategoryStats/for_upset_{ol_type}.{cnv_type}.txt',cnv_type=['DEL','DUP'],ol_type=['bed','bedpe']),
        expand('Figures/upsetplot_sv_validation_{ol_type}.{cnv_type}.png',cnv_type=['DEL','DUP'],ol_type=['bed','bedpe']),
        expand('work/ValidatedCallsByCategoryStats/tabulated_validation_percents_{ol_type}.{cnv_type}.tsv',cnv_type=['DEL','DUP'],ol_type=['bed','bedpe']),
        expand('Figures/venn_diagram_{ol_type}.{cnv_type}.png',cnv_type=['DEL','DUP'],ol_type=['bed','bedpe']),
        expand('work/MislabeledCalls/{cnv_type}-{end}.mismatched_calls.bed',cnv_type=CNV_TYPES,end=['savvy','cnvkit','gatk']),
        expand('work/MislabeledCalls/mislabeled_triples.{cnv_type}.bed',cnv_type=CNV_TYPES),
        'work/SV_Validation/overall_validation_stats.txt',
        'work/SV_Validation/overlap_validation_stats.txt'


rule get_percent_of_recip_calls_validated_by_sv:
    input:
        bed='work/RecipCategoriesBeds/{cnv_type}-{end}.bed',
        bedpe='work/RecipCategoriesBedPEs/{cnv_type}-{end}.bed',
        real_bed='work/SV_Validation/real_calls.{cnv_type}.bed'
    output:
        bed='work/ValidatedCallsByCategory/{cnv_type}-{end}.bed.txt',
        bedpe='work/ValidatedCallsByCategoryBedPE/{cnv_type}-{end}.bed.txt',
        stats='work/ValidatedCallsByCategoryStats/{cnv_type}-{end}.stats'
    shell:
        """
        mkdir -p work/ValidatedCallsByCategory
        mkdir -p work/ValidatedCallsByCategoryBedPE
        mkdir -p work/ValidatedCallsByCategoryStats

        # make a list of sample specific calls
        bedtools intersect -a {input.real_bed} -b {input.bed} -f .9 -r -wb | python Scripts/make_sample_specific.py | sort | uniq > {output.bed}
        # unconvert from bedpe to bed
        cat {input.bedpe} | awk '{{print $1"\t"$2+1"\t"$5+1"\t"$11"\t"$12"\t"$13}}' > {input.bedpe}.tmp
        bedtools intersect -a {input.real_bed} -b {input.bedpe}.tmp -f .9 -r -wb  | python Scripts/make_sample_specific.py | sort | uniq  > {output.bedpe}
        rm {input.bedpe}.tmp
        # to a file write the number of lines in {input.bed}, {output.bed} the differences
        python Scripts/calc_bed_difference_stats.py -a {input.bed} -b {output.bed} -l "bed.{wildcards.cnv_type}-{wildcards.end}"  > {output.stats}
        python Scripts/calc_bed_difference_stats.py -a {input.bedpe} -b {output.bedpe} -l "bedpe.{wildcards.cnv_type}-{wildcards.end}"  >> {output.stats}
        """

rule agg_stats:
    input:
        expand('work/ValidatedCallsByCategoryStats/{cnv_type}-{end}.stats',cnv_type=CNV_TYPES,end=ENDS)
    output:
        'work/ValidatedCallsByCategoryStats/agg_stats.txt'
    shell:
        """
        cat {input} > {output}
        """

"""
Savvy: 18423
CNVKit: 3286
GATK: 15196
Savvy-GATK: 757
Savvy-CNVKit: 133
CNVKit-GATK: 22
Triple: 18
"""

rule agg_stats_for_upset_plot:
    input:
        'work/ValidatedCallsByCategoryStats/agg_stats.txt'
    output:
        'work/ValidatedCallsByCategoryStats/for_upset_{ol_type}.{cnv_type}.txt'
    params:
        ol_type = lambda wildcards, output: output[0].split('.')[0].split('_')[-1],
        call_type=lambda wildcards, output: output[0].split('.')[1].replace('.txt','')
    run:
        with open(output[0],'w') as outfile:
            for line in open(input[0]):
                if params.ol_type + '.' in line and params.call_type in line:
                    row = line.strip().split('\t')
                    group = row[0].split('-')[-1].replace('savvy','Savvy').replace('gatk','GATK').replace('cnvkit','CNVKit').replace('triple','Triple').replace('_','-') + ': '
                    percent = row[-1]
                    diff = str(int(row[1]) - int(row[2]))
                    outfile.write(group + percent + '\n')

rule tabulate_agg_stats:
    input:
        'work/ValidatedCallsByCategoryStats/agg_stats.txt'
    output:
        'work/ValidatedCallsByCategoryStats/tabulated_validation_percents_{ol_type}.{cnv_type}.tsv'
    params:
        ol_type = lambda wildcards, output: output[0].split('.')[0].split('_')[-1],
        call_type=lambda wildcards, output: output[0].split('.')[1].replace('.txt','')
    run:
        with open(output[0],'w') as outfile:
            for line in open(input[0]):
                if params.ol_type + '.' in line and params.call_type in line:
                    row = line.strip().split('\t')
                    group = row[0].split('-')[-1].replace('savvy','Savvy').replace('gatk','GATK').replace('cnvkit','CNVKit').replace('triple','Triple').replace('_','-') + '\t'
                    percent = row[-1]
                    diff = str(int(row[1]) - int(row[2]))
                    outfile.write(group + row[1] + '\t' + percent + '\n')

rule upset_plot_percent_of_calls_validated_by_sv:
    input:
        'work/ValidatedCallsByCategoryStats/for_upset_{ol_type}.{cnv_type}.txt' 
    output:
        'Figures/upsetplot_sv_validation_{ol_type}.{cnv_type}.png'
    params: 
        label = lambda wildcards, output: wildcards.ol_type.replace('bedpe','Breakpoint').replace('bed','Reciprocal') + ' ' + wildcards.cnv_type.replace('DEL',' Deletions').replace('DUP',' Duplications')
    shell:
        """
        python Scripts/upset_plot.py -i {input} -l "{params.label}" -o {output} --logscale 0
        """

rule plot_venn_diagrams:
    input:
        'work/ValidatedCallsByCategoryStats/tabulated_validation_percents_{ol_type}.{cnv_type}.tsv'
    output:
        'Figures/venn_diagram_{ol_type}.{cnv_type}.png'
    shell:
        """
        python Scripts/plot_venn_diagram.py -i {input} -o {output}
        """

# how many calls have the wrong label? Were marked dups by SVs but dels by cnvs?
rule find_mismatched_call_types:
    input:
        bed_calls='work/RecipCategoriesBeds/{cnv_type}-{end}.bed',
        validated_beds='work/ValidatedCallsByCategory/{cnv_type}-{end}.bed.txt',
        real_bed='work/SV_Validation/real_calls.{cnv_type}.bed'
    output:
        'work/MislabeledCalls/{cnv_type}-{end}.mismatched_calls.bed'
    params:
        real_bed_of_oppisite_type = lambda wildcards: 'work/SV_Validation/real_calls.' + ('DEL' if wildcards.cnv_type == 'DUP' else 'DUP') + '.bed'
    shell:
        """
        mkdir -p work/MislabeledCalls
        # find the calls not validated by SVs
        python Scripts/get_non_validated_calls.py -a {input.bed_calls} -b {input.validated_beds} -o {output}.tmp.bed

        bedtools intersect -a {params.real_bed_of_oppisite_type} -b {output}.tmp.bed -f .9 -wb | python Scripts/make_sample_specific.py | sort | uniq > {output}
        """

# is there any over lap between the mislabeled calls and the triple overlaps?
rule overlap_mismatched_and_triples:
    input:
        mismatched=expand('work/MislabeledCalls/{cnv_type}-{end}.mismatched_calls.bed',cnv_type='{cnv_type}',end=['savvy','cnvkit','gatk']),
        triples='work/VenDiagramResults/overlap_numbers.{cnv_type}.sample_specific.triples.bed.txt'
    output:
        'work/MislabeledCalls/mislabeled_triples.{cnv_type}.bed'
    shell:
        """
        cat {input.mismatched} > work/MislabeledCalls/all.{wildcards.cnv_type}.mislabled.bed
        bedtools sort -i work/MislabeledCalls/all.{wildcards.cnv_type}.mislabled.bed > work/MislabeledCalls/all.{wildcards.cnv_type}.mislabled.sorted.bed
        bedtools intersect -a work/MislabeledCalls/all.{wildcards.cnv_type}.mislabled.sorted.bed -b {input.triples} -f .99 -wb | python Scripts/make_sample_specific.py | sort | uniq > {output}
        """

def count_lines_in_file(file,caller=None):
    count = 0
    for line in open(file,'r'):
        if caller is None:
            count += 1
        elif caller in line:
            count += 1
    return count

rule overall_validation_stats:
    input:
        validated_dels='work/SV_Validation/real_calls.DEL.bed',
        validated_dups='work/SV_Validation/real_calls.DUP.bed',
        og_dups='work/CNVTypes/cnv_calls.DEL.bed',
        og_dels='work/CNVTypes/cnv_calls.DEL.bed'
    output:
        'work/SV_Validation/overall_validation_stats.txt'
    run:
        with open(output[0],'w') as outfile:
            outfile.write('\t'.join(['','DEL','DUP']) + '\n')
            for caller in ['GATK','Savvy','CNVkit']:
                outfile.write(caller + '\t')
                for cnv_type in ['DEL','DUP']:
                    validated = count_lines_in_file('work/SV_Validation/real_calls.{}.bed'.format(cnv_type),caller=caller)
                    og = count_lines_in_file('work/CNVTypes/cnv_calls.{}.bed'.format(cnv_type),caller=caller)
                    print(caller, cnv_type)
                    print('work/SV_Validation/real_calls.{}.bed'.format(cnv_type))
                    print(validated)
                    print('work/CNVTypes/cnv_calls.{}.bed'.format(cnv_type))
                    print(og)
                    outfile.write('\t' + str((og - validated)/og))
                outfile.write('\n')

def load_double_overlap_file_as_set(file):
    double_overlaps = set()
    for line in open(file,'r'):
        row = line.strip().split('\t')
        call1 = '\t'.join(row[:6])
        call2 = '\t'.join(row[6:-1])
        if row[4] != row[10]:
            continue
        double_overlaps.add(call1)
        double_overlaps.add(call2)
    return double_overlaps

def load_bed_file_as_set(file):
    bed = set()
    for line in open(file,'r'):
        bed.add(line.strip())
    return bed

# what is the overall percentage of double overlap calls that are validated by SVs?
rule validation_stats_for_overlap:
    input:
        del_overlaps=expand('work/CallerSpecificOverlapsBed/venn.DEL.{caller}.bed.txt',caller=['Savvy','CNVkit','GATK']),
        dup_overlaps=expand('work/CallerSpecificOverlapsBed/venn.DUP.{caller}.bed.txt',caller=['Savvy','CNVkit','GATK']),
        validated_dels='work/SV_Validation/real_calls.DEL.bed',
        validated_dups='work/SV_Validation/real_calls.DUP.bed',
        del_triples='work/VenDiagramResults/overlap_numbers.DEL.sample_specific.triples.bed.txt',
        dup_triples='work/VenDiagramResults/overlap_numbers.DUP.sample_specific.triples.bed.txt'
    output:
        'work/SV_Validation/overlap_validation_stats.txt'
    run:
        # get a set of call doubles
        del_double_overlaps = set()
        for f in input.del_overlaps:
            tmp = load_double_overlap_file_as_set(f)
            del_double_overlaps = del_double_overlaps.union(tmp)
        print(list(del_double_overlaps)[:10])
        # get a set of all duplication doubles
        dup_double_overlaps = set()
        for f in input.dup_overlaps:
            tmp = load_double_overlap_file_as_set(f)
            dup_double_overlaps = dup_double_overlaps.union(tmp)
        # get the set of all triple dels
        del_triples = load_bed_file_as_set(input.del_triples)
        # get the set of all triple dups
        dup_triples = load_bed_file_as_set(input.dup_triples)
        # get the set of all validated dels
        validated_dels = load_bed_file_as_set(input.validated_dels)
        # get the set of all validated dups
        validated_dups = load_bed_file_as_set(input.validated_dups)
        # print the percent of double overlaps that are validated
        percent_double_dels = len(del_double_overlaps.intersection(validated_dels))/len(del_double_overlaps)
        percent_double_dups = len(dup_double_overlaps.intersection(validated_dups))/len(dup_double_overlaps)
        # print the percent of triple overlaps that are validated
        percent_triple_dels = len(del_triples.intersection(validated_dels))/len(del_triples)
        percent_triple_dups = len(dup_triples.intersection(validated_dups))/len(dup_triples)
        with open(output[0],'w') as outfile:
            outfile.write('\t'.join(['','DEL','DUP']) + '\n')
            outfile.write('\t'.join(['Double Overlaps',str(percent_double_dels),str(percent_double_dups)]) + '\n')
            outfile.write('\t'.join(['Triple Overlaps',str(percent_triple_dels),str(percent_triple_dups)]) + '\n')
        with open('double.del.bed','w') as outfile:
            for call in del_double_overlaps:
                outfile.write(call + '\n')
        with open('double.dup.bed','w') as outfile:
            for call in dup_double_overlaps:
                outfile.write(call + '\n')


        
