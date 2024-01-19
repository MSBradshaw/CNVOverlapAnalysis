cnv_calls = 'Data/all_calls.chco.bed'
SV_calls = 'Data/all_SVs.bed'

rule all:
    input:
        'workchco/triple_overlaps.DEL.bed',
        'workchco/triple_overlaps.DUP.bed',
        'workchco/double_overlaps.bed.txt',
        'workchco/overlap_stats.txt',
        'workchco/upset_plot.DEL.png',
        'workchco/upset_plot.DUP.png',
        'Figures/chco_calls_per_sample_and_sizes.png'
        
rule sort:
    input:
        cnv_calls = cnv_calls
    output:
        cnv_calls_sorted = 'workchco/all_calls.sorted.bed'
    shell:
        """
        mkdir -p workchco
        bedtools sort -i {input.cnv_calls} > {output.cnv_calls_sorted}
        """

rule get_double_overlaps:
    input:
        cnv_calls = 'workchco/all_calls.sorted.bed',
    output:
        double_overlaps = 'workchco/double_overlaps.bed.txt'
    shell:
        """
        # overlap cnvs with cnvs with 60% reciprocal overlap
        bedtools intersect -a {input.cnv_calls} -b {input.cnv_calls} -wao -f .6 -r | cut -f 1-12 > workchco/all_cnv_x_cnv.bed.tmp
        cat workchco/all_cnv_x_cnv.bed.tmp | python Scripts/remove_call_to_call_overlaps.py > workchco/all_cnv_x_cnv.bed
        cat workchco/all_cnv_x_cnv.bed | python Scripts/make_sample_specific_and_caller_disjoint.py > {output.double_overlaps}        
        """

rule split_doubles:
    input:
        double_overlaps = 'workchco/double_overlaps.bed.txt'
    output:
        split_doubles = 'workchco/split_doubles.bed'
    shell:
        """
        cat {input.double_overlaps} | python Scripts/split_doubles.py > {output.split_doubles}
        """

def get_calls_from_bed(filename):
    cnvs = set()
    for line in open(filename):
        cnvs.add(line.strip())
    return cnvs

rule generate_stats_about_single_overlaps:
    input:
        valid_singles = 'workchco/single_overlaps_and_sv.bed',
        all_doubles = 'workchco/split_doubles.bed',
        all_singles = 'workchco/all_calls.sorted.bed'
    output:
        stats = 'workchco/single_stats.txt'
    run:
        with open(output[0],'w') as out:
            out.write('caller\tcnv_type\percent_validated\n')
            for caller in ['GATK','CNVkit','Savvy']:
                for cnv_type in ['DEL','DUP']:
                    count = 0
                    specific_valid_singles = [ x.strip() for x in open(input.valid_singles) if caller in x and cnv_type in x]
                    specific_all_doubles = [x.strip() for x in open(input.all_doubles) if caller in x and cnv_type in x]
                    specific_all_singles = [x.strip() for x in open(input.all_singles) if caller in x and cnv_type in x]
                    # assert that every specific_valid_singles is in specific_all_singles
                    for x in specific_valid_singles:
                        assert x in specific_all_singles, x
                    print(cnv_type,caller)
                    print(len(specific_valid_singles),len(specific_all_singles))
                    
                    specific_valid_singles_only = [x for x in specific_valid_singles if x not in specific_all_doubles]
                    specific_all_singles_only = [x for x in specific_all_singles if x not in specific_all_doubles]

                    if len(specific_valid_singles_only)  == 0 and len(specific_all_singles_only) == 0:
                        print('no singles')
                        continue

                    print(len(specific_valid_singles_only),len(specific_all_singles_only))
                    print(caller,cnv_type,len(specific_valid_singles_only)/len(specific_all_singles_only),sep='\t')
                    out.write('\t'.join([caller,cnv_type,str(len(specific_valid_singles_only)/len(specific_all_singles_only))])+'\n')

rule triple_overlap:
    input:
        cnv_calls_sorted = 'workchco/all_calls.sorted.bed'
    output:
        'workchco/triple_overlaps.DEL.bed',
        'workchco/triple_overlaps.DUP.bed'
    shell:
        """
            grep GATK workchco/all_calls.sorted.bed > workchco/all_calls.GATK.sorted.bed
            grep CNVkit workchco/all_calls.sorted.bed > workchco/all_calls.CNVkit.sorted.bed
            grep Savvy workchco/all_calls.sorted.bed > workchco/all_calls.Savvy.sorted.bed

            grep DEL workchco/all_calls.GATK.sorted.bed > workchco/all_calls.GATK.DEL.sorted.bed
            grep DUP workchco/all_calls.GATK.sorted.bed > workchco/all_calls.GATK.DUP.sorted.bed
            grep DEL workchco/all_calls.CNVkit.sorted.bed > workchco/all_calls.CNVkit.DEL.sorted.bed
            grep DUP workchco/all_calls.CNVkit.sorted.bed > workchco/all_calls.CNVkit.DUP.sorted.bed
            grep DEL workchco/all_calls.Savvy.sorted.bed > workchco/all_calls.Savvy.DEL.sorted.bed
            grep DUP workchco/all_calls.Savvy.sorted.bed > workchco/all_calls.Savvy.DUP.sorted.bed

            bedtools intersect -a workchco/all_calls.GATK.DEL.sorted.bed -b workchco/all_calls.Savvy.DEL.sorted.bed -f .6 -r -wao | cut -f 1-12 | python Scripts/make_sample_specific.py > workchco/all_calls.GATK.DEL.x_Savvy.DEL.bed
            bedtools intersect -a workchco/all_calls.GATK.DEL.sorted.bed -b workchco/all_calls.CNVkit.DEL.sorted.bed -f .6 -r -wao | cut -f 1-12 | python Scripts/make_sample_specific.py > workchco/all_calls.GATK.DEL.x_CNVkit.DEL.bed
            bedtools intersect -a workchco/all_calls.Savvy.DEL.sorted.bed -b workchco/all_calls.CNVkit.DEL.sorted.bed -f .6 -r -wao | cut -f 1-12 | python Scripts/make_sample_specific.py > workchco/all_calls.Savvy.DEL.x_CNVkit.DEL.bed

            bedtools intersect -a workchco/all_calls.GATK.DUP.sorted.bed -b workchco/all_calls.Savvy.DUP.sorted.bed -f .6 -r -wao | cut -f 1-12 | python Scripts/make_sample_specific.py > workchco/all_calls.GATK.DUP.x_Savvy.DUP.bed
            bedtools intersect -a workchco/all_calls.GATK.DUP.sorted.bed -b workchco/all_calls.CNVkit.DUP.sorted.bed -f .6 -r -wao | cut -f 1-12 | python Scripts/make_sample_specific.py > workchco/all_calls.GATK.DUP.x_CNVkit.DUP.bed
            bedtools intersect -a workchco/all_calls.Savvy.DUP.sorted.bed -b workchco/all_calls.CNVkit.DUP.sorted.bed -f .6 -r -wao | cut -f 1-12 | python Scripts/make_sample_specific.py > workchco/all_calls.Savvy.DUP.x_CNVkit.DUP.bed

            python Scripts/get_triples.py -a workchco/all_calls.GATK.DEL.x_Savvy.DEL.bed -b workchco/all_calls.GATK.DEL.x_CNVkit.DEL.bed -c workchco/all_calls.Savvy.DEL.x_CNVkit.DEL.bed -o workchco/triple_overlaps.DEL.bed
            python Scripts/get_triples.py -a workchco/all_calls.GATK.DUP.x_Savvy.DUP.bed -b workchco/all_calls.GATK.DUP.x_CNVkit.DUP.bed -c workchco/all_calls.Savvy.DUP.x_CNVkit.DUP.bed -o workchco/triple_overlaps.DUP.bed
        """

rule generate_stats_about_overlaps:
    input:
        all_overlaps = 'workchco/all_calls.sorted.bed',
        doubles = 'workchco/double_overlaps.bed.txt',
        triples_del = 'workchco/triple_overlaps.DEL.bed',
        triples_dup = 'workchco/triple_overlaps.DUP.bed'
    output:
        stats = 'workchco/overlap_stats.txt'
    shell:
        """
        # remove doubles that overlap a DEL and a DUP
        grep -v 'DUP' {input.doubles} > workchco/double_overlaps.DEL.bed
        grep -v 'DEL' {input.doubles} > workchco/double_overlaps.DUP.bed
        cat workchco/double_overlaps.DEL.bed > workchco/double_overlaps.type_specific.bed
        cat workchco/double_overlaps.DUP.bed >> workchco/double_overlaps.type_specific.bed

        # combine triples
        cat {input.triples_del} > workchco/triple_overlaps.type_specific.bed
        cat {input.triples_dup} >> workchco/triple_overlaps.type_specific.bed

        python Scripts/generate_stats.py --singles {input.all_overlaps} --double workchco/double_overlaps.type_specific.bed --triples workchco/triple_overlaps.type_specific.bed --output {output.stats}
        """

rule upset_plots:
    input:
        stats = 'workchco/overlap_stats.txt'
    output:
        upset_plot_del = 'workchco/upset_plot.DEL.png',
        upset_plot_dup = 'workchco/upset_plot.DUP.png'
    shell:
        """
        # replace CNVkit with CNVKit in stats file
        sed -i 's/CNVkit/CNVKit/g' {input.stats}

        grep DEL {input.stats} | cut -f2 > workchco/DEL.stats.txt
        python Scripts/upset_plot.py -i workchco/DEL.stats.txt -l Deletion -o {output.upset_plot_del}
        grep DUP {input.stats} | cut -f2 > workchco/DUP.stats.txt
        python Scripts/upset_plot.py -i workchco/DUP.stats.txt -l Duplication -o {output.upset_plot_dup}
        """

rule plot_sizes:
    input:
        all_calls=cnv_calls
    output:
        calls='Figures/chco_calls_per_sample_and_sizes.png'
    shell:
        """
        mkdir -p workchco/CallerSpecificCNVTypes
        mkdir -p Figures/
        # cat inputs into single caller specific input tmp files
        grep 'GATK' {input.all_calls} > workchco/CallerSpecificCNVTypes/tmp.GATK.bed
        grep 'CNVkit' {input.all_calls} > workchco/CallerSpecificCNVTypes/tmp.CNVkit.bed
        grep 'Savvy' {input.all_calls} > workchco/CallerSpecificCNVTypes/tmp.Savvy.bed

        # cat all dels into one file
        grep DEL {input.all_calls} > workchco/CallerSpecificCNVTypes/tmp.DEL.bed
        grep DUP {input.all_calls} > workchco/CallerSpecificCNVTypes/tmp.DUP.bed

        # replace GATK with gCNV in workchco/CallerSpecificCNVTypes/tmp.DEL.bed  and workchco/CallerSpecificCNVTypes/tmp.DUP.bed 
        sed -i 's/GATK/gCNV/g' workchco/CallerSpecificCNVTypes/tmp.DEL.bed
        sed -i 's/GATK/gCNV/g' workchco/CallerSpecificCNVTypes/tmp.DUP.bed

        python Scripts/plot_sizes_and_calls_per_sample.py --dels workchco/CallerSpecificCNVTypes/tmp.DEL.bed --dups workchco/CallerSpecificCNVTypes/tmp.DUP.bed --output {output.calls}
        """