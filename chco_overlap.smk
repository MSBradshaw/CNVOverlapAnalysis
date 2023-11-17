cnv_calls = 'Data/all_calls.chco.bed'
SV_calls = 'Data/all_SVs.bed'

rule all:
    input:
        'debugworkchco/triple_overlaps.DEL.bed',
        'debugworkchco/triple_overlaps.DUP.bed',
        'debugworkchco/double_overlaps.bed.txt',
        'debugworkchco/overlap_stats.txt',
        'debugworkchco/upset_plot.DEL.png',
        'debugworkchco/upset_plot.DUP.png',
        'CHCOFigures/chco_size_distribution.png',
        'CHCOFigures/chco_calls_per_sample.png'
        
rule sort:
    input:
        cnv_calls = cnv_calls
    output:
        cnv_calls_sorted = 'debugworkchco/all_calls.sorted.bed'
    shell:
        """
        mkdir -p debugworkchco
        bedtools sort -i {input.cnv_calls} > {output.cnv_calls_sorted}
        """

rule get_double_overlaps:
    input:
        cnv_calls = 'debugworkchco/all_calls.sorted.bed',
    output:
        double_overlaps = 'debugworkchco/double_overlaps.bed.txt'
    shell:
        """
        # overlap cnvs with cnvs with 60% reciprocal overlap
        bedtools intersect -a {input.cnv_calls} -b {input.cnv_calls} -wao -f .6 -r | cut -f 1-12 > debugworkchco/all_cnv_x_cnv.bed.tmp
        cat debugworkchco/all_cnv_x_cnv.bed.tmp | python Scripts/remove_call_to_call_overlaps.py > debugworkchco/all_cnv_x_cnv.bed
        cat debugworkchco/all_cnv_x_cnv.bed | python Scripts/make_sample_specific_and_caller_disjoint.py > {output.double_overlaps}        
        """

rule split_doubles:
    input:
        double_overlaps = 'debugworkchco/double_overlaps.bed.txt'
    output:
        split_doubles = 'debugworkchco/split_doubles.bed'
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
        valid_singles = 'debugworkchco/single_overlaps_and_sv.bed',
        all_doubles = 'debugworkchco/split_doubles.bed',
        all_singles = 'debugworkchco/all_calls.sorted.bed'
    output:
        stats = 'debugworkchco/single_stats.txt'
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
        cnv_calls_sorted = 'debugworkchco/all_calls.sorted.bed'
    output:
        'debugworkchco/triple_overlaps.DEL.bed',
        'debugworkchco/triple_overlaps.DUP.bed'
    shell:
        """
            grep GATK debugworkchco/all_calls.sorted.bed > debugworkchco/all_calls.GATK.sorted.bed
            grep CNVkit debugworkchco/all_calls.sorted.bed > debugworkchco/all_calls.CNVkit.sorted.bed
            grep Savvy debugworkchco/all_calls.sorted.bed > debugworkchco/all_calls.Savvy.sorted.bed

            grep DEL debugworkchco/all_calls.GATK.sorted.bed > debugworkchco/all_calls.GATK.DEL.sorted.bed
            grep DUP debugworkchco/all_calls.GATK.sorted.bed > debugworkchco/all_calls.GATK.DUP.sorted.bed
            grep DEL debugworkchco/all_calls.CNVkit.sorted.bed > debugworkchco/all_calls.CNVkit.DEL.sorted.bed
            grep DUP debugworkchco/all_calls.CNVkit.sorted.bed > debugworkchco/all_calls.CNVkit.DUP.sorted.bed
            grep DEL debugworkchco/all_calls.Savvy.sorted.bed > debugworkchco/all_calls.Savvy.DEL.sorted.bed
            grep DUP debugworkchco/all_calls.Savvy.sorted.bed > debugworkchco/all_calls.Savvy.DUP.sorted.bed

            bedtools intersect -a debugworkchco/all_calls.GATK.DEL.sorted.bed -b debugworkchco/all_calls.Savvy.DEL.sorted.bed -f .6 -r -wao | cut -f 1-12 | python Scripts/make_sample_specific.py > debugworkchco/all_calls.GATK.DEL.x_Savvy.DEL.bed
            bedtools intersect -a debugworkchco/all_calls.GATK.DEL.sorted.bed -b debugworkchco/all_calls.CNVkit.DEL.sorted.bed -f .6 -r -wao | cut -f 1-12 | python Scripts/make_sample_specific.py > debugworkchco/all_calls.GATK.DEL.x_CNVkit.DEL.bed
            bedtools intersect -a debugworkchco/all_calls.Savvy.DEL.sorted.bed -b debugworkchco/all_calls.CNVkit.DEL.sorted.bed -f .6 -r -wao | cut -f 1-12 | python Scripts/make_sample_specific.py > debugworkchco/all_calls.Savvy.DEL.x_CNVkit.DEL.bed

            bedtools intersect -a debugworkchco/all_calls.GATK.DUP.sorted.bed -b debugworkchco/all_calls.Savvy.DUP.sorted.bed -f .6 -r -wao | cut -f 1-12 | python Scripts/make_sample_specific.py > debugworkchco/all_calls.GATK.DUP.x_Savvy.DUP.bed
            bedtools intersect -a debugworkchco/all_calls.GATK.DUP.sorted.bed -b debugworkchco/all_calls.CNVkit.DUP.sorted.bed -f .6 -r -wao | cut -f 1-12 | python Scripts/make_sample_specific.py > debugworkchco/all_calls.GATK.DUP.x_CNVkit.DUP.bed
            bedtools intersect -a debugworkchco/all_calls.Savvy.DUP.sorted.bed -b debugworkchco/all_calls.CNVkit.DUP.sorted.bed -f .6 -r -wao | cut -f 1-12 | python Scripts/make_sample_specific.py > debugworkchco/all_calls.Savvy.DUP.x_CNVkit.DUP.bed

            python Scripts/get_triples.py -a debugworkchco/all_calls.GATK.DEL.x_Savvy.DEL.bed -b debugworkchco/all_calls.GATK.DEL.x_CNVkit.DEL.bed -c debugworkchco/all_calls.Savvy.DEL.x_CNVkit.DEL.bed -o debugworkchco/triple_overlaps.DEL.bed
            python Scripts/get_triples.py -a debugworkchco/all_calls.GATK.DUP.x_Savvy.DUP.bed -b debugworkchco/all_calls.GATK.DUP.x_CNVkit.DUP.bed -c debugworkchco/all_calls.Savvy.DUP.x_CNVkit.DUP.bed -o debugworkchco/triple_overlaps.DUP.bed
        """

rule generate_stats_about_overlaps:
    input:
        all_overlaps = 'debugworkchco/all_calls.sorted.bed',
        doubles = 'debugworkchco/double_overlaps.bed.txt',
        triples_del = 'debugworkchco/triple_overlaps.DEL.bed',
        triples_dup = 'debugworkchco/triple_overlaps.DUP.bed'
    output:
        stats = 'debugworkchco/overlap_stats.txt'
    shell:
        """
        # remove doubles that overlap a DEL and a DUP
        grep -v 'DUP' {input.doubles} > debugworkchco/double_overlaps.DEL.bed
        grep -v 'DEL' {input.doubles} > debugworkchco/double_overlaps.DUP.bed
        cat debugworkchco/double_overlaps.DEL.bed > debugworkchco/double_overlaps.type_specific.bed
        cat debugworkchco/double_overlaps.DUP.bed >> debugworkchco/double_overlaps.type_specific.bed

        # combine triples
        cat {input.triples_del} > debugworkchco/triple_overlaps.type_specific.bed
        cat {input.triples_dup} >> debugworkchco/triple_overlaps.type_specific.bed

        python Scripts/generate_stats.py --singles {input.all_overlaps} --double debugworkchco/double_overlaps.type_specific.bed --triples debugworkchco/triple_overlaps.type_specific.bed --output {output.stats}
        """

rule upset_plots:
    input:
        stats = 'debugworkchco/overlap_stats.txt'
    output:
        upset_plot_del = 'debugworkchco/upset_plot.DEL.png',
        upset_plot_dup = 'debugworkchco/upset_plot.DUP.png'
    shell:
        """
        # replace CNVkit with CNVKit in stats file
        sed -i 's/CNVkit/CNVKit/g' {input.stats}

        grep DEL {input.stats} | cut -f2 > debugworkchco/DEL.stats.txt
        python Scripts/upset_plot.py -i debugworkchco/DEL.stats.txt -l Deletion -o {output.upset_plot_del}
        grep DUP {input.stats} | cut -f2 > debugworkchco/DUP.stats.txt
        python Scripts/upset_plot.py -i debugworkchco/DUP.stats.txt -l Duplication -o {output.upset_plot_dup}
        """

rule plot_sizes:
    input:
        all_calls=cnv_calls
    output:
        size='CHCOFigures/chco_size_distribution.png',
        calls='CHCOFigures/chco_calls_per_sample.png'
    shell:
        """
        mkdir -p workchco/CallerSpecificCNVTypes
        mkdir -p Figures/
        # cat inputs into single caller specific input tmp files
        grep 'GATK' {input.all_calls} > workchco/CallerSpecificCNVTypes/tmp.GATK.bed
        grep 'CNVkit' {input.all_calls} > workchco/CallerSpecificCNVTypes/tmp.CNVkit.bed
        grep 'Savvy' {input.all_calls} > workchco/CallerSpecificCNVTypes/tmp.Savvy.bed

        python Scripts/plot_sizes.py -i workchco/CallerSpecificCNVTypes/tmp.GATK.bed \
            -i workchco/CallerSpecificCNVTypes/tmp.Savvy.bed \
            -i workchco/CallerSpecificCNVTypes/tmp.CNVkit.bed \
            -l gCNV \
            -l Savvy \
            -l CNVkit \
            -o {output.size}


        # cat all dels into one file
        grep DEL {input.all_calls} > workchco/CallerSpecificCNVTypes/tmp.DEL.bed
        grep DUP {input.all_calls} > workchco/CallerSpecificCNVTypes/tmp.DUP.bed

        # replace GATK with gCNV in workchco/CallerSpecificCNVTypes/tmp.DEL.bed  and workchco/CallerSpecificCNVTypes/tmp.DUP.bed 
        sed -i 's/GATK/gCNV/g' workchco/CallerSpecificCNVTypes/tmp.DEL.bed
        sed -i 's/GATK/gCNV/g' workchco/CallerSpecificCNVTypes/tmp.DUP.bed

        python Scripts/plot_calls_per_sample.py --dels workchco/CallerSpecificCNVTypes/tmp.DEL.bed --dups workchco/CallerSpecificCNVTypes/tmp.DUP.bed --output {output.calls}
        """