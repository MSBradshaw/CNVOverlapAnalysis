cnv_calls = 'Data/all_calls.bed'
SV_calls = 'Data/all_SVs.bed'
percentages = ['0.9','0.6','0.5','0.4','0.3','0.2','0.1']

rule all:
    input:
        'debugwork/upset_plot.DEL.png',
        'debugwork/upset_plot.DUP.png',
        expand('debugwork/validated_upset_plot.DEL.{percent}.png',percent=percentages),
        expand('debugwork/validated_upset_plot.DUP.{percent}.png',percent=percentages),
        expand('debugwork/validated_percent_upset_plot.DEL.{percent}.png',percent=percentages),
        expand('debugwork/validated_percent_upset_plot.DUP.{percent}.png',percent=percentages),

rule sort:
    input:
        cnv_calls = cnv_calls,
        SV_calls = SV_calls
    output:
        cnv_calls_sorted = 'debugwork/all_calls.sorted.bed',
        SV_calls_sorted = 'debugwork/all_SVs.sorted.bed'
    shell:
        """
        mkdir -p debugwork
        bedtools sort -i {input.cnv_calls} > {output.cnv_calls_sorted}
        bedtools sort -i {input.SV_calls} > {output.SV_calls_sorted}.tmp
        # strip chr from SVs
        cat {output.SV_calls_sorted}.tmp | sed 's/^chr//' > {output.SV_calls_sorted}
        """

rule get_double_overlaps:
    input:
        cnv_calls = 'debugwork/all_calls.sorted.bed',
    output:
        double_overlaps = 'debugwork/double_overlaps.bed.txt'
    shell:
        """
        # overlap cnvs with cnvs with 60% reciprocal overlap
        bedtools intersect -a {input.cnv_calls} -b {input.cnv_calls} -wao -f .6 -r | cut -f 1-12 > debugwork/all_cnv_x_cnv.bed.tmp
        cat debugwork/all_cnv_x_cnv.bed.tmp | python Scripts/remove_call_to_call_overlaps.py > debugwork/all_cnv_x_cnv.bed
        cat debugwork/all_cnv_x_cnv.bed | python Scripts/make_sample_specific_and_caller_disjoint.py > {output.double_overlaps}        
        """

rule split_doubles:
    input:
        double_overlaps = 'debugwork/double_overlaps.bed.txt'
    output:
        split_doubles = 'debugwork/split_doubles.bed'
    shell:
        """
        cat {input.double_overlaps} | python Scripts/split_doubles.py > {output.split_doubles}
        """

def get_calls_from_bed(filename):
    cnvs = set()
    for line in open(filename):
        cnvs.add(line.strip())
    return cnvs

rule overlap_singles_and_sv:
    input:
        cnv_calls_sorted = 'debugwork/all_calls.sorted.bed',
        SV_calls_sorted = 'debugwork/all_SVs.sorted.bed',
    output:
        'debugwork/single_overlaps_and_sv.{percent}.bed'
    shell:
        """
        bedtools intersect -a {input.cnv_calls_sorted} -b {input.SV_calls_sorted} -wao -f {wildcards.percent} | cut -f 1-12 | python Scripts/make_sample_specific.py | python Scripts/make_type_specific.py | cut -f1,2,3,4,5,6 | sort | uniq > {output}.tmp.bed
        bedtools sort -i  {output}.tmp.bed > {output}
        """

rule overlap_singles_and_sv_all:
    input:
        expand('debugwork/single_overlaps_and_sv.{percent}.bed',percent=percentages)

rule generate_stats_about_single_overlaps:
    input:
        valid_singles = 'debugwork/single_overlaps_and_sv.{percent}.bed',
        all_doubles = 'debugwork/split_doubles.bed',
        all_singles = 'debugwork/all_calls.sorted.bed'
    output:
        stats = 'debugwork/single_stats.{percent}.txt'
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
        cnv_calls_sorted = 'debugwork/all_calls.sorted.bed',
        validated_cnvs = 'debugwork/single_overlaps_and_sv.{percent}.bed'

    output:
        'debugwork/triple_overlaps.DEL.validated.{percent}.bed',
        'debugwork/triple_overlaps.DUP.validated.{percent}.bed',
        'debugwork/triple_overlaps.DUP.{percent}.bed',
        'debugwork/triple_overlaps.DEL.{percent}.bed',
    threads: 64
    shell:
        """
            grep GATK {input.cnv_calls_sorted} > debugwork/all_calls.GATK.sorted.bed
            grep CNVkit {input.cnv_calls_sorted}  > debugwork/all_calls.CNVkit.sorted.bed
            grep Savvy {input.cnv_calls_sorted}  > debugwork/all_calls.Savvy.sorted.bed

            grep DEL debugwork/all_calls.GATK.sorted.bed > debugwork/all_calls.GATK.DEL.sorted.bed
            grep DUP debugwork/all_calls.GATK.sorted.bed > debugwork/all_calls.GATK.DUP.sorted.bed
            grep DEL debugwork/all_calls.CNVkit.sorted.bed > debugwork/all_calls.CNVkit.DEL.sorted.bed
            grep DUP debugwork/all_calls.CNVkit.sorted.bed > debugwork/all_calls.CNVkit.DUP.sorted.bed
            grep DEL debugwork/all_calls.Savvy.sorted.bed > debugwork/all_calls.Savvy.DEL.sorted.bed
            grep DUP debugwork/all_calls.Savvy.sorted.bed > debugwork/all_calls.Savvy.DUP.sorted.bed

            bedtools intersect -a debugwork/all_calls.GATK.DEL.sorted.bed -b debugwork/all_calls.Savvy.DEL.sorted.bed -f .6 -r -wao | cut -f 1-12 | python Scripts/make_sample_specific.py > debugwork/all_calls.GATK.DEL.x_Savvy.DEL.bed
            bedtools intersect -a debugwork/all_calls.GATK.DEL.sorted.bed -b debugwork/all_calls.CNVkit.DEL.sorted.bed -f .6 -r -wao | cut -f 1-12 | python Scripts/make_sample_specific.py > debugwork/all_calls.GATK.DEL.x_CNVkit.DEL.bed
            bedtools intersect -a debugwork/all_calls.Savvy.DEL.sorted.bed -b debugwork/all_calls.CNVkit.DEL.sorted.bed -f .6 -r -wao | cut -f 1-12 | python Scripts/make_sample_specific.py > debugwork/all_calls.Savvy.DEL.x_CNVkit.DEL.bed

            bedtools intersect -a debugwork/all_calls.GATK.DUP.sorted.bed -b debugwork/all_calls.Savvy.DUP.sorted.bed -f .6 -r -wao | cut -f 1-12 | python Scripts/make_sample_specific.py > debugwork/all_calls.GATK.DUP.x_Savvy.DUP.bed
            bedtools intersect -a debugwork/all_calls.GATK.DUP.sorted.bed -b debugwork/all_calls.CNVkit.DUP.sorted.bed -f .6 -r -wao | cut -f 1-12 | python Scripts/make_sample_specific.py > debugwork/all_calls.GATK.DUP.x_CNVkit.DUP.bed
            bedtools intersect -a debugwork/all_calls.Savvy.DUP.sorted.bed -b debugwork/all_calls.CNVkit.DUP.sorted.bed -f .6 -r -wao | cut -f 1-12 | python Scripts/make_sample_specific.py > debugwork/all_calls.Savvy.DUP.x_CNVkit.DUP.bed

            python Scripts/get_triples.py -a debugwork/all_calls.GATK.DEL.x_Savvy.DEL.bed -b debugwork/all_calls.GATK.DEL.x_CNVkit.DEL.bed -c debugwork/all_calls.Savvy.DEL.x_CNVkit.DEL.bed -o debugwork/triple_overlaps.DEL.{wildcards.percent}.bed
            python Scripts/get_triples.py -a debugwork/all_calls.GATK.DUP.x_Savvy.DUP.bed -b debugwork/all_calls.GATK.DUP.x_CNVkit.DUP.bed -c debugwork/all_calls.Savvy.DUP.x_CNVkit.DUP.bed -o debugwork/triple_overlaps.DUP.{wildcards.percent}.bed

            # intersect triples with real calls
            bedtools intersect -a debugwork/triple_overlaps.DEL.{wildcards.percent}.bed -b {input.validated_cnvs} -wao -f .99 -r | cut -f 1-12 | python Scripts/make_sample_specific.py | python Scripts/make_caller_specific.py | cut -f 1-6 | sort | uniq > debugwork/triple_overlaps.DEL.validated.{wildcards.percent}.bed
            bedtools intersect -a debugwork/triple_overlaps.DUP.{wildcards.percent}.bed -b {input.validated_cnvs} -wao -f .99 | cut -f 1-12 | python Scripts/make_sample_specific.py | python Scripts/make_caller_specific.py | cut -f 1-6 | sort | uniq > debugwork/triple_overlaps.DUP.validated.{wildcards.percent}.bed
        """


rule generate_stats_about_overlaps:
    input:
        all_overlaps = 'debugwork/all_calls.sorted.bed',
        doubles = 'debugwork/double_overlaps.bed.txt',
    output:
        stats = 'debugwork/overlap_stats.txt'
    shell:
        """
        # make the triple overlaps
        grep GATK {input.all_overlaps} > debugwork/all_calls.GATK.sorted.bed
        grep CNVkit {input.all_overlaps}  > debugwork/all_calls.CNVkit.sorted.bed
        grep Savvy {input.all_overlaps}  > debugwork/all_calls.Savvy.sorted.bed

        grep DEL debugwork/all_calls.GATK.sorted.bed > debugwork/all_calls.GATK.DEL.sorted.bed
        grep DUP debugwork/all_calls.GATK.sorted.bed > debugwork/all_calls.GATK.DUP.sorted.bed
        grep DEL debugwork/all_calls.CNVkit.sorted.bed > debugwork/all_calls.CNVkit.DEL.sorted.bed
        grep DUP debugwork/all_calls.CNVkit.sorted.bed > debugwork/all_calls.CNVkit.DUP.sorted.bed
        grep DEL debugwork/all_calls.Savvy.sorted.bed > debugwork/all_calls.Savvy.DEL.sorted.bed
        grep DUP debugwork/all_calls.Savvy.sorted.bed > debugwork/all_calls.Savvy.DUP.sorted.bed

        bedtools intersect -a debugwork/all_calls.GATK.DEL.sorted.bed -b debugwork/all_calls.Savvy.DEL.sorted.bed -f .6 -r -wao | cut -f 1-12 | python Scripts/make_sample_specific.py > debugwork/all_calls.GATK.DEL.x_Savvy.DEL.bed
        bedtools intersect -a debugwork/all_calls.GATK.DEL.sorted.bed -b debugwork/all_calls.CNVkit.DEL.sorted.bed -f .6 -r -wao | cut -f 1-12 | python Scripts/make_sample_specific.py > debugwork/all_calls.GATK.DEL.x_CNVkit.DEL.bed
        bedtools intersect -a debugwork/all_calls.Savvy.DEL.sorted.bed -b debugwork/all_calls.CNVkit.DEL.sorted.bed -f .6 -r -wao | cut -f 1-12 | python Scripts/make_sample_specific.py > debugwork/all_calls.Savvy.DEL.x_CNVkit.DEL.bed

        bedtools intersect -a debugwork/all_calls.GATK.DUP.sorted.bed -b debugwork/all_calls.Savvy.DUP.sorted.bed -f .6 -r -wao | cut -f 1-12 | python Scripts/make_sample_specific.py > debugwork/all_calls.GATK.DUP.x_Savvy.DUP.bed
        bedtools intersect -a debugwork/all_calls.GATK.DUP.sorted.bed -b debugwork/all_calls.CNVkit.DUP.sorted.bed -f .6 -r -wao | cut -f 1-12 | python Scripts/make_sample_specific.py > debugwork/all_calls.GATK.DUP.x_CNVkit.DUP.bed
        bedtools intersect -a debugwork/all_calls.Savvy.DUP.sorted.bed -b debugwork/all_calls.CNVkit.DUP.sorted.bed -f .6 -r -wao | cut -f 1-12 | python Scripts/make_sample_specific.py > debugwork/all_calls.Savvy.DUP.x_CNVkit.DUP.bed

        python Scripts/get_triples.py -a debugwork/all_calls.GATK.DEL.x_Savvy.DEL.bed -b debugwork/all_calls.GATK.DEL.x_CNVkit.DEL.bed -c debugwork/all_calls.Savvy.DEL.x_CNVkit.DEL.bed -o debugwork/triple_overlaps.DEL.bed
        python Scripts/get_triples.py -a debugwork/all_calls.GATK.DUP.x_Savvy.DUP.bed -b debugwork/all_calls.GATK.DUP.x_CNVkit.DUP.bed -c debugwork/all_calls.Savvy.DUP.x_CNVkit.DUP.bed -o debugwork/triple_overlaps.DUP.bed


        # remove doubles that overlap a DEL and a DUP
        grep -v 'DUP' {input.doubles} > debugwork/double_overlaps.DEL.bed
        grep -v 'DEL' {input.doubles} > debugwork/double_overlaps.DUP.bed
        cat debugwork/double_overlaps.DEL.bed > debugwork/double_overlaps.type_specific.bed
        cat debugwork/double_overlaps.DUP.bed >> debugwork/double_overlaps.type_specific.bed

        # combine triples
        cat debugwork/triple_overlaps.DEL.bed > debugwork/triple_overlaps.type_specific.bed
        cat debugwork/triple_overlaps.DUP.bed >> debugwork/triple_overlaps.type_specific.bed

        python Scripts/generate_stats.py --singles {input.all_overlaps} --double debugwork/double_overlaps.type_specific.bed --triples debugwork/triple_overlaps.type_specific.bed --output {output.stats}
        """

rule upset_plots:
    input:
        stats = 'debugwork/overlap_stats.txt'
    output:
        upset_plot_del = 'debugwork/upset_plot.DEL.png',
        upset_plot_dup = 'debugwork/upset_plot.DUP.png'
    shell:
        """
        # replace CNVkit with CNVKit in stats file
        sed -i 's/CNVkit/CNVKit/g' {input.stats}

        grep DEL {input.stats} | cut -f2 > debugwork/DEL.stats.txt
        python Scripts/upset_plot.py -i debugwork/DEL.stats.txt -l Deletion -o {output.upset_plot_del}
        grep DUP {input.stats} | cut -f2 > debugwork/DUP.stats.txt
        python Scripts/upset_plot.py -i debugwork/DUP.stats.txt -l Duplication -o {output.upset_plot_dup}
        """

rule generate_stats_about_validated_overlaps:
    input:
        real_calls= 'debugwork/single_overlaps_and_sv.{percent}.bed',

        triples_del = 'debugwork/triple_overlaps.DEL.validated.{percent}.bed',
        triples_dup = 'debugwork/triple_overlaps.DUP.validated.{percent}.bed'
    output:
        stats = 'debugwork/validated_overlap_stats.{percent}.txt'
    shell:
        """
        # get the double overlap real calls
        bedtools intersect -a {input.real_calls} -b {input.real_calls} -wao -f .6 -r | cut -f 1-12 | python Scripts/make_sample_specific.py | python Scripts/make_type_specific.py > debugwork/validated_double_overlaps.{wildcards.percent}.bed.txt

        # remove doubles that overlap a DEL and a DUP
        grep -v 'DUP' debugwork/validated_double_overlaps.{wildcards.percent}.bed.txt > debugwork/validated_double_overlaps.{wildcards.percent}.DEL.bed
        grep -v 'DEL' debugwork/validated_double_overlaps.{wildcards.percent}.bed.txt > debugwork/validated_double_overlaps.{wildcards.percent}.DUP.bed
        cat debugwork/validated_double_overlaps.{wildcards.percent}.DEL.bed > debugwork/validated_double_overlaps.type_specific.{wildcards.percent}.bed
        cat debugwork/validated_double_overlaps.{wildcards.percent}.DUP.bed >> debugwork/validated_double_overlaps.type_specific.{wildcards.percent}.bed

        # combine triples
        cat {input.triples_del} > debugwork/validated_triple_overlaps.type_specific.{wildcards.percent}.bed
        cat {input.triples_dup} >> debugwork/validated_triple_overlaps.type_specific.{wildcards.percent}.bed

        python Scripts/generate_stats.py --singles {input.real_calls} --double debugwork/validated_double_overlaps.type_specific.{wildcards.percent}.bed --triples debugwork/validated_triple_overlaps.type_specific.{wildcards.percent}.bed --output {output.stats}
        """

rule calc_stats_as_portions:
    input:
        stats = 'debugwork/overlap_stats.txt',
        valid_stats = 'debugwork/validated_overlap_stats.{percent}.txt'
    output:
        stats = 'debugwork/validated_percent_overlap_stats.{percent}.txt'
    run:
        # load stats as a dict, split on ': '
        stats = {}
        for line in open(input.stats):
            line = line.strip().split(': ')
            stats[line[0]] = float(line[1])
        # load valid stats in the same way
        valid_stats = {}
        for line in open(input.valid_stats):
            line = line.strip().split(': ')
            valid_stats[line[0]] = float(line[1])
        # divide each valid stat by the stat
        with open(output[0],'w') as out:
            for key in stats:
                out.write(key+': '+str(valid_stats[key]/stats[key])+'\n')

rule upset_plots_valdiated:
    input:
        stats = 'debugwork/validated_overlap_stats.{percent}.txt'
    output:
        upset_plot_del = 'debugwork/validated_upset_plot.DEL.{percent}.png',
        upset_plot_dup = 'debugwork/validated_upset_plot.DUP.{percent}.png'
    shell:
        """
        # replace CNVkit with CNVKit in stats file
        sed -i 's/CNVkit/CNVKit/g' {input.stats}

        grep DEL {input.stats} | cut -f2 > debugwork/DEL.{wildcards.percent}.stats.txt
        python Scripts/upset_plot.py -i debugwork/DEL.{wildcards.percent}.stats.txt -l 'Deletion {wildcards.percent}' -o {output.upset_plot_del}
        grep DUP {input.stats} | cut -f2 > debugwork/DUP.{wildcards.percent}.stats.txt
        python Scripts/upset_plot.py -i debugwork/DUP.{wildcards.percent}.stats.txt -l 'Duplication {wildcards.percent}' -o {output.upset_plot_dup}
        """

rule upset_plots_valdiated_percents:
    input:
        stats = 'debugwork/validated_percent_overlap_stats.{percent}.txt'
    output:
        upset_plot_del = 'debugwork/validated_percent_upset_plot.DEL.{percent}.png',
        upset_plot_dup = 'debugwork/validated_percent_upset_plot.DUP.{percent}.png'
    shell:
        """
        # replace CNVkit with CNVKit in stats file
        sed -i 's/CNVkit/CNVKit/g' {input.stats}

        grep DEL {input.stats} | cut -f2 > debugwork/DEL.stats.percent.{wildcards.percent}.txt
        python Scripts/upset_plot.py -i debugwork/DEL.stats.percent.{wildcards.percent}.txt -l 'Deletion {wildcards.percent}' -o {output.upset_plot_del} --logscale 0 --percent
        grep DUP {input.stats} | cut -f2 > debugwork/DUP.stats.percent.{wildcards.percent}.txt
        python Scripts/upset_plot.py -i debugwork/DUP.stats.percent.{wildcards.percent}.txt -l 'Duplication {wildcards.percent}' -o {output.upset_plot_dup} --logscale 0 --percent
        """
