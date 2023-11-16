cnv_calls = 'Data/all_calls.bed'
SV_calls = 'Data/all_SVs.bed'

# this is the parameter for the percentage of overlap the CNVs and SVs to say a CNV is real. 0.3 is the default and optimal according to my experimentation
if 'percentage_param' in config:
    percentage_param = config['percentage_param']
else:
    percentage_param = '0.3'

outputdir = f'1kg_sv_validation{percentage_param}'

rule all:
    input:
        f'{outputdir}/double_overlaps.bed.txt',
        f'{outputdir}/split_doubles.bed',
        f'{outputdir}/double_overlaps_and_sv.bed',
        f'{outputdir}/double_stats.txt',
        f'{outputdir}/single_stats.txt',
        f'{outputdir}/triple_overlaps.DEL.validated.bed',
        f'{outputdir}/triple_overlaps.DUP.validated.bed',
        f'{outputdir}/upset_plot.DEL.png',
        f'{outputdir}/upset_plot.DUP.png',
        f'{outputdir}/validated_overlap_stats.txt',
        f'{outputdir}/validated_upset_plot.DEL.png',
        f'{outputdir}/validated_upset_plot.DUP.png',
        f'{outputdir}/validated_percent_upset_plot.DEL.png',
        f'{outputdir}/validated_percent_upset_plot.DUP.png'

rule sort:
    input:
        cnv_calls = cnv_calls,
        SV_calls = SV_calls
    output:
        cnv_calls_sorted = f'{outputdir}/all_calls.sorted.bed',
        SV_calls_sorted = f'{outputdir}/all_SVs.sorted.bed'
    params:
        outdir = outputdir
    shell:
        """
        mkdir -p {params.outdir}
        bedtools sort -i {input.cnv_calls} > {output.cnv_calls_sorted}
        bedtools sort -i {input.SV_calls} > {output.SV_calls_sorted}.tmp
        # strip chr from SVs
        cat {output.SV_calls_sorted}.tmp | sed 's/^chr//' > {output.SV_calls_sorted}
        """

rule get_double_overlaps:
    input:
        cnv_calls = f'{outputdir}/all_calls.sorted.bed',
    output:
        double_overlaps = f'{outputdir}/double_overlaps.bed.txt'
    params:
        outdir = outputdir
    shell:
        """
        # overlap cnvs with cnvs with 60% reciprocal overlap
        bedtools intersect -a {input.cnv_calls} -b {input.cnv_calls} -wao -f .6 -r | cut -f 1-12 > {params.outdir}/all_cnv_x_cnv.bed.tmp
        cat {params.outdir}/all_cnv_x_cnv.bed.tmp | python Scripts/remove_call_to_call_overlaps.py > {params.outdir}/all_cnv_x_cnv.bed
        cat {params.outdir}/all_cnv_x_cnv.bed | python Scripts/make_sample_specific_and_caller_disjoint.py > {output.double_overlaps}        
        """

rule split_doubles:
    input:
        double_overlaps = f'{outputdir}/double_overlaps.bed.txt'
    output:
        split_doubles = f'{outputdir}/split_doubles.bed'
    params:
        outdir = outputdir
    shell:
        """
        cat {input.double_overlaps} | python Scripts/split_doubles.py > {output.split_doubles}
        """

rule overlap_doubles_and_sv:
    input:
        split_doubles = f'{outputdir}/split_doubles.bed',
        SV_calls = f'{outputdir}/all_SVs.sorted.bed'
    output:
        double_overlaps_and_sv = f'{outputdir}/double_overlaps_and_sv.bed'
    params:
        outdir = outputdir
    shell:
        """
        bedtools intersect -a {input.split_doubles} -b {input.SV_calls} -wao -f .6 | cut -f 1-12 | python Scripts/make_sample_specific.py | python Scripts/make_type_specific.py > {output.double_overlaps_and_sv}.tmp
        cat {output.double_overlaps_and_sv}.tmp | cut -f 1-6 | sort | uniq > {output.double_overlaps_and_sv}
        """

def get_calls_from_bed(filename):
    cnvs = set()
    for line in open(filename):
        cnvs.add(line.strip())
    return cnvs

rule generate_stats_about_double_overlaps:
    input:
        double_overlaps_and_sv = f'{outputdir}/double_overlaps_and_sv.bed',
        non_validated_ol = f'{outputdir}/split_doubles.bed'
    output:
        stats = f'{outputdir}/double_stats.txt'
    params:
        outdir = outputdir
    run:
        non_validated_doubles = get_calls_from_bed(input.non_validated_ol)
        with open(output[0],'w') as out:
            out.write('caller\tcnv_type\percent_validated\n')
            for caller in ['GATK','CNVkit','Savvy']:
                for cnv_type in ['DEL','DUP']:
                    count = 0
                    for line in open(input.double_overlaps_and_sv):
                        if cnv_type in line and caller in line:
                            count += 1
                    specific_non_validated_doubles = [x for x in non_validated_doubles if caller in x and cnv_type in x]
                    print(cnv_type,caller,count, len(specific_non_validated_doubles))
                    print(caller,cnv_type,count/len(specific_non_validated_doubles),sep='\t')
                    out.write('\t'.join([caller,cnv_type,str(count/len(specific_non_validated_doubles))])+'\n')

rule overlap_singles_and_sv:
    input:
        cnv_calls_sorted = f'{outputdir}/all_calls.sorted.bed',
        SV_calls_sorted = f'{outputdir}/all_SVs.sorted.bed',
    output:
        f'{outputdir}/single_overlaps_and_sv.bed'
    params:
        percentage = percentage_param
    params:
        outdir = outputdir
    shell:
        """
        bedtools intersect -a {input.cnv_calls_sorted} -b {input.SV_calls_sorted} -wao -f {params.percentage} | cut -f 1-12 | python Scripts/make_sample_specific.py | python Scripts/make_type_specific.py | cut -f1,2,3,4,5,6 | sort | uniq > {output}
        """

rule generate_stats_about_single_overlaps:
    input:
        valid_singles = f'{outputdir}/single_overlaps_and_sv.bed',
        all_doubles = f'{outputdir}/split_doubles.bed',
        all_singles = f'{outputdir}/all_calls.sorted.bed'
    output:
        stats = f'{outputdir}/single_stats.txt'
    params:
        outdir = outputdir
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
        cnv_calls_sorted = f'{outputdir}/all_calls.sorted.bed',
        validated_cnvs = f'{outputdir}/single_overlaps_and_sv.bed'

    output:
        f'{outputdir}/triple_overlaps.DEL.validated.bed',
        f'{outputdir}/triple_overlaps.DUP.validated.bed',
        f'{outputdir}/triple_overlaps.DUP.bed',
        f'{outputdir}/triple_overlaps.DEL.bed',
    params:
        outdir = outputdir
    shell:
        """
            grep GATK {params.outdir}/all_calls.sorted.bed > {params.outdir}/all_calls.GATK.sorted.bed
            grep CNVkit {params.outdir}/all_calls.sorted.bed > {params.outdir}/all_calls.CNVkit.sorted.bed
            grep Savvy {params.outdir}/all_calls.sorted.bed > {params.outdir}/all_calls.Savvy.sorted.bed

            grep DEL {params.outdir}/all_calls.GATK.sorted.bed > {params.outdir}/all_calls.GATK.DEL.sorted.bed
            grep DUP {params.outdir}/all_calls.GATK.sorted.bed > {params.outdir}/all_calls.GATK.DUP.sorted.bed
            grep DEL {params.outdir}/all_calls.CNVkit.sorted.bed > {params.outdir}/all_calls.CNVkit.DEL.sorted.bed
            grep DUP {params.outdir}/all_calls.CNVkit.sorted.bed > {params.outdir}/all_calls.CNVkit.DUP.sorted.bed
            grep DEL {params.outdir}/all_calls.Savvy.sorted.bed > {params.outdir}/all_calls.Savvy.DEL.sorted.bed
            grep DUP {params.outdir}/all_calls.Savvy.sorted.bed > {params.outdir}/all_calls.Savvy.DUP.sorted.bed

            bedtools intersect -a {params.outdir}/all_calls.GATK.DEL.sorted.bed -b {params.outdir}/all_calls.Savvy.DEL.sorted.bed -f .6 -r -wao | cut -f 1-12 | python Scripts/make_sample_specific.py > {params.outdir}/all_calls.GATK.DEL.x_Savvy.DEL.bed
            bedtools intersect -a {params.outdir}/all_calls.GATK.DEL.sorted.bed -b {params.outdir}/all_calls.CNVkit.DEL.sorted.bed -f .6 -r -wao | cut -f 1-12 | python Scripts/make_sample_specific.py > {params.outdir}/all_calls.GATK.DEL.x_CNVkit.DEL.bed
            bedtools intersect -a {params.outdir}/all_calls.Savvy.DEL.sorted.bed -b {params.outdir}/all_calls.CNVkit.DEL.sorted.bed -f .6 -r -wao | cut -f 1-12 | python Scripts/make_sample_specific.py > {params.outdir}/all_calls.Savvy.DEL.x_CNVkit.DEL.bed

            bedtools intersect -a {params.outdir}/all_calls.GATK.DUP.sorted.bed -b {params.outdir}/all_calls.Savvy.DUP.sorted.bed -f .6 -r -wao | cut -f 1-12 | python Scripts/make_sample_specific.py > {params.outdir}/all_calls.GATK.DUP.x_Savvy.DUP.bed
            bedtools intersect -a {params.outdir}/all_calls.GATK.DUP.sorted.bed -b {params.outdir}/all_calls.CNVkit.DUP.sorted.bed -f .6 -r -wao | cut -f 1-12 | python Scripts/make_sample_specific.py > {params.outdir}/all_calls.GATK.DUP.x_CNVkit.DUP.bed
            bedtools intersect -a {params.outdir}/all_calls.Savvy.DUP.sorted.bed -b {params.outdir}/all_calls.CNVkit.DUP.sorted.bed -f .6 -r -wao | cut -f 1-12 | python Scripts/make_sample_specific.py > {params.outdir}/all_calls.Savvy.DUP.x_CNVkit.DUP.bed

            python Scripts/get_triples.py -a {params.outdir}/all_calls.GATK.DEL.x_Savvy.DEL.bed -b {params.outdir}/all_calls.GATK.DEL.x_CNVkit.DEL.bed -c {params.outdir}/all_calls.Savvy.DEL.x_CNVkit.DEL.bed -o {params.outdir}/triple_overlaps.DEL.bed
            python Scripts/get_triples.py -a {params.outdir}/all_calls.GATK.DUP.x_Savvy.DUP.bed -b {params.outdir}/all_calls.GATK.DUP.x_CNVkit.DUP.bed -c {params.outdir}/all_calls.Savvy.DUP.x_CNVkit.DUP.bed -o {params.outdir}/triple_overlaps.DUP.bed

            # intersect triples with real calls
            bedtools intersect -a {params.outdir}/triple_overlaps.DEL.bed -b {input.validated_cnvs} -wao -f .99 -r | cut -f 1-12 | python Scripts/make_sample_specific.py | python Scripts/make_caller_specific.py | cut -f 1-6 | sort | uniq > {params.outdir}/triple_overlaps.DEL.validated.bed
            bedtools intersect -a {params.outdir}/triple_overlaps.DUP.bed -b {input.validated_cnvs} -wao -f .99 -r | cut -f 1-12 | python Scripts/make_sample_specific.py | python Scripts/make_caller_specific.py | cut -f 1-6 | sort | uniq > {params.outdir}/triple_overlaps.DUP.validated.bed
        """


rule generate_stats_about_overlaps:
    input:
        all_overlaps = f'{outputdir}/all_calls.sorted.bed',
        doubles = f'{outputdir}/double_overlaps.bed.txt',
        triples_del = f'{outputdir}/triple_overlaps.DEL.bed',
        triples_dup = f'{outputdir}/triple_overlaps.DUP.bed'
    output:
        stats = f'{outputdir}/overlap_stats.txt'
    params:
        outdir = outputdir
    shell:
        """
        # remove doubles that overlap a DEL and a DUP
        grep -v 'DUP' {input.doubles} > {params.outdir}/double_overlaps.DEL.bed
        grep -v 'DEL' {input.doubles} > {params.outdir}/double_overlaps.DUP.bed
        cat {params.outdir}/double_overlaps.DEL.bed > {params.outdir}/double_overlaps.type_specific.bed
        cat {params.outdir}/double_overlaps.DUP.bed >> {params.outdir}/double_overlaps.type_specific.bed

        # combine triples
        cat {input.triples_del} > {params.outdir}/triple_overlaps.type_specific.bed
        cat {input.triples_dup} >> {params.outdir}/triple_overlaps.type_specific.bed

        python Scripts/generate_stats.py --singles {input.all_overlaps} --double {params.outdir}/double_overlaps.type_specific.bed --triples {params.outdir}/triple_overlaps.type_specific.bed --output {output.stats}
        """

rule upset_plots:
    input:
        stats = f'{outputdir}/overlap_stats.txt'
    output:
        upset_plot_del = f'{outputdir}/upset_plot.DEL.png',
        upset_plot_dup = f'{outputdir}/upset_plot.DUP.png'
    params:
        outdir = outputdir
    shell:
        """
        # replace CNVkit with CNVKit in stats file
        sed -i 's/CNVkit/CNVKit/g' {input.stats}

        grep DEL {input.stats} | cut -f2 > {params.outdir}/DEL.stats.txt
        python Scripts/upset_plot.py -i {params.outdir}/DEL.stats.txt -l Deletion -o {output.upset_plot_del}
        grep DUP {input.stats} | cut -f2 > {params.outdir}/DUP.stats.txt
        python Scripts/upset_plot.py -i {params.outdir}/DUP.stats.txt -l Duplication -o {output.upset_plot_dup}
        """

rule generate_stats_about_validated_overlaps:
    input:
        real_calls= f'{outputdir}/single_overlaps_and_sv.bed',

        triples_del = f'{outputdir}/triple_overlaps.DEL.validated.bed',
        triples_dup = f'{outputdir}/triple_overlaps.DUP.validated.bed'
    output:
        stats = f'{outputdir}/validated_overlap_stats.txt'
    params:
        outdir = outputdir
    shell:
        """
        # get the double overlap real calls
        bedtools intersect -a {input.real_calls} -b {input.real_calls} -wao -f .6 -r | cut -f 1-12 | python Scripts/make_sample_specific.py | python Scripts/make_type_specific.py > {params.outdir}/validated_double_overlaps.bed.txt


        
        # remove doubles that overlap a DEL and a DUP
        grep -v 'DUP' {params.outdir}/validated_double_overlaps.bed.txt > {params.outdir}/validated_double_overlaps.DEL.bed
        grep -v 'DEL' {params.outdir}/validated_double_overlaps.bed.txt > {params.outdir}/validated_double_overlaps.DUP.bed
        cat {params.outdir}/validated_double_overlaps.DEL.bed > {params.outdir}/validated_double_overlaps.type_specific.bed
        cat {params.outdir}/validated_double_overlaps.DUP.bed >> {params.outdir}/validated_double_overlaps.type_specific.bed

        # combine triples
        cat {input.triples_del} > {params.outdir}/validated_triple_overlaps.type_specific.bed
        cat {input.triples_dup} >> {params.outdir}/validated_triple_overlaps.type_specific.bed

        python Scripts/generate_stats.py --singles {input.real_calls} --double {params.outdir}/validated_double_overlaps.type_specific.bed --triples {params.outdir}/triple_overlaps.type_specific.bed --output {output.stats}
        """

rule calc_stats_as_portions:
    input:
        stats = f'{outputdir}/overlap_stats.txt',
        valid_stats = f'{outputdir}/validated_overlap_stats.txt'
    output:
        stats = f'{outputdir}/validated_percent_overlap_stats.txt'
    params:
        outdir = outputdir
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
        stats = f'{outputdir}/validated_overlap_stats.txt'
    output:
        upset_plot_del = f'{outputdir}/validated_upset_plot.DEL.png',
        upset_plot_dup = f'{outputdir}/validated_upset_plot.DUP.png'
    params:
        outdir = outputdir,
        percent = percentage_param
    shell:
        """
        # replace CNVkit with CNVKit in stats file
        sed -i 's/CNVkit/CNVKit/g' {input.stats}

        grep DEL {input.stats} | cut -f2 > {params.outdir}/DEL.stats.txt
        python Scripts/upset_plot.py -i {params.outdir}/DEL.stats.txt -l 'Deletion {params.percent}' -o {output.upset_plot_del}
        grep DUP {input.stats} | cut -f2 > {params.outdir}/DUP.stats.txt
        python Scripts/upset_plot.py -i {params.outdir}/DUP.stats.txt -l 'Duplication {params.percent}' -o {output.upset_plot_dup}
        """

rule upset_plots_valdiated_percents:
    input:
        stats = f'{outputdir}/validated_percent_overlap_stats.txt'
    output:
        upset_plot_del = f'{outputdir}/validated_percent_upset_plot.DEL.png',
        upset_plot_dup = f'{outputdir}/validated_percent_upset_plot.DUP.png'
    params:
        outdir = outputdir,
        percent = percentage_param
    shell:
        """
        # replace CNVkit with CNVKit in stats file
        sed -i 's/CNVkit/CNVKit/g' {input.stats}

        grep DEL {input.stats} | cut -f2 > {params.outdir}/DEL.stats.percent.txt
        python Scripts/upset_plot.py -i {params.outdir}/DEL.stats.percent.txt -l 'Deletion {params.percent}' -o {output.upset_plot_del} --logscale 0 --percent
        grep DUP {input.stats} | cut -f2 > {params.outdir}/DUP.stats.percent.txt
        python Scripts/upset_plot.py -i {params.outdir}/DUP.stats.percent.txt -l 'Duplication {params.percent}' -o {output.upset_plot_dup} --logscale 0 --percent
        """