cnv_calls = 'Data/all_calls.bed'
SV_calls = 'Data/all_SVs.bed'

rule all:
    input:
        'debugwork/double_overlaps.bed.txt',
        'debugwork/split_doubles.bed',
        'debugwork/double_overlaps_and_sv.bed',
        'debugwork/double_stats.txt',
        'debugwork/single_stats.txt',
        'debugwork/triple_overlaps.DEL.validated.bed',
        'debugwork/triple_overlaps.DUP.validated.bed'
        

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

rule overlap_doubles_and_sv:
    input:
        split_doubles = 'debugwork/split_doubles.bed',
        SV_calls = 'debugwork/all_SVs.sorted.bed'
    output:
        double_overlaps_and_sv = 'debugwork/double_overlaps_and_sv.bed'
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
        double_overlaps_and_sv = 'debugwork/double_overlaps_and_sv.bed',
        non_validated_ol = 'debugwork/split_doubles.bed'
    output:
        stats = 'debugwork/double_stats.txt'
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
        cnv_calls_sorted = 'debugwork/all_calls.sorted.bed',
        SV_calls_sorted = 'debugwork/all_SVs.sorted.bed',
    output:
        'debugwork/single_overlaps_and_sv.bed'
    shell:
        """
        bedtools intersect -a {input.cnv_calls_sorted} -b {input.SV_calls_sorted} -wao -f .9 | cut -f 1-12 | python Scripts/make_sample_specific.py | python Scripts/make_type_specific.py tmp.bed | cut -f1,2,3,4,5,6 | sort | uniq > {output}
        """

rule generate_stats_about_single_overlaps:
    input:
        valid_singles = 'debugwork/single_overlaps_and_sv.bed',
        all_doubles = 'debugwork/split_doubles.bed',
        all_singles = 'debugwork/all_calls.sorted.bed'
    output:
        stats = 'debugwork/single_stats.txt'
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
        validated_cnvs = 'debugwork/single_overlaps_and_sv.bed'

    output:
        'debugwork/triple_overlaps.DEL.validated.bed',
        'debugwork/triple_overlaps.DUP.validated.bed'
    shell:
        """
            grep GATK debugwork/all_calls.sorted.bed > debugwork/all_calls.GATK.sorted.bed
            grep CNVkit debugwork/all_calls.sorted.bed > debugwork/all_calls.CNVkit.sorted.bed
            grep Savvy debugwork/all_calls.sorted.bed > debugwork/all_calls.Savvy.sorted.bed

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

            # intersect triples with real calls
            bedtools intersect -a debugwork/triple_overlaps.DEL.bed -b {input.validated_cnvs} -wao -f .99 -r | cut -f 1-12 | python Scripts/make_sample_specific.py | python Scripts/make_caller_specific.py | cut -f 1-6 | sort | uniq > debugwork/triple_overlaps.DEL.validated.bed
            bedtools intersect -a debugwork/triple_overlaps.DUP.bed -b {input.validated_cnvs} -wao -f .99 -r | cut -f 1-12 | python Scripts/make_sample_specific.py | python Scripts/make_caller_specific.py | cut -f 1-6 | sort | uniq > debugwork/triple_overlaps.DUP.validated.bed
        """