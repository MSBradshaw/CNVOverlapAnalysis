
# The list of calls made in external repos
cnv_calls = 'Data/all_calls.bed'
SV_calls = 'Data/all_SVs.bed'
PERCENTS = ['0.10','0.15','0.20','0.25','0.30','0.35','0.40','0.45','0.50','0.55','0.60','0.65','0.70','0.75','0.80','0.85','0.90','0.95']
bedpePercents = [str(i) for i in range(10,100,5)]
ENDS=['savvy','cnvkit','gatk','triple','savvy_gatk','savvy_cnvkit','cnvkit_gatk']

# list of numbers from 00 to 99
SPLITS = list( str(i) if i > 9 else '0' + str(i) for i in range(100))
CNV_TYPES = ['DEL', 'DUP']
CALLERS = ['Savvy','CNVkit','GATK']

rule all:
    input:
        'Figures/reciprocal_overlap_4_panel.png',

# split cnv_calls by DEL and DUP
rule split_cnv_calls:
    input:
        cnv_calls=cnv_calls,
    params:
        cnv_type=lambda wildcards, output: output[0].split('.')[1]
    output:
        'work/CNVTypes/cnv_calls.{cnv_type}.bed',
    shell:
        """
        grep "{params.cnv_type}" {input.cnv_calls} > {output[0]}.tmp

        # sort the cnv calls
        bedtools sort -i {output[0]}.tmp > {output[0]}
        bedtools sort -i {output[1]}.tmp > {output[1]}
        # remove prefix chr from chromosome names in {output[1]}
        cat {output[1]} | sed -e "s/^chr//" > {output[1]}.tmp
        cat {output[1]}.tmp | sed -e "s/GATK Cohort Mode/1KG_SV/" > {output[1]}
        mv {output[1]}.tmp {output[1]}
        """

rule sort_bgzip_tabix:
    input:
        cnv_bed = 'work/CNVTypes/cnv_calls.{cnv_type}.bed',
    output:
        cnv_gz='work/CNVTypes/cnv_calls.{cnv_type}.bed.gz',
        cnv_tbi='work/CNVTypes/cnv_calls.{cnv_type}.bed.gz.tbi',
    shell:
        """
        bedtools sort -i {input.cnv_bed} | bgzip -c > {output.cnv_gz}
        tabix -p bed {output.cnv_gz}
        """
# keep
rule reciprocal_overlap:
    input:
        cnv_gz='work/CNVTypes/cnv_calls.{cnv_type}.bed.gz',
        cnv_tbi='work/CNVTypes/cnv_calls.{cnv_type}.bed.gz.tbi',
    params:
        percent=lambda wildcards, output: '0.' + output[0].split('.')[2]
    output:
        'work/Recip/reciprocal_overlap.{percent}.{cnv_type}.txt'
    shell:
        """
        mkdir -p work/Recip/
        bedtools intersect -a {input.cnv_gz} -b {input.cnv_gz} -wa -wb -f {params.percent} | wc -l > {output}
        """

# split the cnv calls into 100 files with ~620 calls each
rule sample_level_reciprocal_overlap_split:
    input:
        cnv_bed = 'work/CNVTypes/cnv_calls.{cnv_type}.bed',
    output:
        expand('work/SPLITS/cnv_split_{split}.{cnv_type}.bed.gz', split=SPLITS, cnv_type='{cnv_type}'),
        expand('work/SPLITS/cnv_split_{split}.{cnv_type}.bed', split=SPLITS, cnv_type='{cnv_type}')
    shell:
        """
        mkdir -p work/SPLITS
        # sort the cnv calls
        bedtools sort -i {input.cnv_bed} > work/cnv_sorted_{wildcards.cnv_type}.bed
        split  --additional-suffix=.{wildcards.cnv_type}.bed -d -n l/100 work/cnv_sorted_{wildcards.cnv_type}.bed work/SPLITS/cnv_split_
        files=$(ls work/SPLITS/cnv_split_*.{wildcards.cnv_type}.bed | grep -v gz)
        for f in $files;
        do
            echo $f
            cat $f  | bgzip -c > "$f.gz";
            tabix -p bed "$f.gz";
        done
        """

#  intersect each line from work/SPLITS/cnv_split_{split}.bed with work/SPLITS/cnv_split_{split}.bed.gz and count the number of lines
rule sample_level_reciprocal_overlap:
    input:
        cnv_bed='work/SPLITS/cnv_split_{split}.{cnv_type}.bed',
        cnv_gz='work/SPLITS/cnv_split_{split}.{cnv_type}.bed.gz',
        all_cnv_bed = 'work/CNVTypes/cnv_calls.{cnv_type}.bed.gz',
    params:
        percent=lambda wildcards, output: '0.' + output[0].split('.')[2]
    output:
        'work/PercentSplits/call_level_reciprocal_overlap.{percent}.{split}.{cnv_type}.txt'
    threads: 2
    log: 
        'logs/call_level_reciprocal_overlap.{percent}.{split}.{cnv_type}.log'
    shell:
        """
        mkdir -p work/PercentSplits
        # remove output, if it exists already
        rm -f {output} 2>> {log}
        # for each line in input.cnv_bed
        while IFS= read -r line; do
            echo '-' >> {log}
            echo "$line"  >> {log}
            echo '-' >> {log}
            echo "$line" > work/SPLITS/cnv_split_{wildcards.split}_{wildcards.percent}_{wildcards.cnv_type}.tmp.bed
            cat work/SPLITS/cnv_split_{wildcards.split}_{wildcards.percent}_{wildcards.cnv_type}.tmp.bed >> {log}
            # sleep 1s
            bgzip -f work/SPLITS/cnv_split_{wildcards.split}_{wildcards.percent}_{wildcards.cnv_type}.tmp.bed  2>> {log}
            tabix -f work/SPLITS/cnv_split_{wildcards.split}_{wildcards.percent}_{wildcards.cnv_type}.tmp.bed.gz 2>> {log}
            bedtools intersect -a work/SPLITS/cnv_split_{wildcards.split}_{wildcards.percent}_{wildcards.cnv_type}.tmp.bed.gz -b {input.all_cnv_bed} -wa -wb -r -f {params.percent} | wc -l >> {output} 2>> {log}
            echo 'line done' >> {log}
        done < {input.cnv_bed}
        """

# keep
# aggregate sample_level_reciprocal_overlap by percentage, cat all files for a % into one file
rule aggregate_sample_level_reciprocal_overlap:
    input:
        expand('work/PercentSplits/call_level_reciprocal_overlap.{percent}.{split}.{cnv_type}.txt', split=SPLITS, percent='{percent}', cnv_type='{cnv_type}')
    output:
        'work/PercentSplits/agg_call_level_reciprocal_overlap.{percent}.{cnv_type}.txt'
    shell:
        """
        cat {input} > {output}
        """
# keep
# add a second column with the % and a third with CNV type to the agg_call_level_reciprocal_overlap files so they are easier to plot later
rule annotate_agg_sample_level_reciprocal_overlap:
    input:
        'work/PercentSplits/agg_call_level_reciprocal_overlap.{percent}.{cnv_type}.txt'
    output:
        'work/AnnotatedPercentSplits/agg_call_level_reciprocal_overlap.{percent}.{cnv_type}.annotated.txt'
    shell:
        """
        mkdir -p work/AnnotatedPercentSplits
        cat {input} | awk -v percent={wildcards.percent} '{{print $0"\t"percent}}' | awk -v cnv_type={wildcards.cnv_type} '{{print $0"\t"cnv_type}}' > {output}
        """
# keep
rule annotate_all_by_all_reciprocal_overlap:
    input:
        'work/Recip/reciprocal_overlap.{percent}.{cnv_type}.txt'
    output:
        'work/AnnotatedRecip/reciprocal_overlap.{percent}.{cnv_type}.annotated.txt'
    shell:
        """
        mkdir -p work/AnnotatedRecip
        cat {input} | awk -v percent={wildcards.percent} '{{print $0"\t"percent}}' | awk -v cnv_type={wildcards.cnv_type} '{{print $0"\t"cnv_type}}' > {output}
        """
# keep
rule plot_reciprocal_overlap:
    # take all the agg_call_level_reciprocal_overlap files and plot them
    input:
        single_files = expand('work/AnnotatedPercentSplits/agg_call_level_reciprocal_overlap.{percent}.{cnv_type}.annotated.txt', percent=PERCENTS, cnv_type=CNV_TYPES),
        all_by_all_files = expand('work/AnnotatedRecip/reciprocal_overlap.{percent}.{cnv_type}.annotated.txt', percent=PERCENTS, cnv_type=CNV_TYPES)
    output:
        'Figures/reciprocal_overlap_4_panel.png'
    shell:
        """
        mkdir -p Figures/
        cat {input.single_files} > work/agg_call_level_reciprocal_overlap.all.txt
        cat {input.all_by_all_files} > work/agg_reciprocal_overlap.all.txt
        python Scripts/plot_overlap.py work/agg_call_level_reciprocal_overlap.all.txt work/agg_reciprocal_overlap.all.txt {output} reciprocal
        """

