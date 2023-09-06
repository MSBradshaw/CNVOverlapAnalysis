
# The list of calls made in external repos
cnv_calls = 'Data/all_calls.bed'
SV_calls = 'Data/all_SVs.bed'
PERCENTS = ['0.10','0.15','0.20','0.25','0.30','0.35','0.40','0.45','0.50','0.55','0.60','0.65','0.70','0.75','0.80','0.85','0.90','0.95']
bedpePercents = [str(i) for i in range(10,100,5)]

# list of numbers from 00 to 99
SPLITS = list( str(i) if i > 9 else '0' + str(i) for i in range(100))
CNV_TYPES = ['DEL', 'DUP']
CALLERS = ['Savvy','CNVkit','GATK']

rule all:
    input:
        expand('work/CNVTypes/cnv_calls.{cnv_type}.bed', cnv_type=CNV_TYPES),
        expand('work/CNVTypes/svs_calls.{cnv_type}.bed', cnv_type=CNV_TYPES),
        expand('work/CNVTypes/cnv_calls.{cnv_type}.bed.gz', cnv_type=CNV_TYPES),
        expand('work/CNVTypes/cnv_calls.{cnv_type}.bed.gz.tbi', cnv_type=CNV_TYPES),
        expand('work/CNVTypes/svs_calls.{cnv_type}.bed.gz', cnv_type=CNV_TYPES),
        expand('work/CNVTypes/svs_calls.{cnv_type}.bed.gz.tbi', cnv_type=CNV_TYPES),
        expand('work/Recip/reciprocal_overlap.{percent}.{cnv_type}.txt', percent=PERCENTS, cnv_type=CNV_TYPES),
        expand('work/SPLITS/cnv_split_{split}.{cnv_type}.bed.gz', split=SPLITS, cnv_type=CNV_TYPES),
        expand('work/SPLITS/cnv_split_{split}.{cnv_type}.bed', split=SPLITS, cnv_type=CNV_TYPES),
        expand('work/PercentSplits/call_level_reciprocal_overlap.{percent}.{split}.{cnv_type}.txt', split=SPLITS, percent=PERCENTS, cnv_type=CNV_TYPES),
        expand('work/PercentSplits/agg_call_level_reciprocal_overlap.{percent}.{cnv_type}.txt', percent=PERCENTS, cnv_type=CNV_TYPES),
        'Figures/reciprocal_overlap_4_panel.png',

        # --------------------- BedPE Outputs ---------------------
        expand('work/BedPEOverlap/cnv_calls.{cnv_type}.0p{percent}.bedpe.overlap.txt', percent=bedpePercents, cnv_type=CNV_TYPES),
        expand('work/BedPESPLITS/cnv_split_{split}.{cnv_type}.bedpe', split=SPLITS, cnv_type=CNV_TYPES),
        expand('work/BedPEPercentSplits/call_level_reciprocal_overlap.0p{percent}.{split}.{cnv_type}.txt', split=SPLITS, percent=bedpePercents, cnv_type=CNV_TYPES),
        'Figures/breakpoint_bedpe_overlap_4_panel.png',

        # --------------------- SLOP Outputs ---------------------
        expand('work/BedPESlopPercentSplits/call_level_slop_breakpoint_overlap.0p{percent}.{split}.{cnv_type}.txt', split=SPLITS, percent=bedpePercents, cnv_type=CNV_TYPES),
        'work/BedPESlopPercentSplits/agg_call_level_slop_breakpoint_overlap.all.txt',
        'Figures/slop_breakpoint_bedpe_overlap_2_panel.png',

        # --------------------- Ven Diagrams Bed ---------------------
        expand('work/CallerSpecificCNVTypes/cnv_calls.{cnv_type}.{caller}.bed',cnv_type=CNV_TYPES, caller=CALLERS),
        expand('work/CallerSpecificBedPECNVTypes/cnv_calls.{cnv_type}.{caller}.bedpe',cnv_type=CNV_TYPES, caller=CALLERS),
        expand('work/CallerSpecificOverlapsBed/venn.{cnv_type}.{caller}.bed.txt',cnv_type=CNV_TYPES, caller=CALLERS),
        expand('work/CallerSpecificOverTripleOverlap/venn.{cnv_type}.triple.bed.txt',cnv_type=CNV_TYPES),
        expand('work/VenDiagramResults/overlap_numbers.{cnv_type}.txt',cnv_type=CNV_TYPES),
        expand('work/VenDiagramResults/overlap_numbers.{cnv_type}.sample_specific.txt',cnv_type=CNV_TYPES),
        expand('work/VenDiagramResults/overlapping_calls.{cnv_type}.sample_specific.triples.bed.txt',cnv_type=CNV_TYPES),

        
        # --------------------- Ven Diagrams BedPE ---------------------
        # expand('work/CallerSpecificOverlapsBedPE/venn.{cnv_type}.{caller}.bedpe.txt',cnv_type=CNV_TYPES, caller=CALLERS),
        expand('work/CallerSpecificTripleOverlapBedPE/venn.{cnv_type}.triple.bedpe.txt',cnv_type=CNV_TYPES),
        expand('work/VenDiagramResultsBedPE/overlap_numbers.{cnv_type}.txt',cnv_type=CNV_TYPES),
        expand('work/VenDiagramResultsBedPE/overlap_numbers.{cnv_type}.sample_specific.txt',cnv_type=CNV_TYPES),
        expand('work/VenDiagramResultsBoth/venn.bed_and_bedpe.{cnv_type}.txt',cnv_type=CNV_TYPES),
        expand('work/VenDiagramResultsBoth/venn.bed_and_bedpe.{cnv_type}.sample_specific.txt',cnv_type=CNV_TYPES),

        # --------------------- Triple calls ---------------------
        expand('work/VenDiagramResultsBoth/venn.bed_and_bedpe.{cnv_type}.sample_specific.triples.bed.txt', cnv_type=CNV_TYPES),
        expand('work/VenDiagramResults/overlap_numbers.{cnv_type}.sample_specific.triples.bed.txt', cnv_type=CNV_TYPES),
        expand('work/VenDiagramResultsBedPE/overlap_numbers.{cnv_type}.sample_specific.triples.bed.txt', cnv_type=CNV_TYPES),
        'work/AllSampleSpecificTripleCalls/venn.bed_and_bedpe.all.sample_specific.triples.bed.txt',

        # --------------------- CNV validation with SVs ---------------------
        expand('work/SV_Validation/real_calls.{cnv_type}.bed',cnv_type=CNV_TYPES),
        expand('work/SV_Validation/real_calls.{cnv_type}.report.txt',cnv_type=CNV_TYPES)
        

# split cnv_calls by DEL and DUP
rule split_cnv_calls:
    input:
        cnv_calls=cnv_calls,
        sv_bed = SV_calls
    params:
        cnv_type=lambda wildcards, output: output[0].split('.')[1]
    output:
        'work/CNVTypes/cnv_calls.{cnv_type}.bed',
        'work/CNVTypes/svs_calls.{cnv_type}.bed',
    shell:
        """
        grep "{params.cnv_type}" {input.cnv_calls} > {output[0]}.tmp
        grep "{params.cnv_type}" {input.sv_bed} > {output[1]}.tmp

        # sort the cnv calls
        bedtools sort -i {output[0]}.tmp > {output[0]}
        bedtools sort -i {output[1]}.tmp > {output[1]}
        """

rule sort_bgzip_tabix:
    input:
        cnv_bed = 'work/CNVTypes/cnv_calls.{cnv_type}.bed',
        sv_bed = 'work/CNVTypes/svs_calls.{cnv_type}.bed'
    output:
        cnv_gz='work/CNVTypes/cnv_calls.{cnv_type}.bed.gz',
        cnv_tbi='work/CNVTypes/cnv_calls.{cnv_type}.bed.gz.tbi',
        sv_gz='work/CNVTypes/svs_calls.{cnv_type}.bed.gz',
        sv_tbi='work/CNVTypes/svs_calls.{cnv_type}.bed.gz.tbi',
    shell:
        """
        bedtools sort -i {input.cnv_bed} | bgzip -c > {output.cnv_gz}
        tabix -p bed {output.cnv_gz}
        
        bedtools sort -i {input.sv_bed} | bgzip -c > {output.sv_gz}
        tabix -p bed {output.sv_gz}
        """

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


"""
----- Breakpoint Overlap -----
"""

rule convert_to_bedpe:
    input:
        cnv='work/CNVTypes/cnv_calls.{cnv_type}.bed',
        sv='work/CNVTypes/svs_calls.{cnv_type}.bed'
    output:
        cnv='work/BedPECNVTypes/cnv_calls.{cnv_type}.bedpe',
        sv='work/BedPECNVTypes/svs_calls.{cnv_type}.bedpe'
    shell:
        """
        mkdir -p work/BedPECNVTypes
        # columns for bedpe
        # chrom1 start1 end1 chrom2 start2 end2 name score strand1 strand2 extras...
        # make sure strand 1 and 2 are filled in with .
        cat {input.cnv} | awk '{{print $1"\t"$2-1"\t"$2+1"\t"$1"\t"$3-1"\t"$3+1"\t.\t.\t.\t.\t"$4"\t"$5"\t"$6"\t"}}' > {output.cnv}
        cat {input.sv} | awk '{{print $1"\t"$2-1"\t"$2+1"\t"$1"\t"$3-1"\t"$3+1"\t.\t.\t.\t.\t"$4"\t"$5"\t"$6"\t"}}' > {output.sv}
        """

rule bedpe_breakpoint_percent_overlap_all_v_all:
    input:
        cnv='work/BedPECNVTypes/cnv_calls.{cnv_type}.bedpe'
    params:
        percent=lambda wildcards, output: output[0].split('0p')[1].split('.')[0]
    output:
        'work/BedPEOverlap/cnv_calls.{cnv_type}.0p{percent}.bedpe.overlap.txt'
    shell:
        """
        mkdir -p work/BedPEOverlap
        bedtools pairtopair -f 0.{params.percent} -a {input.cnv} -b {input.cnv} -type both | wc -l > {output}.tmp
        # annotate with percent and cnv type
        cat {output}.tmp | awk -v percent={wildcards.percent} '{{print $0"\t0."percent}}' | awk -v cnv_type={wildcards.cnv_type} '{{print $0"\t"cnv_type}}' > {output}
        """


# split the bedpe cnv calls into 100 files with ~620 calls each
rule sample_level_bedpe_overlap_split:
    input:
        cnv='work/BedPECNVTypes/cnv_calls.{cnv_type}.bedpe',
    output:
        cnv=expand('work/BedPESPLITS/cnv_split_{split}.{cnv_type}.bedpe', split=SPLITS, cnv_type='{cnv_type}')
    shell:
        """
        mkdir -p work/BedPESPLITS
        split  --additional-suffix=.{wildcards.cnv_type}.bedpe -d -n l/100 {input} work/BedPESPLITS/cnv_split_
        """

rule sample_level_bedpe_breakpoint_overlap:
    input:
        cnv_bed='work/BedPESPLITS/cnv_split_{split}.{cnv_type}.bedpe',
        all_cnv='work/BedPECNVTypes/cnv_calls.{cnv_type}.bedpe',
    params:
        percent=lambda wildcards, output: output[0].split('0p')[1].split('.')[0]
    output:
        'work/BedPEPercentSplits/call_level_reciprocal_overlap.0p{percent}.{split}.{cnv_type}.txt'
    threads: 2
    log: 
        'logs/BedPEPercentSplits/call_level_reciprocal_overlap.0p{percent}.{split}.{cnv_type}.log'
    shell:
        """
        echo 'Percent: {wildcards.percent}'
        mkdir -p work/BedPEPercentSplits
        mkdir -p logs/BedPEPercentSplits

        # remove output, if it exists already
        rm -f {output} 2>> {log}
        # for each line in input.cnv_bed
        while IFS= read -r line; do
            echo '-' >> {log}
            echo "$line"  >> {log}
            echo '-' >> {log}
            echo "$line" > work/BedPEPercentSplits/call_level_reciprocal_overlap_{wildcards.split}_0p{wildcards.percent}_{wildcards.cnv_type}.tmp.bedpe
            cat work/BedPEPercentSplits/call_level_reciprocal_overlap_{wildcards.split}_0p{wildcards.percent}_{wildcards.cnv_type}.tmp.bedpe >> {log}
            bedtools pairtopair -a work/BedPEPercentSplits/call_level_reciprocal_overlap_{wildcards.split}_0p{wildcards.percent}_{wildcards.cnv_type}.tmp.bedpe -b {input.all_cnv} -type both -f 0.{params.percent} | wc -l >> {output} 2>> {log}
            echo 'line done' >> {log}
        done < {input.cnv_bed}
        """

rule sample_level_bedpe_slop_overlap:
    input:
        cnv_bed='work/BedPESPLITS/cnv_split_{split}.{cnv_type}.bedpe',
        all_cnv='work/BedPECNVTypes/cnv_calls.{cnv_type}.bedpe',
    params:
        percent=lambda wildcards, output: output[0].split('0p')[1].split('.')[0]
    output:
        'work/BedPESlopPercentSplits/call_level_slop_breakpoint_overlap.0p{percent}.{split}.{cnv_type}.txt'
    threads: 2
    log: 
        'logs/BedPESlopPercentSplits/call_level_slop_breakpoint_overlap.0p{percent}.{split}.{cnv_type}.log'
    shell:
        """
        echo 'Percent: {wildcards.percent}'
        mkdir -p work/BedPESlopPercentSplits
        mkdir -p logs/BedPESlopPercentSplits

        # remove output, if it exists already
        rm -f {output} 2>> {log}
        # for each line in input.cnv_bed
        while IFS= read -r line; do
            echo '-' >> {log}
            echo "$line"  >> {log}
            echo '-' >> {log}
            echo "$line" > work/BedPESlopPercentSplits/call_level_slop_breakpoint_overlap_{wildcards.split}_0p{wildcards.percent}_{wildcards.cnv_type}.tmp.bedpe
            slop=$(python Scripts/get_slop_from_percent.py work/BedPESlopPercentSplits/call_level_slop_breakpoint_overlap_{wildcards.split}_0p{wildcards.percent}_{wildcards.cnv_type}.tmp.bedpe 0.{params.percent})
            cat work/BedPESlopPercentSplits/call_level_slop_breakpoint_overlap_{wildcards.split}_0p{wildcards.percent}_{wildcards.cnv_type}.tmp.bedpe >> {log}
            bedtools pairtopair -a work/BedPESlopPercentSplits/call_level_slop_breakpoint_overlap_{wildcards.split}_0p{wildcards.percent}_{wildcards.cnv_type}.tmp.bedpe -b {input.all_cnv} -type both -slop $slop | wc -l | awk -v percent={params.percent} '{{print $0"\t0."percent}}' | awk -v cnvtype={wildcards.cnv_type} '{{print $0"\t"cnvtype}}' | awk -v slop=$slop '{{print $0"\t"slop}}' >> {output} 2>> {log}
            echo 'line done' >> {log}
        done < {input.cnv_bed}
        rm work/BedPESlopPercentSplits/call_level_slop_breakpoint_overlap_{wildcards.split}_0p{wildcards.percent}_{wildcards.cnv_type}.tmp.bedpe
        """

# aggregate all slop
rule agg_sample_level_bedpe_slop_overlap:
    input:
        expand('work/BedPESlopPercentSplits/call_level_slop_breakpoint_overlap.0p{percent}.{split}.{cnv_type}.txt', split=SPLITS, percent=bedpePercents, cnv_type=CNV_TYPES)
    output:
        'work/BedPESlopPercentSplits/agg_call_level_slop_breakpoint_overlap.all.txt'
    shell:
        """
        mkdir -p work/BedPESlopPercentSplits
        cat {input} > {output}
        """

# plot slop
rule plot_slop_overlap:
    input:
        'work/BedPESlopPercentSplits/agg_call_level_slop_breakpoint_overlap.all.txt'
    output:
        'Figures/slop_breakpoint_bedpe_overlap_2_panel.png'
    shell:
        """
        python Scripts/plot_slop_overlap.py {input} {output} slop breakpoint
        """

# aggregate sample_level_reciprocal_overlap by percentage, cat all files for a % into one file
rule annotate_and_aggregate_sample_level_bedpe_overlap:
    input:
        expand('work/BedPEPercentSplits/call_level_reciprocal_overlap.0p{percent}.{split}.{cnv_type}.txt', split=SPLITS, percent='{percent}', cnv_type='{cnv_type}')
    output:
        'work/BedPEPercentSplitsAnnotated/agg_call_level_reciprocal_overlap.{percent}.{cnv_type}.txt'
    shell:
        """
        mkdir -p work/BedPEPercentSplitsAnnotated
        cat {input} > {output}.tmp
        # annotate with percent and cnv type
        cat {output}.tmp | awk -v percent={wildcards.percent} '{{print $0"\t0."percent}}' | awk -v cnv_type={wildcards.cnv_type} '{{print $0"\t"cnv_type}}' > {output}
        """

rule plot_bedpe_overlap:
    # take all the agg_call_level_reciprocal_overlap files and plot them
    input:
        single_files = expand('work/BedPEPercentSplitsAnnotated/agg_call_level_reciprocal_overlap.{percent}.{cnv_type}.txt', percent=bedpePercents, cnv_type=CNV_TYPES),
        all_by_all_files = expand('work/BedPEOverlap/cnv_calls.{cnv_type}.0p{percent}.bedpe.overlap.txt', percent=bedpePercents, cnv_type=CNV_TYPES)
    output:
        'Figures/breakpoint_bedpe_overlap_4_panel.png'
    shell:
        """
        mkdir -p Figures/
        cat {input.single_files} > work/agg_call_level_breakpoint_bedpe_overlap.all.txt
        cat {input.all_by_all_files} > work/agg_breakpoint_bedpe_overlap.all.txt
        python Scripts/plot_overlap.py work/agg_call_level_breakpoint_bedpe_overlap.all.txt work/agg_breakpoint_bedpe_overlap.all.txt {output} breakpoint
        """

# ----------------------------- Ven Diagrams ----------------------------- #

# # divide calls based on caller
rule divide_by_caller:
    input:
        bed='work/CNVTypes/cnv_calls.{cnv_type}.bed',
        bedpe='work/BedPECNVTypes/cnv_calls.{cnv_type}.bedpe'
    params:
        caller=lambda wildcards, output: output[0].split('.')[2]
    output:
        bed='work/CallerSpecificCNVTypes/cnv_calls.{cnv_type}.{caller}.bed',
        bedpe='work/CallerSpecificBedPECNVTypes/cnv_calls.{cnv_type}.{caller}.bedpe'
    shell:
        """
        mkdir -p work/CallerSpecificCNVTypes
        mkdir -p work/CallerSpecificBedPECNVTypes

        grep {params.caller} {input.bed} > {output.bed}
        grep {params.caller} {input.bedpe} > {output.bedpe}
        """

# # do the intersections for bed
rule intersect_callers_bed:
    input:
        bed='work/CallerSpecificCNVTypes/cnv_calls.{cnv_type}.{caller}.bed',
        all_beds=expand('work/CallerSpecificCNVTypes/cnv_calls.{cnv_type}.{caller}.bed', caller=CALLERS, cnv_type=CNV_TYPES)
    params:
        savvy = 'Savvy',
        cnvkit = 'CNVkit',
        gatk = 'GATK'
    output:
        'work/CallerSpecificOverlapsBed/venn.{cnv_type}.{caller}.bed.txt'
    shell:
        """
        mkdir -p work/CallerSpecificOverlapsBed
        rm -f {output}

        # ---- savvy ---- #
        STR='work/CallerSpecificCNVTypes/cnv_calls.{wildcards.cnv_type}.Savvy.bed'
        SUB='{wildcards.caller}'
        if [[ "$STR" != *"$SUB"* ]]; then
            bedtools intersect -a {input.bed} -b $STR -wo -f 0.60 -r >> {output}
        fi
        
        # ---- cnvkit ---- #
        STR='work/CallerSpecificCNVTypes/cnv_calls.{wildcards.cnv_type}.CNVkit.bed'
        SUB='{wildcards.caller}'
        if [[ "$STR" != *"$SUB"* ]]; then
            bedtools intersect -a {input.bed} -b $STR -wo -f 0.60 -r >> {output}
        fi

        # ---- gatk ---- #
        STR='work/CallerSpecificCNVTypes/cnv_calls.{wildcards.cnv_type}.GATK.bed'
        SUB='{wildcards.caller}'
        if [[ "$STR" != *"$SUB"* ]]; then
            bedtools intersect -a {input.bed} -b $STR -wo -f 0.60 -r >> {output}
        fi
        """

# find the triple overlaps
rule triple_overlap_bed:
    input:
        expand('work/CallerSpecificOverlapsBed/venn.{cnv_type}.{caller}.bed.txt', caller=CALLERS, cnv_type='{cnv_type}')
    output:
        'work/CallerSpecificOverTripleOverlap/venn.{cnv_type}.triple.bed.txt'
    shell:
        """
        mkdir -p work/CallerSpecificOverTripleOverlap
        python Scripts/triple_overlap.py -s work/CallerSpecificOverlapsBed/venn.{wildcards.cnv_type}.Savvy.bed.txt \
                -g work/CallerSpecificOverlapsBed/venn.{wildcards.cnv_type}.GATK.bed.txt \
                -c work/CallerSpecificOverlapsBed/venn.{wildcards.cnv_type}.CNVkit.bed.txt > {output}
        """

# report the number of single, double and triple overlaps
rule report_overlap_numbers:
    input:
        triple='work/CallerSpecificOverTripleOverlap/venn.{cnv_type}.triple.bed.txt',
        double=expand('work/CallerSpecificOverlapsBed/venn.{cnv_type}.{caller}.bed.txt', caller=CALLERS, cnv_type='{cnv_type}'),
        single=expand('work/CallerSpecificCNVTypes/cnv_calls.{cnv_type}.{caller}.bed', caller=CALLERS, cnv_type='{cnv_type}')
    output:
        non_specific='work/VenDiagramResults/overlap_numbers.{cnv_type}.txt',
        sample_specific='work/VenDiagramResults/overlap_numbers.{cnv_type}.sample_specific.txt',
        sample_specific_triples='work/VenDiagramResults/overlap_numbers.{cnv_type}.sample_specific.triples.bed.txt',
        sample_specific_triple_calls='work/VenDiagramResults/overlapping_calls.{cnv_type}.sample_specific.triples.bed.txt'
    shell:
        """
        mkdir -p work/VenDiagramResults
        # number of calls in each caller
        python Scripts/count_overlaps.py -s work/CallerSpecificCNVTypes/cnv_calls.{wildcards.cnv_type}.Savvy.bed \
            -c work/CallerSpecificCNVTypes/cnv_calls.{wildcards.cnv_type}.CNVkit.bed \
            -g work/CallerSpecificCNVTypes/cnv_calls.{wildcards.cnv_type}.GATK.bed \
            --sx work/CallerSpecificOverlapsBed/venn.{wildcards.cnv_type}.Savvy.bed.txt \
            --cx work/CallerSpecificOverlapsBed/venn.{wildcards.cnv_type}.CNVkit.bed.txt \
            --gx work/CallerSpecificOverlapsBed/venn.{wildcards.cnv_type}.GATK.bed.txt \
            -t {input.triple} \
            -o {output.non_specific}

        # force calls to be from the same sample
         python Scripts/count_overlaps.py -s work/CallerSpecificCNVTypes/cnv_calls.{wildcards.cnv_type}.Savvy.bed \
            -c work/CallerSpecificCNVTypes/cnv_calls.{wildcards.cnv_type}.CNVkit.bed \
            -g work/CallerSpecificCNVTypes/cnv_calls.{wildcards.cnv_type}.GATK.bed \
            --sx work/CallerSpecificOverlapsBed/venn.{wildcards.cnv_type}.Savvy.bed.txt \
            --cx work/CallerSpecificOverlapsBed/venn.{wildcards.cnv_type}.CNVkit.bed.txt \
            --gx work/CallerSpecificOverlapsBed/venn.{wildcards.cnv_type}.GATK.bed.txt \
            -t {input.triple} \
            -o {output.sample_specific} --fs --to {output.sample_specific_triples} \
            --sample_specific_triples_out {output.sample_specific_triple_calls}

        """

# ------------------------ BEDPE Ven Diagrams ------------------------ #
# do the intersections for bed
rule intersect_callers_bedpe:
    input:
        bedpe='work/CallerSpecificBedPECNVTypes/cnv_calls.{cnv_type}.{caller}.bedpe',
        all_bedpes=expand('work/CallerSpecificBedPECNVTypes/cnv_calls.{cnv_type}.{caller}.bedpe', caller=CALLERS, cnv_type=CNV_TYPES)
    params:
        savvy = 'Savvy',
        cnvkit = 'CNVkit',
        gatk = 'GATK'
    output:
        'work/CallerSpecificOverlapsBedPE/venn.{cnv_type}.{caller}.bedpe.txt'
    shell:
        """
        mkdir -p work/CallerSpecificOverlapsBedPE
        rm -f {output}

        # ---- savvy ---- #
        STR='work/CallerSpecificBedPECNVTypes/cnv_calls.{wildcards.cnv_type}.Savvy.bedpe'
        SUB='{wildcards.caller}'
        if [[ "$STR" != *"$SUB"* ]]; then
            while IFS= read -r line; do
                echo "$line" > work/CallerSpecificOverlapsBedPE/tmp.{wildcards.cnv_type}.{wildcards.caller}.bedpe
                slop=$(python Scripts/get_slop_from_percent.py work/CallerSpecificOverlapsBedPE/tmp.{wildcards.cnv_type}.{wildcards.caller}.bedpe 0.50)
                bedtools pairtopair -a work/CallerSpecificOverlapsBedPE/tmp.{wildcards.cnv_type}.{wildcards.caller}.bedpe -b $STR -slop $slop >> {output}
            done < {input.bedpe}
        fi
        
        # ---- cnvkit ---- #
        STR='work/CallerSpecificBedPECNVTypes/cnv_calls.{wildcards.cnv_type}.CNVkit.bedpe'
        SUB='{wildcards.caller}'
        if [[ "$STR" != *"$SUB"* ]]; then
            while IFS= read -r line; do
                echo "$line" > work/CallerSpecificOverlapsBedPE/tmp.{wildcards.cnv_type}.{wildcards.caller}.bedpe
                slop=$(python Scripts/get_slop_from_percent.py work/CallerSpecificOverlapsBedPE/tmp.{wildcards.cnv_type}.{wildcards.caller}.bedpe 0.50)
                bedtools pairtopair -a work/CallerSpecificOverlapsBedPE/tmp.{wildcards.cnv_type}.{wildcards.caller}.bedpe -b $STR -slop $slop >> {output}
            done < {input.bedpe}
        fi

        # ---- gatk ---- #
        STR='work/CallerSpecificBedPECNVTypes/cnv_calls.{wildcards.cnv_type}.GATK.bedpe'
        SUB='{wildcards.caller}'
        if [[ "$STR" != *"$SUB"* ]]; then
            while IFS= read -r line; do
                echo "$line" > work/CallerSpecificOverlapsBedPE/tmp.{wildcards.cnv_type}.{wildcards.caller}.bedpe
                slop=$(python Scripts/get_slop_from_percent.py work/CallerSpecificOverlapsBedPE/tmp.{wildcards.cnv_type}.{wildcards.caller}.bedpe 0.50)
                bedtools pairtopair -a work/CallerSpecificOverlapsBedPE/tmp.{wildcards.cnv_type}.{wildcards.caller}.bedpe -b $STR -slop $slop >> {output}
            done < {input.bedpe}
        fi
        """

rule bedpe_triple_overlap:
    input:
        expand('work/CallerSpecificOverlapsBedPE/venn.{cnv_type}.{caller}.bedpe.txt', caller=CALLERS, cnv_type='{cnv_type}')
    output:
        'work/CallerSpecificTripleOverlapBedPE/venn.{cnv_type}.triple.bedpe.txt'
    shell:
        """
        mkdir -p work/CallerSpecificTripleOverlapBedPE

        python Scripts/triple_overlap.py -s work/CallerSpecificOverlapsBedPE/venn.{wildcards.cnv_type}.Savvy.bedpe.txt \
                    -g work/CallerSpecificOverlapsBedPE/venn.{wildcards.cnv_type}.GATK.bedpe.txt \
                    -c work/CallerSpecificOverlapsBedPE/venn.{wildcards.cnv_type}.CNVkit.bedpe.txt -p > {output}
        """

# report the number of single, double and triple overlaps
rule report_overlap_numbers_bedpe:
    input:
        triple='work/CallerSpecificTripleOverlapBedPE/venn.{cnv_type}.triple.bedpe.txt',
        double=expand('work/CallerSpecificOverlapsBedPE/venn.{cnv_type}.{caller}.bedpe.txt', caller=CALLERS, cnv_type='{cnv_type}'),
        single=expand('work/CallerSpecificBedPECNVTypes/cnv_calls.{cnv_type}.{caller}.bedpe', caller=CALLERS, cnv_type='{cnv_type}')
    output:
        non_specific='work/VenDiagramResultsBedPE/overlap_numbers.{cnv_type}.txt',
        sample_specific='work/VenDiagramResultsBedPE/overlap_numbers.{cnv_type}.sample_specific.txt',
        sample_specific_triples='work/VenDiagramResultsBedPE/overlap_numbers.{cnv_type}.sample_specific.triples.bed.txt'
    shell:
        """
        mkdir -p work/VenDiagramResultsBedPE
        # number of calls in each caller
        python Scripts/count_overlaps.py -s work/CallerSpecificBedPECNVTypes/cnv_calls.{wildcards.cnv_type}.Savvy.bedpe \
            -c work/CallerSpecificBedPECNVTypes/cnv_calls.{wildcards.cnv_type}.CNVkit.bedpe \
            -g work/CallerSpecificBedPECNVTypes/cnv_calls.{wildcards.cnv_type}.GATK.bedpe \
            --sx work/CallerSpecificOverlapsBedPE/venn.{wildcards.cnv_type}.Savvy.bedpe.txt \
            --cx work/CallerSpecificOverlapsBedPE/venn.{wildcards.cnv_type}.CNVkit.bedpe.txt \
            --gx work/CallerSpecificOverlapsBedPE/venn.{wildcards.cnv_type}.GATK.bedpe.txt \
            -t {input.triple} \
            -o {output.non_specific} -p
        
        python Scripts/count_overlaps.py -s work/CallerSpecificBedPECNVTypes/cnv_calls.{wildcards.cnv_type}.Savvy.bedpe \
            -c work/CallerSpecificBedPECNVTypes/cnv_calls.{wildcards.cnv_type}.CNVkit.bedpe \
            -g work/CallerSpecificBedPECNVTypes/cnv_calls.{wildcards.cnv_type}.GATK.bedpe \
            --sx work/CallerSpecificOverlapsBedPE/venn.{wildcards.cnv_type}.Savvy.bedpe.txt \
            --cx work/CallerSpecificOverlapsBedPE/venn.{wildcards.cnv_type}.CNVkit.bedpe.txt \
            --gx work/CallerSpecificOverlapsBedPE/venn.{wildcards.cnv_type}.GATK.bedpe.txt \
            -t {input.triple} \
            -o {output.sample_specific} --fs -p --to {output.sample_specific_triples}
        """


# ------------------------------ BED + BEDPE Overlaps Venn Diagrams ------------------------------ #

rule report_all_bed_and_bedpe_overlap:
    input:
        bedpe_triple='work/CallerSpecificTripleOverlapBedPE/venn.{cnv_type}.triple.bedpe.txt',
        bedpe_double=expand('work/CallerSpecificOverlapsBedPE/venn.{cnv_type}.{caller}.bedpe.txt', caller=CALLERS, cnv_type='{cnv_type}'),
        triple='work/CallerSpecificOverTripleOverlap/venn.{cnv_type}.triple.bed.txt',
        double=expand('work/CallerSpecificOverlapsBed/venn.{cnv_type}.{caller}.bed.txt', caller=CALLERS, cnv_type='{cnv_type}'),
    output:
        non_specific='work/VenDiagramResultsBoth/venn.bed_and_bedpe.{cnv_type}.txt',
        sample_specific='work/VenDiagramResultsBoth/venn.bed_and_bedpe.{cnv_type}.sample_specific.txt',
        sample_specific_triples='work/VenDiagramResultsBoth/venn.bed_and_bedpe.{cnv_type}.sample_specific.triples.bed.txt',
    shell:
        """
        mkdir -p work/VenDiagramResultsBoth
        python Scripts/get_bed_bedpe_overlaps.py \
            --sb work/CallerSpecificOverlapsBed/venn.{wildcards.cnv_type}.Savvy.bed.txt \
            --sp work/CallerSpecificOverlapsBedPE/venn.{wildcards.cnv_type}.Savvy.bedpe.txt \
            --cb work/CallerSpecificOverlapsBed/venn.{wildcards.cnv_type}.CNVkit.bed.txt \
            --cp  work/CallerSpecificOverlapsBedPE/venn.{wildcards.cnv_type}.CNVkit.bedpe.txt \
            --gb work/CallerSpecificOverlapsBed/venn.{wildcards.cnv_type}.GATK.bed.txt \
            --gp  work/CallerSpecificOverlapsBedPE/venn.{wildcards.cnv_type}.GATK.bedpe.txt \
            --tb {input.triple} \
            --tp {input.bedpe_triple} \
            --sx work/VenDiagramResultsBoth/tmp.{wildcards.cnv_type}.Savvy.bed.txt \
            --cx work/VenDiagramResultsBoth/tmp.{wildcards.cnv_type}.CNVkit.bed.txt \
            --gx work/VenDiagramResultsBoth/tmp.{wildcards.cnv_type}.GATK.bed.txt \
            --tx work/VenDiagramResultsBoth/tmp.{wildcards.cnv_type}.triple.bed.txt
        
        # overlap again but force the samples to be the same
        python Scripts/count_overlaps.py -s work/CallerSpecificCNVTypes/cnv_calls.{wildcards.cnv_type}.Savvy.bed \
            -c work/CallerSpecificCNVTypes/cnv_calls.{wildcards.cnv_type}.CNVkit.bed \
            -g work/CallerSpecificCNVTypes/cnv_calls.{wildcards.cnv_type}.GATK.bed \
            --sx work/VenDiagramResultsBoth/tmp.{wildcards.cnv_type}.Savvy.bed.txt \
            --cx work/VenDiagramResultsBoth/tmp.{wildcards.cnv_type}.CNVkit.bed.txt \
            --gx work/VenDiagramResultsBoth/tmp.{wildcards.cnv_type}.GATK.bed.txt \
            -t  work/VenDiagramResultsBoth/tmp.{wildcards.cnv_type}.triple.bed.txt \
            -o {output.non_specific}
        
        # overlap again but force the samples to be the same
        python Scripts/count_overlaps.py -s work/CallerSpecificCNVTypes/cnv_calls.{wildcards.cnv_type}.Savvy.bed \
            -c work/CallerSpecificCNVTypes/cnv_calls.{wildcards.cnv_type}.CNVkit.bed \
            -g work/CallerSpecificCNVTypes/cnv_calls.{wildcards.cnv_type}.GATK.bed \
            --sx work/VenDiagramResultsBoth/tmp.{wildcards.cnv_type}.Savvy.bed.txt \
            --cx work/VenDiagramResultsBoth/tmp.{wildcards.cnv_type}.CNVkit.bed.txt \
            --gx work/VenDiagramResultsBoth/tmp.{wildcards.cnv_type}.GATK.bed.txt \
            -t  work/VenDiagramResultsBoth/tmp.{wildcards.cnv_type}.triple.bed.txt \
            -o {output.sample_specific} --fs --to {output.sample_specific_triples}

        """

rule merge_all_triple_calls:
    input:
        expand('work/VenDiagramResultsBoth/venn.bed_and_bedpe.{cnv_type}.sample_specific.triples.bed.txt', cnv_type=CNV_TYPES),
        expand('work/VenDiagramResults/overlap_numbers.{cnv_type}.sample_specific.triples.bed.txt', cnv_type=CNV_TYPES),
        expand('work/VenDiagramResultsBedPE/overlap_numbers.{cnv_type}.sample_specific.triples.bed.txt', cnv_type=CNV_TYPES)
    output:
        'work/AllSampleSpecificTripleCalls/venn.bed_and_bedpe.all.sample_specific.triples.bed.txt'
    shell:
        """
        cat {input}  | sort | uniq > {output}
        """

rule find_real_calls:
    input:
        svs='work/CNVTypes/svs_calls.{cnv_type}.bed.gz',
        cnvs='work/CNVTypes/cnv_calls.{cnv_type}.bed.gz'
    output:
        bed='work/SV_Validation/real_calls.{cnv_type}.bed',
        txt='work/SV_Validation/real_calls.{cnv_type}.report.txt'
    shell:
        """
        mkdir -p work/SV_Validation
        # this requires the SV to overlap 90% of the CNV
        bedtools intersect -a {input.svs} -b {input.cnvs} -wb -F .9 > {output.bed}

        # count the number of GATK calls in CNVs
        echo 'GATK' > {output.txt}
        echo -n 'CNVs: ' >> {output.txt}
        zcat {input.cnvs} | grep GATK | wc -l >> {output.txt}
        echo -n 'Validated with SV: ' >> {output.txt}
        zcat {output.bed} | grep GATK | wc -l >> {output.txt}

        # count the number of CNVkit calls in CNVs
        echo 'CNVkit' >> {output.txt}
        echo -n 'CNVs: ' >> {output.txt}
        zcat {input.cnvs} | grep CNVkit | wc -l >> {output.txt}
        echo -n 'Validated with SV: ' >> {output.txt}
        zcat {output.bed} | grep CNVkit | wc -l >> {output.txt}

        # count the number of Savvy calls in CNVs
        echo 'Savvy' >> {output.txt}
        echo -n 'CNVs: ' >> {output.txt}
        zcat {input.cnvs} | grep Savvy | wc -l >> {output.txt}
        echo -n 'Validated with SV: ' >> {output.txt}
        zcat {output.bed} | grep Savvy | wc -l >> {output.txt}
        """


    
