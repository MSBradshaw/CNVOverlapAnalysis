# The list of calls made in external repos
cnv_calls = 'Data/all_calls.chco.bed'
PERCENTS = ['0.10','0.15','0.20','0.25','0.30','0.35','0.40','0.45','0.50','0.55','0.60','0.65','0.70','0.75','0.80','0.85','0.90','0.95']
bedpePercents = [str(i) for i in range(10,100,5)]
ENDS=['savvy','cnvkit','gatk','triple','savvy_gatk','savvy_cnvkit','cnvkit_gatk']

# list of numbers from 00 to 99
SPLITS = list( str(i) if i > 9 else '0' + str(i) for i in range(100))
CNV_TYPES = ['DEL', 'DUP']
CALLERS = ['Savvy','CNVkit','GATK']

rule all:
    input:
        expand('workchco/CNVTypes/cnv_calls.{cnv_type}.bed', cnv_type=CNV_TYPES),
        expand('workchco/CNVTypes/cnv_calls.{cnv_type}.bed.gz', cnv_type=CNV_TYPES),
        expand('workchco/CNVTypes/cnv_calls.{cnv_type}.bed.gz.tbi', cnv_type=CNV_TYPES),
        expand('workchco/Recip/reciprocal_overlap.{percent}.{cnv_type}.txt', percent=PERCENTS, cnv_type=CNV_TYPES),
        expand('workchco/SPLITS/cnv_split_{split}.{cnv_type}.bed.gz', split=SPLITS, cnv_type=CNV_TYPES),
        expand('workchco/SPLITS/cnv_split_{split}.{cnv_type}.bed', split=SPLITS, cnv_type=CNV_TYPES),
        expand('workchco/PercentSplits/call_level_reciprocal_overlap.{percent}.{split}.{cnv_type}.txt', split=SPLITS, percent=PERCENTS, cnv_type=CNV_TYPES),
        expand('workchco/PercentSplits/agg_call_level_reciprocal_overlap.{percent}.{cnv_type}.txt', percent=PERCENTS, cnv_type=CNV_TYPES),
        'CHCOFigures/reciprocal_overlap_4_panel.png',

        # --------------------- BedPE Outputs ---------------------
        expand('workchco/BedPEOverlap/cnv_calls.{cnv_type}.0p{percent}.bedpe.overlap.txt', percent=bedpePercents, cnv_type=CNV_TYPES),
        expand('workchco/BedPESPLITS/cnv_split_{split}.{cnv_type}.bedpe', split=SPLITS, cnv_type=CNV_TYPES),
        expand('workchco/BedPEPercentSplits/call_level_reciprocal_overlap.0p{percent}.{split}.{cnv_type}.txt', split=SPLITS, percent=bedpePercents, cnv_type=CNV_TYPES),
        'CHCOFigures/breakpoint_bedpe_overlap_4_panel.png',

        # --------------------- SLOP Outputs ---------------------
        expand('workchco/BedPESlopPercentSplits/call_level_slop_breakpoint_overlap.0p{percent}.{split}.{cnv_type}.txt', split=SPLITS, percent=bedpePercents, cnv_type=CNV_TYPES),
        'workchco/BedPESlopPercentSplits/agg_call_level_slop_breakpoint_overlap.all.txt',
        'CHCOFigures/slop_breakpoint_bedpe_overlap_2_panel.png',

        # --------------------- Ven Diagrams Bed ---------------------
        expand('workchco/CallerSpecificCNVTypes/cnv_calls.{cnv_type}.{caller}.bed',cnv_type=CNV_TYPES, caller=CALLERS),
        expand('workchco/CallerSpecificBedPECNVTypes/cnv_calls.{cnv_type}.{caller}.bedpe',cnv_type=CNV_TYPES, caller=CALLERS),
        expand('workchco/CallerSpecificOverlapsBed/venn.{cnv_type}.{caller}.bed.txt',cnv_type=CNV_TYPES, caller=CALLERS),
        expand('workchco/CallerSpecificOverTripleOverlap/venn.{cnv_type}.triple.bed.txt',cnv_type=CNV_TYPES),
        expand('workchco/VenDiagramResults/overlap_numbers.{cnv_type}.txt',cnv_type=CNV_TYPES),
        expand('workchco/VenDiagramResults/overlap_numbers.{cnv_type}.sample_specific.txt',cnv_type=CNV_TYPES),
        expand('workchco/VenDiagramResults/overlapping_calls.{cnv_type}.sample_specific.triples.bed.txt',cnv_type=CNV_TYPES),

        
        # --------------------- Ven Diagrams BedPE ---------------------
        # expand('workchco/CallerSpecificOverlapsBedPE/venn.{cnv_type}.{caller}.bedpe.txt',cnv_type=CNV_TYPES, caller=CALLERS),
        expand('workchco/CallerSpecificTripleOverlapBedPE/venn.{cnv_type}.triple.bedpe.txt',cnv_type=CNV_TYPES),
        expand('workchco/VenDiagramResultsBedPE/overlap_numbers.{cnv_type}.txt',cnv_type=CNV_TYPES),
        expand('workchco/VenDiagramResultsBedPE/overlap_numbers.{cnv_type}.sample_specific.txt',cnv_type=CNV_TYPES),
        # expand('workchco/VenDiagramResultsBoth/venn.bed_and_bedpe.{cnv_type}.txt',cnv_type=CNV_TYPES),
        # expand('workchco/VenDiagramResultsBoth/venn.bed_and_bedpe.{cnv_type}.sample_specific.txt',cnv_type=CNV_TYPES),

        # --------------------- Triple calls ---------------------
        # expand('workchco/VenDiagramResultsBoth/venn.bed_and_bedpe.{cnv_type}.sample_specific.triples.bed.txt', cnv_type=CNV_TYPES),
        expand('workchco/VenDiagramResults/overlap_numbers.{cnv_type}.sample_specific.triples.bed.txt', cnv_type=CNV_TYPES),
        expand('workchco/VenDiagramResultsBedPE/overlap_numbers.{cnv_type}.sample_specific.triples.bed.txt', cnv_type=CNV_TYPES),
        'workchco/AllSampleSpecificTripleCalls/venn.bed_and_bedpe.all.sample_specific.triples.bed.txt',

        # --------------------- Additional plotting ---------------------
        'CHCOFigures/reciprocal_and_slop_breakpoint_boxplots.png',
        'CHCOFigures/size_distribution.png',
        'CHCOFigures/calls_per_sample.png',
        # expand('CHCOFigures/upsetplot.{cnv_type}.png',cnv_type=CNV_TYPES),
        # expand('CHCOFigures/upsetplot.slop.{cnv_type}.png',cnv_type=CNV_TYPES),
        expand('workchco/RecipCategoriesBeds/{cnv_type}{end}.bed',cnv_type=CNV_TYPES,end=['-savvy','-cnvkit','-gatk','-triple','-savvy_gatk','-savvy_cnvkit','-cnvkit_gatk']),
        expand('workchco/RecipCategoriesBedPEs/{cnv_type}{end}.bed',cnv_type=CNV_TYPES,end=['-savvy','-cnvkit','-gatk','-triple','-savvy_gatk','-savvy_cnvkit','-cnvkit_gatk']),
        

rule read_data:
    input:
        savvy='CHCO_CNV_calling/SavvycnvResults/cnv_list_reformatted.bed',
        cnvkit='CHCO_CNV_calling/CNVkitResults/aggregate_calls.bed',
        gatk='Data/GATK_gCNV_all_nonref_calls.bed'
    output:
        'Data/all_calls.chco.bed'
    # params:
        # callers = lambda output, wildcards: output.split('_')[0],
        # cnv_types = lambda output, wildcards: output.split('_')[1]
    shell:
        """
        mkdir -p work/CHCOCalls
        # re order the last 2 columns of cnvkit
        cat {input.cnvkit} | awk -v OFS='\t' '{{print $1,$2,$3,$5,$4}}' > work/CHCOCalls/cnvkit.tmp.bed

        cat work/CHCOCalls/cnvkit.tmp.bed | grep 'DEL' > work/CHCOCalls/CNVKit_DEL.tmp.bed || true
        cat work/CHCOCalls/cnvkit.tmp.bed | grep 'DUP' > work/CHCOCalls/CNVKit_DUP.tmp.bed || true
        cat {input.savvy} | grep 'DEL' > work/CHCOCalls/Savvy_DEL.tmp.bed || true
        cat {input.savvy} | grep 'DUP' > work/CHCOCalls/Savvy_DUP.tmp.bed || true

        cat {input.gatk} | awk -v OFS='\t' '{{print $0,"GATK"}}' | grep 'DEL' > work/CHCOCalls/GATK_DEL.tmp.bed || true
        cat {input.gatk} | awk -v OFS='\t' '{{print $0,"GATK"}}' | grep 'DUP' > work/CHCOCalls/GATK_DUP.tmp.bed || true

        # add a column to both files that has the caller name in it
        cat work/CHCOCalls/CNVKit_DEL.tmp.bed | awk -v OFS='\t' '{{print $0,"CNVkit"}}' > work/CHCOCalls/CNVKit_DEL.tmp2.bed
        cat work/CHCOCalls/CNVKit_DUP.tmp.bed | awk -v OFS='\t' '{{print $0,"CNVkit"}}' > work/CHCOCalls/CNVKit_DUP.tmp2.bed
        cat work/CHCOCalls/Savvy_DEL.tmp.bed | awk -v OFS='\t' '{{print $0,"Savvy"}}' > work/CHCOCalls/Savvy_DEL.tmp2.bed
        cat work/CHCOCalls/Savvy_DUP.tmp.bed | awk -v OFS='\t' '{{print $0,"Savvy"}}' > work/CHCOCalls/Savvy_DUP.tmp2.bed

        bedtools sort -i work/CHCOCalls/CNVKit_DEL.tmp2.bed > work/CHCOCalls/CNVKit_DEL.bed
        bedtools sort -i work/CHCOCalls/CNVKit_DUP.tmp2.bed > work/CHCOCalls/CNVKit_DUP.bed
        bedtools sort -i work/CHCOCalls/Savvy_DEL.tmp2.bed > work/CHCOCalls/Savvy_DEL.bed
        bedtools sort -i work/CHCOCalls/Savvy_DUP.tmp2.bed > work/CHCOCalls/Savvy_DUP.bed
        bedtools sort -i work/CHCOCalls/GATK_DEL.tmp.bed > work/CHCOCalls/GATK_DEL.bed
        bedtools sort -i work/CHCOCalls/GATK_DUP.tmp.bed > work/CHCOCalls/GATK_DUP.bed

    
        cat work/CHCOCalls/CNVKit_DEL.bed > {output}
        cat work/CHCOCalls/CNVKit_DUP.bed >> {output}
        cat work/CHCOCalls/Savvy_DEL.bed >> {output}
        cat work/CHCOCalls/Savvy_DUP.bed >> {output}
        cat work/CHCOCalls/GATK_DEL.bed >> {output}
        cat work/CHCOCalls/GATK_DUP.bed >> {output}
        """

# split cnv_calls by DEL and DUP
rule split_cnv_calls:
    input:
        cnv_calls='Data/all_calls.chco.bed',
    params:
        cnv_type=lambda wildcards, output: output[0].split('.')[1]
    output:
        'workchco/CNVTypes/cnv_calls.{cnv_type}.bed',
    shell:
        """
        grep "{params.cnv_type}" {input.cnv_calls} > {output[0]}.tmp

        # sort the cnv calls
        bedtools sort -i {output[0]}.tmp > {output[0]}
        """

rule sort_bgzip_tabix:
    input:
        cnv_bed = 'workchco/CNVTypes/cnv_calls.{cnv_type}.bed',
    output:
        cnv_gz='workchco/CNVTypes/cnv_calls.{cnv_type}.bed.gz',
        cnv_tbi='workchco/CNVTypes/cnv_calls.{cnv_type}.bed.gz.tbi',

    shell:
        """
        bedtools sort -i {input.cnv_bed} | bgzip -c > {output.cnv_gz}
        tabix -p bed {output.cnv_gz}
        """

rule reciprocal_overlap:
    input:
        cnv_gz='workchco/CNVTypes/cnv_calls.{cnv_type}.bed.gz',
        cnv_tbi='workchco/CNVTypes/cnv_calls.{cnv_type}.bed.gz.tbi',
    params:
        percent=lambda wildcards, output: '0.' + output[0].split('.')[2]
    output:
        'workchco/Recip/reciprocal_overlap.{percent}.{cnv_type}.txt'
    shell:
        """
        mkdir -p workchco/Recip/
        bedtools intersect -a {input.cnv_gz} -b {input.cnv_gz} -wa -wb -f {params.percent} | wc -l > {output}
        """

# split the cnv calls into 100 files with ~620 calls each
rule sample_level_reciprocal_overlap_split:
    input:
        cnv_bed = 'workchco/CNVTypes/cnv_calls.{cnv_type}.bed',
    output:
        expand('workchco/SPLITS/cnv_split_{split}.{cnv_type}.bed.gz', split=SPLITS, cnv_type='{cnv_type}'),
        expand('workchco/SPLITS/cnv_split_{split}.{cnv_type}.bed', split=SPLITS, cnv_type='{cnv_type}')
    shell:
        """
        mkdir -p workchco/SPLITS
        # sort the cnv calls
        bedtools sort -i {input.cnv_bed} > workchco/cnv_sorted_{wildcards.cnv_type}.bed
        split  --additional-suffix=.{wildcards.cnv_type}.bed -d -n l/100 workchco/cnv_sorted_{wildcards.cnv_type}.bed workchco/SPLITS/cnv_split_
        files=$(ls workchco/SPLITS/cnv_split_*.{wildcards.cnv_type}.bed | grep -v gz)
        for f in $files;
        do
            echo $f
            cat $f  | bgzip -c > "$f.gz";
            tabix -p bed "$f.gz";
        done
        """

#  intersect each line from workchco/SPLITS/cnv_split_{split}.bed with workchco/SPLITS/cnv_split_{split}.bed.gz and count the number of lines
rule sample_level_reciprocal_overlap:
    input:
        cnv_bed='workchco/SPLITS/cnv_split_{split}.{cnv_type}.bed',
        cnv_gz='workchco/SPLITS/cnv_split_{split}.{cnv_type}.bed.gz',
        all_cnv_bed = 'workchco/CNVTypes/cnv_calls.{cnv_type}.bed.gz',
    params:
        percent=lambda wildcards, output: '0.' + output[0].split('.')[2]
    output:
        'workchco/PercentSplits/call_level_reciprocal_overlap.{percent}.{split}.{cnv_type}.txt'
    threads: 2
    log: 
        'logs/call_level_reciprocal_overlap.{percent}.{split}.{cnv_type}.log'
    shell:
        """
        mkdir -p workchco/PercentSplits
        # remove output, if it exists already
        rm -f {output} 2>> {log}
        # for each line in input.cnv_bed
        while IFS= read -r line; do
            echo '-' >> {log}
            echo "$line"  >> {log}
            echo '-' >> {log}
            echo "$line" > workchco/SPLITS/cnv_split_{wildcards.split}_{wildcards.percent}_{wildcards.cnv_type}.tmp.bed
            cat workchco/SPLITS/cnv_split_{wildcards.split}_{wildcards.percent}_{wildcards.cnv_type}.tmp.bed >> {log}
            # sleep 1s
            bgzip -f workchco/SPLITS/cnv_split_{wildcards.split}_{wildcards.percent}_{wildcards.cnv_type}.tmp.bed  2>> {log}
            tabix -f workchco/SPLITS/cnv_split_{wildcards.split}_{wildcards.percent}_{wildcards.cnv_type}.tmp.bed.gz 2>> {log}
            bedtools intersect -a workchco/SPLITS/cnv_split_{wildcards.split}_{wildcards.percent}_{wildcards.cnv_type}.tmp.bed.gz -b {input.all_cnv_bed} -wa -wb -r -f {params.percent} | wc -l >> {output} 2>> {log}
            echo 'line done' >> {log}
        done < {input.cnv_bed}
        """

# aggregate sample_level_reciprocal_overlap by percentage, cat all files for a % into one file
rule aggregate_sample_level_reciprocal_overlap:
    input:
        expand('workchco/PercentSplits/call_level_reciprocal_overlap.{percent}.{split}.{cnv_type}.txt', split=SPLITS, percent='{percent}', cnv_type='{cnv_type}')
    output:
        'workchco/PercentSplits/agg_call_level_reciprocal_overlap.{percent}.{cnv_type}.txt'
    shell:
        """
        cat {input} > {output}
        """

# add a second column with the % and a third with CNV type to the agg_call_level_reciprocal_overlap files so they are easier to plot later
rule annotate_agg_sample_level_reciprocal_overlap:
    input:
        'workchco/PercentSplits/agg_call_level_reciprocal_overlap.{percent}.{cnv_type}.txt'
    output:
        'workchco/AnnotatedPercentSplits/agg_call_level_reciprocal_overlap.{percent}.{cnv_type}.annotated.txt'
    shell:
        """
        mkdir -p workchco/AnnotatedPercentSplits
        cat {input} | awk -v percent={wildcards.percent} '{{print $0"\t"percent}}' | awk -v cnv_type={wildcards.cnv_type} '{{print $0"\t"cnv_type}}' > {output}
        """
rule annotate_all_by_all_reciprocal_overlap:
    input:
        'workchco/Recip/reciprocal_overlap.{percent}.{cnv_type}.txt'
    output:
        'workchco/AnnotatedRecip/reciprocal_overlap.{percent}.{cnv_type}.annotated.txt'
    shell:
        """
        mkdir -p workchco/AnnotatedRecip
        cat {input} | awk -v percent={wildcards.percent} '{{print $0"\t"percent}}' | awk -v cnv_type={wildcards.cnv_type} '{{print $0"\t"cnv_type}}' > {output}
        """

rule plot_reciprocal_overlap:
    # take all the agg_call_level_reciprocal_overlap files and plot them
    input:
        single_files = expand('workchco/AnnotatedPercentSplits/agg_call_level_reciprocal_overlap.{percent}.{cnv_type}.annotated.txt', percent=PERCENTS, cnv_type=CNV_TYPES),
        all_by_all_files = expand('workchco/AnnotatedRecip/reciprocal_overlap.{percent}.{cnv_type}.annotated.txt', percent=PERCENTS, cnv_type=CNV_TYPES)
    output:
        'CHCOFigures/reciprocal_overlap_4_panel.png'
    shell:
        """
        mkdir -p CHCOFigures/
        cat {input.single_files} > workchco/agg_call_level_reciprocal_overlap.all.txt
        cat {input.all_by_all_files} > workchco/agg_reciprocal_overlap.all.txt
        python Scripts/plot_overlap.py workchco/agg_call_level_reciprocal_overlap.all.txt workchco/agg_reciprocal_overlap.all.txt {output} reciprocal
        """


"""
----- Breakpoint Overlap -----
"""

rule convert_to_bedpe:
    input:
        cnv='workchco/CNVTypes/cnv_calls.{cnv_type}.bed'
    output:
        cnv='workchco/BedPECNVTypes/cnv_calls.{cnv_type}.bedpe',
    shell:
        """
        mkdir -p workchco/BedPECNVTypes
        # columns for bedpe
        # chrom1 start1 end1 chrom2 start2 end2 name score strand1 strand2 extras...
        # make sure strand 1 and 2 are filled in with .
        cat {input.cnv} | awk '{{print $1"\t"$2-1"\t"$2+1"\t"$1"\t"$3-1"\t"$3+1"\t.\t.\t.\t.\t"$4"\t"$5"\t"$6"\t"}}' > {output.cnv}
        """

rule bedpe_breakpoint_percent_overlap_all_v_all:
    input:
        cnv='workchco/BedPECNVTypes/cnv_calls.{cnv_type}.bedpe'
    params:
        percent=lambda wildcards, output: output[0].split('0p')[1].split('.')[0]
    output:
        'workchco/BedPEOverlap/cnv_calls.{cnv_type}.0p{percent}.bedpe.overlap.txt'
    shell:
        """
        mkdir -p workchco/BedPEOverlap
        bedtools pairtopair -f 0.{params.percent} -a {input.cnv} -b {input.cnv} -type both | wc -l > {output}.tmp
        # annotate with percent and cnv type
        cat {output}.tmp | awk -v percent={wildcards.percent} '{{print $0"\t0."percent}}' | awk -v cnv_type={wildcards.cnv_type} '{{print $0"\t"cnv_type}}' > {output}
        """


# split the bedpe cnv calls into 100 files with ~620 calls each
rule sample_level_bedpe_overlap_split:
    input:
        cnv='workchco/BedPECNVTypes/cnv_calls.{cnv_type}.bedpe',
    output:
        cnv=expand('workchco/BedPESPLITS/cnv_split_{split}.{cnv_type}.bedpe', split=SPLITS, cnv_type='{cnv_type}')
    shell:
        """
        mkdir -p workchco/BedPESPLITS
        split  --additional-suffix=.{wildcards.cnv_type}.bedpe -d -n l/100 {input} workchco/BedPESPLITS/cnv_split_
        """

rule sample_level_bedpe_breakpoint_overlap:
    input:
        cnv_bed='workchco/BedPESPLITS/cnv_split_{split}.{cnv_type}.bedpe',
        all_cnv='workchco/BedPECNVTypes/cnv_calls.{cnv_type}.bedpe',
    params:
        percent=lambda wildcards, output: output[0].split('0p')[1].split('.')[0]
    output:
        'workchco/BedPEPercentSplits/call_level_reciprocal_overlap.0p{percent}.{split}.{cnv_type}.txt'
    threads: 2
    log: 
        'logs/BedPEPercentSplits/call_level_reciprocal_overlap.0p{percent}.{split}.{cnv_type}.log'
    shell:
        """
        echo 'Percent: {wildcards.percent}'
        mkdir -p workchco/BedPEPercentSplits
        mkdir -p logs/BedPEPercentSplits

        # remove output, if it exists already
        rm -f {output} 2>> {log}
        # for each line in input.cnv_bed
        while IFS= read -r line; do
            echo '-' >> {log}
            echo "$line"  >> {log}
            echo '-' >> {log}
            echo "$line" > workchco/BedPEPercentSplits/call_level_reciprocal_overlap_{wildcards.split}_0p{wildcards.percent}_{wildcards.cnv_type}.tmp.bedpe
            cat workchco/BedPEPercentSplits/call_level_reciprocal_overlap_{wildcards.split}_0p{wildcards.percent}_{wildcards.cnv_type}.tmp.bedpe >> {log}
            bedtools pairtopair -a workchco/BedPEPercentSplits/call_level_reciprocal_overlap_{wildcards.split}_0p{wildcards.percent}_{wildcards.cnv_type}.tmp.bedpe -b {input.all_cnv} -type both -f 0.{params.percent} | wc -l >> {output} 2>> {log}
            echo 'line done' >> {log}
        done < {input.cnv_bed}
        """

rule sample_level_bedpe_slop_overlap:
    input:
        cnv_bed='workchco/BedPESPLITS/cnv_split_{split}.{cnv_type}.bedpe',
        all_cnv='workchco/BedPECNVTypes/cnv_calls.{cnv_type}.bedpe',
    params:
        percent=lambda wildcards, output: output[0].split('0p')[1].split('.')[0]
    output:
        'workchco/BedPESlopPercentSplits/call_level_slop_breakpoint_overlap.0p{percent}.{split}.{cnv_type}.txt'
    threads: 2
    log: 
        'logs/BedPESlopPercentSplits/call_level_slop_breakpoint_overlap.0p{percent}.{split}.{cnv_type}.log'
    shell:
        """
        echo 'Percent: {wildcards.percent}'
        mkdir -p workchco/BedPESlopPercentSplits
        mkdir -p logs/BedPESlopPercentSplits

        # remove output, if it exists already
        rm -f {output} 2>> {log}
        # for each line in input.cnv_bed
        while IFS= read -r line; do
            echo '-' >> {log}
            echo "$line"  >> {log}
            echo '-' >> {log}
            echo "$line" > workchco/BedPESlopPercentSplits/call_level_slop_breakpoint_overlap_{wildcards.split}_0p{wildcards.percent}_{wildcards.cnv_type}.tmp.bedpe
            slop=$(python Scripts/get_slop_from_percent.py workchco/BedPESlopPercentSplits/call_level_slop_breakpoint_overlap_{wildcards.split}_0p{wildcards.percent}_{wildcards.cnv_type}.tmp.bedpe 0.{params.percent})
            cat workchco/BedPESlopPercentSplits/call_level_slop_breakpoint_overlap_{wildcards.split}_0p{wildcards.percent}_{wildcards.cnv_type}.tmp.bedpe >> {log}
            bedtools pairtopair -a workchco/BedPESlopPercentSplits/call_level_slop_breakpoint_overlap_{wildcards.split}_0p{wildcards.percent}_{wildcards.cnv_type}.tmp.bedpe -b {input.all_cnv} -type both -slop $slop | wc -l | awk -v percent={params.percent} '{{print $0"\t0."percent}}' | awk -v cnvtype={wildcards.cnv_type} '{{print $0"\t"cnvtype}}' | awk -v slop=$slop '{{print $0"\t"slop}}' >> {output} 2>> {log}
            echo 'line done' >> {log}
        done < {input.cnv_bed}
        rm workchco/BedPESlopPercentSplits/call_level_slop_breakpoint_overlap_{wildcards.split}_0p{wildcards.percent}_{wildcards.cnv_type}.tmp.bedpe
        """

# aggregate all slop
rule agg_sample_level_bedpe_slop_overlap:
    input:
        expand('workchco/BedPESlopPercentSplits/call_level_slop_breakpoint_overlap.0p{percent}.{split}.{cnv_type}.txt', split=SPLITS, percent=bedpePercents, cnv_type=CNV_TYPES)
    output:
        'workchco/BedPESlopPercentSplits/agg_call_level_slop_breakpoint_overlap.all.txt'
    shell:
        """
        mkdir -p workchco/BedPESlopPercentSplits
        cat {input} > {output}
        """

# plot slop
rule plot_slop_overlap:
    input:
        'workchco/BedPESlopPercentSplits/agg_call_level_slop_breakpoint_overlap.all.txt'
    output:
        'CHCOFigures/slop_breakpoint_bedpe_overlap_2_panel.png'
    shell:
        """
        python Scripts/plot_slop_overlap.py {input} {output} slop breakpoint
        """

# aggregate sample_level_reciprocal_overlap by percentage, cat all files for a % into one file
rule annotate_and_aggregate_sample_level_bedpe_overlap:
    input:
        expand('workchco/BedPEPercentSplits/call_level_reciprocal_overlap.0p{percent}.{split}.{cnv_type}.txt', split=SPLITS, percent='{percent}', cnv_type='{cnv_type}')
    output:
        'workchco/BedPEPercentSplitsAnnotated/agg_call_level_reciprocal_overlap.{percent}.{cnv_type}.txt'
    shell:
        """
        mkdir -p workchco/BedPEPercentSplitsAnnotated
        cat {input} > {output}.tmp
        # annotate with percent and cnv type
        cat {output}.tmp | awk -v percent={wildcards.percent} '{{print $0"\t0."percent}}' | awk -v cnv_type={wildcards.cnv_type} '{{print $0"\t"cnv_type}}' > {output}
        """

rule plot_bedpe_overlap:
    # take all the agg_call_level_reciprocal_overlap files and plot them
    input:
        single_files = expand('workchco/BedPEPercentSplitsAnnotated/agg_call_level_reciprocal_overlap.{percent}.{cnv_type}.txt', percent=bedpePercents, cnv_type=CNV_TYPES),
        all_by_all_files = expand('workchco/BedPEOverlap/cnv_calls.{cnv_type}.0p{percent}.bedpe.overlap.txt', percent=bedpePercents, cnv_type=CNV_TYPES)
    output:
        'CHCOFigures/breakpoint_bedpe_overlap_4_panel.png'
    shell:
        """
        mkdir -p CHCOFigures/
        cat {input.single_files} > workchco/agg_call_level_breakpoint_bedpe_overlap.all.txt
        cat {input.all_by_all_files} > workchco/agg_breakpoint_bedpe_overlap.all.txt
        python Scripts/plot_overlap.py workchco/agg_call_level_breakpoint_bedpe_overlap.all.txt workchco/agg_breakpoint_bedpe_overlap.all.txt {output} breakpoint
        """

# ----------------------------- Ven Diagrams ----------------------------- #

# # divide calls based on caller
rule divide_by_caller:
    input:
        bed='workchco/CNVTypes/cnv_calls.{cnv_type}.bed',
        bedpe='workchco/BedPECNVTypes/cnv_calls.{cnv_type}.bedpe'
    params:
        caller=lambda wildcards, output: output[0].split('.')[2]
    output:
        bed='workchco/CallerSpecificCNVTypes/cnv_calls.{cnv_type}.{caller}.bed',
        bedpe='workchco/CallerSpecificBedPECNVTypes/cnv_calls.{cnv_type}.{caller}.bedpe'
    shell:
        """
        mkdir -p workchco/CallerSpecificCNVTypes
        mkdir -p workchco/CallerSpecificBedPECNVTypes

        grep {params.caller} {input.bed} > {output.bed}
        grep {params.caller} {input.bedpe} > {output.bedpe}
        """

# # do the intersections for bed
rule intersect_callers_bed:
    input:
        bed='workchco/CallerSpecificCNVTypes/cnv_calls.{cnv_type}.{caller}.bed',
        all_beds=expand('workchco/CallerSpecificCNVTypes/cnv_calls.{cnv_type}.{caller}.bed', caller=CALLERS, cnv_type=CNV_TYPES)
    params:
        savvy = 'Savvy',
        cnvkit = 'CNVkit',
        gatk = 'GATK'
    output:
        'workchco/CallerSpecificOverlapsBed/venn.{cnv_type}.{caller}.bed.txt'
    shell:
        """
        mkdir -p workchco/CallerSpecificOverlapsBed
        rm -f {output}

        # ---- savvy ---- #
        STR='workchco/CallerSpecificCNVTypes/cnv_calls.{wildcards.cnv_type}.Savvy.bed'
        SUB='{wildcards.caller}'
        if [[ "$STR" != *"$SUB"* ]]; then
            bedtools intersect -a {input.bed} -b $STR -wo -f 0.60 -r >> {output}
        fi
        
        # ---- cnvkit ---- #
        STR='workchco/CallerSpecificCNVTypes/cnv_calls.{wildcards.cnv_type}.CNVkit.bed'
        SUB='{wildcards.caller}'
        if [[ "$STR" != *"$SUB"* ]]; then
            bedtools intersect -a {input.bed} -b $STR -wo -f 0.60 -r >> {output}
        fi

        # ---- gatk ---- #
        STR='workchco/CallerSpecificCNVTypes/cnv_calls.{wildcards.cnv_type}.GATK.bed'
        SUB='{wildcards.caller}'
        if [[ "$STR" != *"$SUB"* ]]; then
            bedtools intersect -a {input.bed} -b $STR -wo -f 0.60 -r >> {output}
        fi
        """

# find the triple overlaps
rule triple_overlap_bed:
    input:
        expand('workchco/CallerSpecificOverlapsBed/venn.{cnv_type}.{caller}.bed.txt', caller=CALLERS, cnv_type='{cnv_type}')
    output:
        'workchco/CallerSpecificOverTripleOverlap/venn.{cnv_type}.triple.bed.txt'
    shell:
        """
        mkdir -p workchco/CallerSpecificOverTripleOverlap
        python Scripts/triple_overlap.py -s workchco/CallerSpecificOverlapsBed/venn.{wildcards.cnv_type}.Savvy.bed.txt \
                -g workchco/CallerSpecificOverlapsBed/venn.{wildcards.cnv_type}.GATK.bed.txt \
                -c workchco/CallerSpecificOverlapsBed/venn.{wildcards.cnv_type}.CNVkit.bed.txt > {output}
        """

# report the number of single, double and triple overlaps
rule report_overlap_numbers:
    input:
        triple='workchco/CallerSpecificOverTripleOverlap/venn.{cnv_type}.triple.bed.txt',
        double=expand('workchco/CallerSpecificOverlapsBed/venn.{cnv_type}.{caller}.bed.txt', caller=CALLERS, cnv_type='{cnv_type}'),
        single=expand('workchco/CallerSpecificCNVTypes/cnv_calls.{cnv_type}.{caller}.bed', caller=CALLERS, cnv_type='{cnv_type}')
    output:
        non_specific='workchco/VenDiagramResults/overlap_numbers.{cnv_type}.txt',
        sample_specific='workchco/VenDiagramResults/overlap_numbers.{cnv_type}.sample_specific.txt',
        sample_specific_triples='workchco/VenDiagramResults/overlap_numbers.{cnv_type}.sample_specific.triples.bed.txt',
        sample_specific_triple_calls='workchco/VenDiagramResults/overlapping_calls.{cnv_type}.sample_specific.triples.bed.txt',
        sample_specific_all_groupings=expand('workchco/RecipCategoriesBeds/{cnv_type}{end}.bed',cnv_type='{cnv_type}',end=['-savvy','-cnvkit','-gatk','-triple','-savvy_gatk','-savvy_cnvkit','-cnvkit_gatk'])
    shell:
        """
        mkdir -p workchco/VenDiagramResults
        mkdir -p workchco/RecipCategoriesBeds
        # number of calls in each caller
        python Scripts/count_overlaps.py -s workchco/CallerSpecificCNVTypes/cnv_calls.{wildcards.cnv_type}.Savvy.bed \
            -c workchco/CallerSpecificCNVTypes/cnv_calls.{wildcards.cnv_type}.CNVkit.bed \
            -g workchco/CallerSpecificCNVTypes/cnv_calls.{wildcards.cnv_type}.GATK.bed \
            --sx workchco/CallerSpecificOverlapsBed/venn.{wildcards.cnv_type}.Savvy.bed.txt \
            --cx workchco/CallerSpecificOverlapsBed/venn.{wildcards.cnv_type}.CNVkit.bed.txt \
            --gx workchco/CallerSpecificOverlapsBed/venn.{wildcards.cnv_type}.GATK.bed.txt \
            -t {input.triple} \
            -o {output.non_specific}

        # force calls to be from the same sample
         python Scripts/count_overlaps.py -s workchco/CallerSpecificCNVTypes/cnv_calls.{wildcards.cnv_type}.Savvy.bed \
            -c workchco/CallerSpecificCNVTypes/cnv_calls.{wildcards.cnv_type}.CNVkit.bed \
            -g workchco/CallerSpecificCNVTypes/cnv_calls.{wildcards.cnv_type}.GATK.bed \
            --sx workchco/CallerSpecificOverlapsBed/venn.{wildcards.cnv_type}.Savvy.bed.txt \
            --cx workchco/CallerSpecificOverlapsBed/venn.{wildcards.cnv_type}.CNVkit.bed.txt \
            --gx workchco/CallerSpecificOverlapsBed/venn.{wildcards.cnv_type}.GATK.bed.txt \
            -t {input.triple} \
            -o {output.sample_specific} --fs --to {output.sample_specific_triples} \
            --sample_specific_triples_out {output.sample_specific_triple_calls} \
            --prefix workchco/RecipCategoriesBeds/{wildcards.cnv_type}

        """

# ------------------------ BEDPE Ven Diagrams ------------------------ #
# do the intersections for bed
rule intersect_callers_bedpe:
    input:
        bedpe='workchco/CallerSpecificBedPECNVTypes/cnv_calls.{cnv_type}.{caller}.bedpe',
        all_bedpes=expand('workchco/CallerSpecificBedPECNVTypes/cnv_calls.{cnv_type}.{caller}.bedpe', caller=CALLERS, cnv_type=CNV_TYPES)
    params:
        savvy = 'Savvy',
        cnvkit = 'CNVkit',
        gatk = 'GATK'
    output:
        'workchco/CallerSpecificOverlapsBedPE/venn.{cnv_type}.{caller}.bedpe.txt'
    shell:
        """
        mkdir -p workchco/CallerSpecificOverlapsBedPE
        rm -f {output}

        # ---- savvy ---- #
        STR='workchco/CallerSpecificBedPECNVTypes/cnv_calls.{wildcards.cnv_type}.Savvy.bedpe'
        SUB='{wildcards.caller}'
        if [[ "$STR" != *"$SUB"* ]]; then
            while IFS= read -r line; do
                echo "$line" > workchco/CallerSpecificOverlapsBedPE/tmp.{wildcards.cnv_type}.{wildcards.caller}.bedpe
                slop=$(python Scripts/get_slop_from_percent.py workchco/CallerSpecificOverlapsBedPE/tmp.{wildcards.cnv_type}.{wildcards.caller}.bedpe 0.50)
                bedtools pairtopair -a workchco/CallerSpecificOverlapsBedPE/tmp.{wildcards.cnv_type}.{wildcards.caller}.bedpe -b $STR -slop $slop >> {output}
            done < {input.bedpe}
        fi
        
        # ---- cnvkit ---- #
        STR='workchco/CallerSpecificBedPECNVTypes/cnv_calls.{wildcards.cnv_type}.CNVkit.bedpe'
        SUB='{wildcards.caller}'
        if [[ "$STR" != *"$SUB"* ]]; then
            while IFS= read -r line; do
                echo "$line" > workchco/CallerSpecificOverlapsBedPE/tmp.{wildcards.cnv_type}.{wildcards.caller}.bedpe
                slop=$(python Scripts/get_slop_from_percent.py workchco/CallerSpecificOverlapsBedPE/tmp.{wildcards.cnv_type}.{wildcards.caller}.bedpe 0.50)
                bedtools pairtopair -a workchco/CallerSpecificOverlapsBedPE/tmp.{wildcards.cnv_type}.{wildcards.caller}.bedpe -b $STR -slop $slop >> {output}
            done < {input.bedpe}
        fi

        # ---- gatk ---- #
        STR='workchco/CallerSpecificBedPECNVTypes/cnv_calls.{wildcards.cnv_type}.GATK.bedpe'
        SUB='{wildcards.caller}'
        if [[ "$STR" != *"$SUB"* ]]; then
            while IFS= read -r line; do
                echo "$line" > workchco/CallerSpecificOverlapsBedPE/tmp.{wildcards.cnv_type}.{wildcards.caller}.bedpe
                slop=$(python Scripts/get_slop_from_percent.py workchco/CallerSpecificOverlapsBedPE/tmp.{wildcards.cnv_type}.{wildcards.caller}.bedpe 0.50)
                bedtools pairtopair -a workchco/CallerSpecificOverlapsBedPE/tmp.{wildcards.cnv_type}.{wildcards.caller}.bedpe -b $STR -slop $slop >> {output}
            done < {input.bedpe}
        fi
        """

rule bedpe_triple_overlap:
    input:
        expand('workchco/CallerSpecificOverlapsBedPE/venn.{cnv_type}.{caller}.bedpe.txt', caller=CALLERS, cnv_type='{cnv_type}')
    output:
        'workchco/CallerSpecificTripleOverlapBedPE/venn.{cnv_type}.triple.bedpe.txt'
    shell:
        """
        mkdir -p workchco/CallerSpecificTripleOverlapBedPE

        python Scripts/triple_overlap.py -s workchco/CallerSpecificOverlapsBedPE/venn.{wildcards.cnv_type}.Savvy.bedpe.txt \
                    -g workchco/CallerSpecificOverlapsBedPE/venn.{wildcards.cnv_type}.GATK.bedpe.txt \
                    -c workchco/CallerSpecificOverlapsBedPE/venn.{wildcards.cnv_type}.CNVkit.bedpe.txt -p > {output}
        """

# report the number of single, double and triple overlaps
rule report_overlap_numbers_bedpe:
    input:
        triple='workchco/CallerSpecificTripleOverlapBedPE/venn.{cnv_type}.triple.bedpe.txt',
        double=expand('workchco/CallerSpecificOverlapsBedPE/venn.{cnv_type}.{caller}.bedpe.txt', caller=CALLERS, cnv_type='{cnv_type}'),
        single=expand('workchco/CallerSpecificBedPECNVTypes/cnv_calls.{cnv_type}.{caller}.bedpe', caller=CALLERS, cnv_type='{cnv_type}')
    output:
        non_specific='workchco/VenDiagramResultsBedPE/overlap_numbers.{cnv_type}.txt',
        sample_specific='workchco/VenDiagramResultsBedPE/overlap_numbers.{cnv_type}.sample_specific.txt',
        sample_specific_triples='workchco/VenDiagramResultsBedPE/overlap_numbers.{cnv_type}.sample_specific.triples.bed.txt',
        sample_specific_all_groupings=expand('workchco/RecipCategoriesBedPEs/{cnv_type}{end}.bed',cnv_type='{cnv_type}',end=['-savvy','-cnvkit','-gatk','-triple','-savvy_gatk','-savvy_cnvkit','-cnvkit_gatk'])
    shell:
        """
        mkdir -p workchco/VenDiagramResultsBedPE
        mkdir -p workchco/RecipCategoriesBedPEs
        # number of calls in each caller
        python Scripts/count_overlaps.py -s workchco/CallerSpecificBedPECNVTypes/cnv_calls.{wildcards.cnv_type}.Savvy.bedpe \
            -c workchco/CallerSpecificBedPECNVTypes/cnv_calls.{wildcards.cnv_type}.CNVkit.bedpe \
            -g workchco/CallerSpecificBedPECNVTypes/cnv_calls.{wildcards.cnv_type}.GATK.bedpe \
            --sx workchco/CallerSpecificOverlapsBedPE/venn.{wildcards.cnv_type}.Savvy.bedpe.txt \
            --cx workchco/CallerSpecificOverlapsBedPE/venn.{wildcards.cnv_type}.CNVkit.bedpe.txt \
            --gx workchco/CallerSpecificOverlapsBedPE/venn.{wildcards.cnv_type}.GATK.bedpe.txt \
            -t {input.triple} \
            -o {output.non_specific} -p
        
        python Scripts/count_overlaps.py -s workchco/CallerSpecificBedPECNVTypes/cnv_calls.{wildcards.cnv_type}.Savvy.bedpe \
            -c workchco/CallerSpecificBedPECNVTypes/cnv_calls.{wildcards.cnv_type}.CNVkit.bedpe \
            -g workchco/CallerSpecificBedPECNVTypes/cnv_calls.{wildcards.cnv_type}.GATK.bedpe \
            --sx workchco/CallerSpecificOverlapsBedPE/venn.{wildcards.cnv_type}.Savvy.bedpe.txt \
            --cx workchco/CallerSpecificOverlapsBedPE/venn.{wildcards.cnv_type}.CNVkit.bedpe.txt \
            --gx workchco/CallerSpecificOverlapsBedPE/venn.{wildcards.cnv_type}.GATK.bedpe.txt \
            -t {input.triple} \
            -o {output.sample_specific} --fs -p --to {output.sample_specific_triples} \
            --prefix workchco/RecipCategoriesBedPEs/{wildcards.cnv_type}
        """


# ------------------------------ BED + BEDPE Overlaps Venn Diagrams ------------------------------ #

rule report_all_bed_and_bedpe_overlap:
    input:
        bedpe_triple='workchco/CallerSpecificTripleOverlapBedPE/venn.{cnv_type}.triple.bedpe.txt',
        bedpe_double=expand('workchco/CallerSpecificOverlapsBedPE/venn.{cnv_type}.{caller}.bedpe.txt', caller=CALLERS, cnv_type='{cnv_type}'),
        triple='workchco/CallerSpecificOverTripleOverlap/venn.{cnv_type}.triple.bed.txt',
        double=expand('workchco/CallerSpecificOverlapsBed/venn.{cnv_type}.{caller}.bed.txt', caller=CALLERS, cnv_type='{cnv_type}'),
    output:
        non_specific='workchco/VenDiagramResultsBoth/venn.bed_and_bedpe.{cnv_type}.txt',
        sample_specific='workchco/VenDiagramResultsBoth/venn.bed_and_bedpe.{cnv_type}.sample_specific.txt',
        sample_specific_triples='workchco/VenDiagramResultsBoth/venn.bed_and_bedpe.{cnv_type}.sample_specific.triples.bed.txt',
    shell:
        """
        mkdir -p workchco/VenDiagramResultsBoth
        python Scripts/get_bed_bedpe_overlaps.py \
            --sb workchco/CallerSpecificOverlapsBed/venn.{wildcards.cnv_type}.Savvy.bed.txt \
            --sp workchco/CallerSpecificOverlapsBedPE/venn.{wildcards.cnv_type}.Savvy.bedpe.txt \
            --cb workchco/CallerSpecificOverlapsBed/venn.{wildcards.cnv_type}.CNVkit.bed.txt \
            --cp  workchco/CallerSpecificOverlapsBedPE/venn.{wildcards.cnv_type}.CNVkit.bedpe.txt \
            --gb workchco/CallerSpecificOverlapsBed/venn.{wildcards.cnv_type}.GATK.bed.txt \
            --gp  workchco/CallerSpecificOverlapsBedPE/venn.{wildcards.cnv_type}.GATK.bedpe.txt \
            --tb {input.triple} \
            --tp {input.bedpe_triple} \
            --sx workchco/VenDiagramResultsBoth/tmp.{wildcards.cnv_type}.Savvy.bed.txt \
            --cx workchco/VenDiagramResultsBoth/tmp.{wildcards.cnv_type}.CNVkit.bed.txt \
            --gx workchco/VenDiagramResultsBoth/tmp.{wildcards.cnv_type}.GATK.bed.txt \
            --tx workchco/VenDiagramResultsBoth/tmp.{wildcards.cnv_type}.triple.bed.txt
        
        # overlap again but force the samples to be the same
        python Scripts/count_overlaps.py -s workchco/CallerSpecificCNVTypes/cnv_calls.{wildcards.cnv_type}.Savvy.bed \
            -c workchco/CallerSpecificCNVTypes/cnv_calls.{wildcards.cnv_type}.CNVkit.bed \
            -g workchco/CallerSpecificCNVTypes/cnv_calls.{wildcards.cnv_type}.GATK.bed \
            --sx workchco/VenDiagramResultsBoth/tmp.{wildcards.cnv_type}.Savvy.bed.txt \
            --cx workchco/VenDiagramResultsBoth/tmp.{wildcards.cnv_type}.CNVkit.bed.txt \
            --gx workchco/VenDiagramResultsBoth/tmp.{wildcards.cnv_type}.GATK.bed.txt \
            -t  workchco/VenDiagramResultsBoth/tmp.{wildcards.cnv_type}.triple.bed.txt \
            -o {output.non_specific}
        
        # overlap again but force the samples to be the same
        python Scripts/count_overlaps.py -s workchco/CallerSpecificCNVTypes/cnv_calls.{wildcards.cnv_type}.Savvy.bed \
            -c workchco/CallerSpecificCNVTypes/cnv_calls.{wildcards.cnv_type}.CNVkit.bed \
            -g workchco/CallerSpecificCNVTypes/cnv_calls.{wildcards.cnv_type}.GATK.bed \
            --sx workchco/VenDiagramResultsBoth/tmp.{wildcards.cnv_type}.Savvy.bed.txt \
            --cx workchco/VenDiagramResultsBoth/tmp.{wildcards.cnv_type}.CNVkit.bed.txt \
            --gx workchco/VenDiagramResultsBoth/tmp.{wildcards.cnv_type}.GATK.bed.txt \
            -t  workchco/VenDiagramResultsBoth/tmp.{wildcards.cnv_type}.triple.bed.txt \
            -o {output.sample_specific} --fs --to {output.sample_specific_triples}

        """

rule merge_all_triple_calls:
    input:
        expand('workchco/VenDiagramResultsBoth/venn.bed_and_bedpe.{cnv_type}.sample_specific.triples.bed.txt', cnv_type=CNV_TYPES),
        expand('workchco/VenDiagramResults/overlap_numbers.{cnv_type}.sample_specific.triples.bed.txt', cnv_type=CNV_TYPES),
        expand('workchco/VenDiagramResultsBedPE/overlap_numbers.{cnv_type}.sample_specific.triples.bed.txt', cnv_type=CNV_TYPES)
    output:
        'workchco/AllSampleSpecificTripleCalls/venn.bed_and_bedpe.all.sample_specific.triples.bed.txt'
    shell:
        """
        cat {input}  | sort | uniq > {output}
        """

# ------------------------------ Plot Reciprocal & Breakpoint ------------------------------ #

# plot slop
rule plot_slop_and_recipt:
    input:
        slop='workchco/BedPESlopPercentSplits/agg_call_level_slop_breakpoint_overlap.all.txt',
        single_files_r = expand('workchco/AnnotatedPercentSplits/agg_call_level_reciprocal_overlap.{percent}.{cnv_type}.annotated.txt', percent=PERCENTS, cnv_type=CNV_TYPES),
        all_by_all_files_r = expand('workchco/AnnotatedRecip/reciprocal_overlap.{percent}.{cnv_type}.annotated.txt', percent=PERCENTS, cnv_type=CNV_TYPES)
    output:
        'CHCOFigures/reciprocal_and_slop_breakpoint_boxplots.png'
    shell:
        """
        mkdir -p CHCOFigures/
        # aggregate files for reciprocal overlap
        cat {input.single_files_r} > workchco/agg_call_level_reciprocal_overlap.all.txt
        cat {input.all_by_all_files_r} > workchco/agg_reciprocal_overlap.all.txt
        # plot it!
        Scripts/plot_recip_and_slop.py {input.slop} workchco/agg_call_level_reciprocal_overlap.all.txt {output}
        """

rule plot_sizes:
    input:
        del_gatk='workchco/CallerSpecificCNVTypes/cnv_calls.DEL.GATK.bed',
        dup_gatk='workchco/CallerSpecificCNVTypes/cnv_calls.DUP.GATK.bed',
        del_cnvkit='workchco/CallerSpecificCNVTypes/cnv_calls.DEL.CNVkit.bed',
        dup_cnvkit='workchco/CallerSpecificCNVTypes/cnv_calls.DUP.CNVkit.bed',
        del_savvy='workchco/CallerSpecificCNVTypes/cnv_calls.DEL.Savvy.bed',
        dup_savvy='workchco/CallerSpecificCNVTypes/cnv_calls.DUP.Savvy.bed'
    output:
        size='CHCOFigures/size_distribution.png',
        calls='CHCOFigures/calls_per_sample.png'
    shell:
        """
        mkdir -p CHCOFigures/
        # cat inputs into single caller specific input tmp files
        cat {input.del_gatk} > workchco/CallerSpecificCNVTypes/tmp.GATK.bed
        cat {input.dup_gatk} >> workchco/CallerSpecificCNVTypes/tmp.GATK.bed
        cat {input.del_cnvkit} > workchco/CallerSpecificCNVTypes/tmp.CNVkit.bed
        cat {input.dup_cnvkit} >> workchco/CallerSpecificCNVTypes/tmp.CNVkit.bed
        cat {input.del_savvy} > workchco/CallerSpecificCNVTypes/tmp.Savvy.bed
        cat {input.dup_savvy} >> workchco/CallerSpecificCNVTypes/tmp.Savvy.bed
        # plot it!
        python Scripts/plot_sizes.py -i workchco/CallerSpecificCNVTypes/tmp.GATK.bed \
            -i workchco/CallerSpecificCNVTypes/tmp.Savvy.bed \
            -i workchco/CallerSpecificCNVTypes/tmp.CNVkit.bed \
            -l gCNV \
            -l Savvy \
            -l CNVkit \
            -o {output.size}


        # cat all dels into one file
        cat {input.del_gatk} > workchco/CallerSpecificCNVTypes/tmp.DEL.bed
        cat {input.del_cnvkit} >> workchco/CallerSpecificCNVTypes/tmp.DEL.bed
        cat {input.del_savvy} >> workchco/CallerSpecificCNVTypes/tmp.DEL.bed
        # cat all dups into one file
        cat {input.dup_gatk} > workchco/CallerSpecificCNVTypes/tmp.DUP.bed
        cat {input.dup_cnvkit} >> workchco/CallerSpecificCNVTypes/tmp.DUP.bed
        cat {input.dup_savvy} >> workchco/CallerSpecificCNVTypes/tmp.DUP.bed

        python Scripts/plot_calls_per_sample.py --dels workchco/CallerSpecificCNVTypes/tmp.DEL.bed --dups workchco/CallerSpecificCNVTypes/tmp.DUP.bed --output {output.calls}
        """

rule upsetplot:
    input:
        recip_del = 'workchco/VenDiagramResults/overlap_numbers.DEL.sample_specific.txt',
        recip_dup = 'workchco/VenDiagramResults/overlap_numbers.DUP.sample_specific.txt',
        slop_del = 'workchco/VenDiagramResultsBedPE/overlap_numbers.DEL.sample_specific.txt',
        slop_dup = 'workchco/VenDiagramResultsBedPE/overlap_numbers.DUP.sample_specific.txt'
    output:
        recip_del='CHCOFigures/upsetplot.DEL.png',
        recip_dup='CHCOFigures/upsetplot.DUP.png',
        slop_del='CHCOFigures/upsetplot.slop.DEL.png',
        slop_dup='CHCOFigures/upsetplot.slop.DUP.png'
    shell:
        """
        python Scripts/upset_plot.py -i {input.recip_del} -l Deletion -o {output.recip_del}
        python Scripts/upset_plot.py -i {input.recip_dup} -l Duplication -o {output.recip_dup}

        python Scripts/upset_plot.py -i {input.slop_del} -l Deletion -o {output.slop_del}
        python Scripts/upset_plot.py -i {input.slop_dup} -l Duplication -o {output.slop_dup}
        """
