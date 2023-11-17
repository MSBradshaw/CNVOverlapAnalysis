# # generate stats and figures about each threshold experiment
# snakemake -s 1kg_sv_validation.smk --cores 64 --config percentage_param=0.1
# snakemake -s 1kg_sv_validation.smk --cores 64 --config percentage_param=0.2
# snakemake -s 1kg_sv_validation.smk --cores 64 --config percentage_param=0.3
# snakemake -s 1kg_sv_validation.smk --cores 64 --config percentage_param=0.4
# snakemake -s 1kg_sv_validation.smk --cores 64 --config percentage_param=0.5
# snakemake -s 1kg_sv_validation.smk --cores 64 --config percentage_param=0.6
# snakemake -s 1kg_sv_validation.smk --cores 64 --config percentage_param=0.7
# snakemake -s 1kg_sv_validation.smk --cores 64 --config percentage_param=0.8
# snakemake -s 1kg_sv_validation.smk --cores 64 --config percentage_param=0.9

mkdir -p 1KG_Validation_experiments

cp 1kg_sv_validation0.1/validated_percent_overlap_stats.txt 1KG_Validation_experiments/1kg_sv_validation0.1.txt
cp 1kg_sv_validation0.2/validated_percent_overlap_stats.txt 1KG_Validation_experiments/1kg_sv_validation0.2.txt
cp 1kg_sv_validation0.3/validated_percent_overlap_stats.txt 1KG_Validation_experiments/1kg_sv_validation0.3.txt
cp 1kg_sv_validation0.4/validated_percent_overlap_stats.txt 1KG_Validation_experiments/1kg_sv_validation0.4.txt
cp 1kg_sv_validation0.5/validated_percent_overlap_stats.txt 1KG_Validation_experiments/1kg_sv_validation0.5.txt
cp 1kg_sv_validation0.6/validated_percent_overlap_stats.txt 1KG_Validation_experiments/1kg_sv_validation0.6.txt
cp 1kg_sv_validation0.7/validated_percent_overlap_stats.txt 1KG_Validation_experiments/1kg_sv_validation0.7.txt
cp 1kg_sv_validation0.8/validated_percent_overlap_stats.txt 1KG_Validation_experiments/1kg_sv_validation0.8.txt
cp 1kg_sv_validation0.9/validated_percent_overlap_stats.txt 1KG_Validation_experiments/1kg_sv_validation0.9.txt

python Scripts/plot_1kg_sv_threshold_experiment.py
cp 1kg_sv_validation0.3/validated_percent_upset_plot.DEL.png Figures/1kg_sv_validation0.3_upset_plot.DEL.png
cp 1kg_sv_validation0.5/validated_percent_upset_plot.DUP.png Figures/1kg_sv_validation0.5_upset_plot.DUP.png
cp 1kg_sv_validation0.3/upset_plot.DEL.png Figures/1kg_upset_plot.DEL.png
cp 1kg_sv_validation0.5/upset_plot.DUP.png Figures/1kg_upset_plot.DUP.png