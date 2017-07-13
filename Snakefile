"""
Reproducible analysis for Ouija paper

To build the pandoc files, call:
export PATH=$PATH:/Applications/RStudio.app/Contents/MacOS/pandoc

This uses a ridiculous number of DLLs, so you might need to set
export R_MAX_NUM_DLLS=200
"""


# Configuration --------------

# Marker-vs-transcriptome variables
index = [1,2,3] # Three datasets
rep = list(range(1,51)) # 50 repetitions
num_additional_markers = [1,5,10,20,50,100,500,1000] # Number of additional markers
nam_index = range(1,len(num_additional_markers)+1)

output_csv = expand("data/marker-vs-txome/ouija_fits/ouija_pseudotime_{i}_{r}_{nam}.csv",
			i = index, r = rep, nam = nam_index)

output_de = expand("data/marker-vs-txome/subset_de_fits/sde_fit_{i}_{r}_{nam}.csv",
			i = index, r = rep, nam = nam_index)


# Benchmarking variables 
# Gs = [6,9,12,15]
# condition = ["true", "noninformative"]
# reps = list(range(1,41))
# regimes = ["logit", "probit", "cloglog", "threshold"]

# ouija_csv = expand('data/benchmarking/{regime}/ouija/ouija_{cond}_{G}_{rep}.csv', cond = condition, G = Gs, rep = reps, regime = regimes)
# pca_files = expand('data/benchmarking/{regime}/pca.rds', regime = regimes)
# dpt_files = expand('data/benchmarking/{regime}/dpt.rds', regime = regimes)


R_opts = "--vanilla"


rule all:
    input:
        "figs/marker-vs-txome-de-results.png",
        "figs/marker-vs-txome-pseudotime-subsets.png"


# Construct SCESets ---------------------

rule read_hsc:
    output:
        "data/scesets/hsc-sce.rds"
    shell:
        "Rscript {R_opts} scripts/read-hsc.R"

rule create_chu:
    output:
        "data/scesets/chu-sce.rds"
    shell:
        "Rscript {R_opts} scripts/read-chu.R"

rule create_dulken:
    output:
        "data/scesets/dulken-sce.rds"
    shell:
        "Rscript {R_opts} scripts/read-dulken.R"

rule create_li:
    output:
        "data/scesets/li-sce.rds"
    shell:
        "Rscript {R_opts} scripts/read-li.R"

rule construct_scesets:
    input:
        "data/scesets/hsc-sce.rds",
        "data/raw/waterfall_data.xlsx",
        "data/raw/hsmm.rds"
    output:
        "data/scesets/scesets.rds"
    shell:
        "Rscript {R_opts} scripts/create_scesets.R"

# Analyses ------

rule chu_analysis:
    input:
        "data/scesets/chu-sce.rds"
    output:
        "analysis/datasets/chu.html",
        "data/mvt_csv/chu.csv",
        "figs/fig2.png"
    shell:
        "Rscript -e \"rmarkdown::render('analysis/datasets/chu.Rmd')\""

rule dulken_analysis:
    input:
        "data/scesets/dulken-sce.rds"
    output:
        "analysis/datasets/dulken.html",
        "data/mvt_csv/dulken.csv"
    shell:
        "Rscript -e \"rmarkdown::render('analysis/datasets/dulken.Rmd')\""

rule shin_analysis:
    output:
        "analysis/datasets/shin.html",
        "data/mvt_csv/shin.csv"
    shell:
        "Rscript -e \"rmarkdown::render('analysis/datasets/shin.Rmd')\""

rule trapnell_analysis:
    output:
        "analysis/datasets/trapnell.html",
        "data/mvt_csv/trapnell.csv",
        "figs/trapnell_correlation.rds",
        "figs/trapnell_example_genes.rds"
    shell:
        "Rscript -e \"rmarkdown::render('analysis/datasets/trapnell.Rmd')\""

rule zhou_analysis:
    output:
        "analysis/datasets/zhou.html",
        "data/mvt_csv/zhou.csv",
        "figs/fig3.png"
    shell:
        "Rscript -e \"rmarkdown::render('analysis/datasets/zhou.Rmd')\""

rule li_analysis:
    output:
        "analysis/datasets/li.html",
        "data/mvt_csv/li.csv"
    shell:
        "Rscript -e \"rmarkdown::render('analysis/datasets/li.Rmd')\""


rule marker_correlation_figs:
    input:
        "data/mvt_csv/trapnell.csv",
        "data/mvt_csv/shin.csv",
        "data/mvt_csv/zhou.csv",
        "data/mvt_csv/chu.csv",
        "data/mvt_csv/dulken.csv",
    output:
        "figs/marker_correlations.rds",
        "figs/relative_marker_correlations.rds"
    shell:
        "Rscript {R_opts} analysis/datasets/create-comparison-figure.R"

# Figures --------

rule figure_1:
    input:
        "figs/marker_correlations.rds",
        "figs/trapnell_correlation.rds",
        "figs/trapnell_example_genes.rds"
    output:
        "figs/fig1.png"
    shell:
        "Rscript {R_opts} analysis/datasets/create-figure-1.R"




# Benchmarking switch-like models --------------------

# rule create_synthetic_probit:
#     output:
#         "data/benchmarking/probit_synthetic.h5",
#         "figs/example_genes_probit.png"
#     shell:
#         "Rscript analysis/benchmarking/create-synthetic-probit.R"

# rule create_synthetic_logit:
#     output:
#         "data/benchmarking/logit_synthetic.h5",
#         "figs/example_genes_logit.png"
#     shell:
#         "Rscript analysis/benchmarking/create-synthetic-logit.R"

# rule create_synthetic_threshold:
#     output:
#         "data/benchmarking/threshold_synthetic.h5",
#         "figs/example_genes_threshold.png"
#     shell:
#         "Rscript analysis/benchmarking/create-synthetic-threshold.R"

# rule create_synthetic_cloglog:
#     output:
#         "data/benchmarking/cloglog_synthetic.h5",
#         "figs/example_genes_cloglog.png"
#     shell:
#         "Rscript analysis/benchmarking/create-synthetic-cloglog.R"

# rule benchmark_pca:
# 	input:
# 		"data/benchmarking/logit_synthetic.h5",
#         "data/benchmarking/probit_synthetic.h5",
#         "data/benchmarking/cloglog_synthetic.h5",
#         "data/benchmarking/threshold_synthetic.h5"
# 	output:
# 		"data/benchmarking/{regime}/pca.rds"
# 	shell:
# 		'Rscript {R_opts} analysis/benchmarking/pca.R {input} {output}'

# rule benchmark_dpt:
# 	input:
# 		"data/benchmarking/logit_synthetic.h5",
#         "data/benchmarking/probit_synthetic.h5",
#         "data/benchmarking/cloglog_synthetic.h5",
#         "data/benchmarking/threshold_synthetic.h5"
# 	output:
# 		"data/benchmarking/{regime}/dpt.rds"
# 	shell:
# 		'Rscript {R_opts} analysis/benchmarking/dpt.R {input} {output}'

# rule ouija_probit:
#     input:
#         synthetic_data_file
#     output:
#         'data/benchmarking/probit/ouija_{cond}_{G}_{rep}.csv'
#     shell:
#         "Rscript {R_opts} data/benchmarking/ouija_probit.R --condition {wildcards.cond} --rep {wildcards.rep} --G {wildcards.G} --input_file {input} --output_file {output}"

# rule ouija_cloglog:
#     input:
#         synthetic_data_file
#     output:
#         "data/benchmarking/cloglog/ouija_{cond}_{G}_{rep}.csv"
#     shell:
#         "Rscript {R_opts} data/benchmarking/ouija_cloglog.R --condition {wildcards.cond} --rep {wildcards.rep} --G {wildcards.G} --input_file {input} --output_file {output}"

# rule ouija_threshold:
#     input:
#         synthetic_data_file
#     output:
#         "data/benchmarking/threshold/ouija_{cond}_{G}_{rep}.csv"
#     shell:
#         "Rscript {R_opts} data/benchmarking/ouija_threshold.R --condition {wildcards.cond} --rep {wildcards.rep} --G {wildcards.G}--input_file {input} --output_file {output}"

