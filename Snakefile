"""
Reproducible analysis for Ouija paper

To build the pandoc files, call:
export PATH=$PATH:/Applications/RStudio.app/Contents/MacOS/pandoc

This uses a ridiculous number of DLLs, so you might need to set
export R_MAX_NUM_DLLS=200
"""


# Configuration --------------


# Benchmarking variables 
Gs = [6,9,12,15]
condition = ["true", "noninformative"]

reps = list(range(1,501))
regimes = ["logit", "probit", "cloglog", "threshold"]

ouija_csv = expand('data/benchmarking/{regime}/ouija_{cond}_{G}_{rep}.csv', cond = condition, G = Gs, rep = reps, regime = regimes[1:4])


synthetic_files = expand("data/benchmarking/{regime}_synthetic.h5", regime = regimes)
benchmark_results = expand("data/benchmarking/{regime}_{algorithm}.rds", regime = regimes, algorithm = ["dpt", "monocle", "tscan"])


# Transient variables
Gts = [8, 12, 16, 24]
prop_switch = [0.75, 0.5]
N_rep_transient = list(range(1, 101))

transient_synthetic_files = expand("data/transient/synthetic/transient_sim_{Gt}_{ps}_{rept}.csv",
                                    Gt = Gts, ps = prop_switch, rept = N_rep_transient)
transient_ouija_files = expand("data/transient/results/ouija_cor_{Gt}_{ps}_{rept}.csv",
                                    Gt = Gts, ps = prop_switch, rept = N_rep_transient)



R_opts = "--vanilla"


rule all:
    input:
        #"figs/fig_interpretable.png",
        # "figs/fig_benchmark.png",
        # "figs/fig1.png",
        # "figs/fig2.png",
        # "figs/fig3.png",
        # "figs/fig4.png",
        # synthetic_files,
        transient_synthetic_files, transient_ouija_files


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
        "figs/fig_interpretable.png"
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
        "data/mvt_csv/shin.csv",
        "figs/supp_shin.png"
    shell:
        "Rscript -e \"rmarkdown::render('analysis/datasets/shin.Rmd')\""

rule trapnell_analysis:
    output:
        "analysis/datasets/trapnell.html",
        "data/mvt_csv/trapnell.csv",
        "figs/trapnell_correlation.rds",
        "figs/trapnell_example_genes.rds",
        "figs/supp_trapnell.png"
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
    input:
        "data/scesets/li-sce.rds"
    output:
        "analysis/datasets/li.html",
        "data/mvt_csv/li.csv",
        "figs/li_cor.rds",
        "figs/li_correct_genes.rds",
        "figs/li_incorrect_genes.rds"
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

rule create_synthetic_probit:
    output:
        "data/benchmarking/probit_synthetic.h5",
        "figs/example_genes_probit.png"
    shell:
        "Rscript analysis/benchmarking/create-synthetic-probit.R"

rule create_synthetic_logit:
    output:
        "data/benchmarking/logit_synthetic.h5",
        "figs/example_genes_logit.png"
    shell:
        "Rscript analysis/benchmarking/create-synthetic-logit.R"

rule create_synthetic_threshold:
    output:
        "data/benchmarking/threshold_synthetic.h5",
        "figs/example_genes_threshold.png"
    shell:
        "Rscript analysis/benchmarking/create-synthetic-threshold.R"

rule create_synthetic_cloglog:
    output:
        "data/benchmarking/cloglog_synthetic.h5",
        "figs/example_genes_cloglog.png"
    shell:
        "Rscript analysis/benchmarking/create-synthetic-cloglog.R"


rule benchmark_dpt:
    input:
        "data/benchmarking/{regime}_synthetic.h5"
    output:
        "data/benchmarking/{regime}_dpt.rds"
    shell:
        "Rscript {R_opts} analysis/benchmarking/fit_pseudotime.R dpt {input} {output}"

rule benchmark_monocle:
    input:
        "data/benchmarking/{regime}_synthetic.h5"
    output:
        "data/benchmarking/{regime}_monocle.rds"
    shell:
        "Rscript {R_opts} analysis/benchmarking/fit_pseudotime.R monocle {input} {output}"

rule benchmark_tscan:
    input:
        "data/benchmarking/{regime}_synthetic.h5"
    output:
        "data/benchmarking/{regime}_tscan.rds"
    shell:
        'Rscript {R_opts} analysis/benchmarking/fit_pseudotime.R tscan {input} {output}'

rule ouija_probit:
    input:
        "data/benchmarking/probit_synthetic.h5"
    output:
        'data/benchmarking/probit/ouija_{cond}_{G}_{rep}.csv'
    shell:
        "Rscript {R_opts} analysis/benchmarking/ouija_probit.R --condition {wildcards.cond} --rep {wildcards.rep} --G {wildcards.G} --input_file {input} --output_file {output}"

rule ouija_cloglog:
    input:
        "data/benchmarking/cloglog_synthetic.h5"
    output:
        "data/benchmarking/cloglog/ouija_{cond}_{G}_{rep}.csv"
    shell:
        "Rscript {R_opts} analysis/benchmarking/ouija_cloglog.R --condition {wildcards.cond} --rep {wildcards.rep} --G {wildcards.G} --input_file {input} --output_file {output}"

rule ouija_threshold:
    input:
        "data/benchmarking/threshold_synthetic.h5"
    output:
        "data/benchmarking/threshold/ouija_{cond}_{G}_{rep}.csv"
    shell:
        "Rscript {R_opts} analysis/benchmarking/ouija_threshold.R --condition {wildcards.cond} --rep {wildcards.rep} --G {wildcards.G} --input_file {input} --output_file {output}"

rule ouija_logit:
    input:
        "data/benchmarking/logit_synthetic.h5"
    output:
        "data/benchmarking/logit/ouija_{cond_logit}_{G}_{rep}.csv"
    shell:
        "Rscript {R_opts} analysis/benchmarking/ouija_logit.R --condition {wildcards.cond_logit} --rep {wildcards.rep} --G {wildcards.G} --input_file {input} --output_file {output}"

rule create_benchmark_fig:
    input:
        benchmark_results,
        ouija_csv
    output:
        "figs/fig_benchmark.png",
        "data/benchmarking/summarised_results.csv"
    shell:
        "Rscript {R_opts} analysis/benchmarking/benchmark-fig.R"

# Transient expression --------

rule create_transient_data:
    output:
        transient_synthetic_files
    shell:
        "Rscript {R_opts} analysis/transient/create-synthetic-for-transient.R"

rule transient_ouija:
    input:
        "data/transient/synthetic/transient_sim_{Gt}_{ps}_{rept}.csv"
    output:
        "data/transient/results/ouija_cor_{Gt}_{ps}_{rept}.csv"
    shell:
        "Rscript {R_opts} analysis/transient/ouija_logit.R --rep {wildcards.rept} --G {wildcards.Gt} --prop_switch {wildcards.ps} --input_file {input} --output_file {output}"


rule transient_figure:
    input:
        transient_ouija_files,
    output:
        "figs/fig_transient.png"
    shell:
        "Rscript {R_opts} analysis/transient/transient-fig.R"
