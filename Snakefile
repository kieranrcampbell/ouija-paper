"""
Reproducible analysis for Ouija paper
"""


# Configuration --------------

# Benchmarking variables
index = [1,2,3] # Three datasets
rep = list(range(1,51)) # 50 repetitions
num_additional_markers = [1,5,10,20,50,100,500,1000] # Number of additional markers
nam_index = range(1,len(num_additional_markers)+1)

output_csv = expand("data/marker-vs-txome/ouija_fits/ouija_pseudotime_{i}_{r}_{nam}.csv",
			i = index, r = rep, nam = nam_index)

output_de = expand("data/marker-vs-txome/subset_de_fits/sde_fit_{i}_{r}_{nam}.csv",
			i = index, r = rep, nam = nam_index)




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

rule construct_scesets:
    input:
        "data/scesets/hsc-sce.rds",
        "data/raw/waterfall_data.xlsx",
        "data/raw/hsmm.rds"
    output:
        "data/scesets/scesets.rds"
    shell:
        "Rscript {R_opts} scripts/create_scesets.R"


# Benchmarking --------------------

rule fit_global_pseudotime:
    input:
        "data/scesets/scesets.rds"
    output:
        "data/scesets/scesets_with_pseudotime.rds"
    shell:
        "Rscript {R_opts} analysis/marker-vs-txome/01_fit_global_pseudotime.R"

rule markers_to_use:
    input:
        "data/scesets/scesets_with_pseudotime.rds"
    output:
        "data/marker-vs-txome/markers_to_use.rds"
    shell:
        "Rscript {R_opts} analysis/marker-vs-txome/02_decide_marker_list.R"


rule fit_pseudotimes_across_marker_subsets:
    input:
        "data/scesets/scesets_with_pseudotime.rds",
        "data/marker-vs-txome/markers_to_use.rds"
    output:
        "data/scesets/scesets_with_marker_pseudotime.rds",
        "data/marker-vs-txome/pseudotimes_across_marker_subset_array.rds"
    shell:
        "Rscript {R_opts} analysis/marker-vs-txome/03_fit_pseudotimes_across_marker_subsets.R"

rule fit_ouija_pseudotimes:
    input:
        "data/scesets/scesets_with_pseudotime.rds",
        "data/marker-vs-txome/markers_to_use.rds"
    output:
        "data/marker-vs-txome/ouija_fits/ouija_pseudotime_{i}_{r}_{nam}.csv"
    shell:
        "Rscript {R_opts} analysis/marker-vs-txome/04_fit_ouija_pseudotimes.R {wildcards.i} {wildcards.r} {wildcards.nam}"

rule fit_subset_de:
	input:
		"data/scesets/scesets_with_pseudotime.rds",
		"data/marker-vs-txome/pseudotimes_across_marker_subset_array.rds",
		"data/marker-vs-txome/ouija_fits/ouija_pseudotime_{i}_{r}_{nam}.csv"
	output:
		"data/marker-vs-txome/subset_de_fits/sde_fit_{i}_{r}_{nam}.csv"
	shell:
		"Rscript {R_opts} analysis/marker-vs-txome/05_fit_subset_de.R {wildcards.i} {wildcards.r} {wildcards.nam}"

rule fit_txwide_marker_de:
	input:
		"data/scesets/scesets_with_marker_pseudotime.rds"
	output:
		"data/marker-vs-txome/txwide_marker_de.csv"
	shell:
		"Rscript {R_opts} analysis/marker-vs-txome/06_fit_marker_txome_de.R"


rule tidy_de:
	input:
		output_de,
		"data/marker-vs-txome/txwide_marker_de.csv"
	output:
		"data/marker-vs-txome/tidy_de.csv"
	shell:
		"Rscript {R_opts} analysis/marker-vs-txome/07_tidy_de.R"

rule make_de_plot:
    input:
        "data/scesets/scesets.rds",
        "data/marker-vs-txome/tidy_de.csv"
    output:
        "figs/marker-vs-txome-de-results.png"
    shell:
        "Rscript {R_opts} analysis/marker-vs-txome/08_graph_all_de.R"

rule graph_pseudotime_subsets:
    input:
        "data/scesets/scesets_with_marker_pseudotime.rds",
        "data/marker-vs-txome/pseudotimes_across_marker_subset_array.rds"
    output:
        "figs/marker-vs-txome-pseudotime-subsets.png"
    shell:
        "Rscript {R_opts} analysis/marker-vs-txome/09_graph_pseudotime_subsets.R"


