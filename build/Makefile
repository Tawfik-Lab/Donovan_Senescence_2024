all: figures

download: datasets/North_2019/downloads/North_2019_supplementary_data.zip datasets/GSE249746/downloads/GSE249746_Expression_matrix_raw_counts.csv.gz datasets/GSE155622/downloads/GSE155622_raw_UMI_counts_1.txt.gz datasets/GSE154659/downloads/GSE154659_C57_Raw_counts.RDS.gz datasets/GSE154659/downloads/GSE154659_Atf3_WT_KO_Raw_counts.RDS.gz

datasets/North_2019/downloads/North_2019_supplementary_data.zip datasets/GSE249746/downloads/GSE249746_Expression_matrix_raw_counts.csv.gz datasets/GSE155622/downloads/GSE155622_raw_UMI_counts_1.txt.gz datasets/GSE154659/downloads/GSE154659_C57_Raw_counts.RDS.gz datasets/GSE154659/downloads/GSE154659_Atf3_WT_KO_Raw_counts.RDS.gz: 
	cd ../src/lib && bash fetch_dataset_files.sh


datasets/GSE154659/downloads/GSE154659_Atf3_WT_KO_Raw_counts.mtx: datasets/GSE154659/downloads/GSE154659_Atf3_WT_KO_Raw_counts.RDS.gz
	gunzip -k -f datasets/GSE154659/downloads/GSE154659_Atf3_WT_KO_Raw_counts.RDS.gz
	Rscript ../src/lib/convert_rds_to_mtx.R datasets/GSE154659/downloads/GSE154659_Atf3_WT_KO_Raw_counts.RDS
	
datasets/GSE154659/downloads/GSE154659_C57_Raw_counts.mtx: datasets/GSE154659/downloads/GSE154659_C57_Raw_counts.RDS.gz
	gunzip -k -f datasets/GSE154659/downloads/GSE154659_C57_Raw_counts.RDS.gz
	Rscript ../src/lib/convert_rds_to_mtx.R datasets/GSE154659/downloads/GSE154659_C57_Raw_counts.RDS

build_h5ad: datasets/GSE154659/GSE154659.h5ad datasets/GSE249746/GSE249746.h5ad datasets/GSE155622/GSE155622.h5ad datasets/North_2019/North_2019.h5ad

datasets/GSE154659/GSE154659.h5ad: datasets/GSE154659/downloads/GSE154659_C57_Raw_counts.mtx datasets/GSE154659/downloads/GSE154659_Atf3_WT_KO_Raw_counts.mtx
	cd ../src/notebooks && PYTHONPATH=../lib jupyter execute --timeout=-1 GSE154659_dataset_builder.ipynb

datasets/GSE249746/GSE249746.h5ad: datasets/GSE249746/downloads/GSE249746_Expression_matrix_raw_counts.csv.gz
	cd ../src/notebooks && PYTHONPATH=../lib jupyter execute --timeout=-1 GSE249746_dataset_builder.ipynb

datasets/GSE155622/GSE155622.h5ad: datasets/GSE155622/downloads/GSE155622_raw_UMI_counts_1.txt.gz
	cd ../src/notebooks && PYTHONPATH=../lib jupyter execute --timeout=-1 GSE155622_dataset_builder.ipynb

datasets/North_2019/North_2019.h5ad: datasets/North_2019/downloads/North_2019_supplementary_data.zip
	cd ../src/notebooks && PYTHONPATH=../lib jupyter execute --timeout=-1 North_2019_dataset_builder.ipynb

pre_figures:
	mkdir -p figures
	mkdir -p supplement


figures/dotplot_renthal_senescence.png figures/dotplot_renthal_atf3_split_senescence.png figures/renthal_senmayo_neuron_heatmap.svg figures/renthal_senmayo_neuron_heatmap.png supplement/renthal_senescence.xlsx supplement/renthal_differential_expression.xlsx supplement/renthal_atf3_split_senescence.xlsx supplement/renthal_senmayo_neuron_scores_stats.xlsx: build_h5ad pre_figures
	cd ../src/notebooks && PYTHONPATH=../lib jupyter execute --timeout=-1 GSE154659_figures_supplement.ipynb

figures/dotplot_wang_senescence.png supplement/wang_senescence.xlsx supplement/wang_differential_expression.xlsx: build_h5ad pre_figures
	cd ../src/notebooks && PYTHONPATH=../lib jupyter execute --timeout=-1 GSE155622_figures_supplement.ipynb

supplement/human_single_soma_sen_signatures.xlsx: build_h5ad pre_figures
	cd ../src/notebooks && PYTHONPATH=../lib jupyter execute --timeout=-1 GSE249746_figures_supplement.ipynb

supplement/north_2019_tests.xlsx: build_h5ad pre_figures
	cd ../src/notebooks && PYTHONPATH=../lib jupyter execute --timeout=-1 North_2019_figures_supplement.ipynb

figures: supplement/north_2019_tests.xlsx supplement/human_single_soma_sen_signatures.xlsx supplement/wang_differential_expression.xlsx supplement/renthal_differential_expression.xlsx

clean:
	rm -r datasets
	rm -r figures
	rm -r supplement
