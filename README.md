*This GitHub repo accompanies the manuscript “iNovo479: metabolic modeling provides a roadmap to optimize bioproduct yield from deconstructed lignin aromatics by Novosphingobium aromaticivorans”
Alexandra M. Linz,Yanjun Ma, Samuel Scholz, Daniel R. Noguera, and Timothy J. Donohue*

*iNovo479 is a genome-scale metabolic model intended to assist in designing pipelines for converting depolymerized lignin to bioproducts. This site contains built models, the code and input files to build models, code to run the published simulations on iNovo479, and results from the model runs.*

HOW TO BUILD A MODEL

There are three possible model builds provide here allowing alternative aromatic demethylation pathways. Each model is built using a file of compound IDs and a file of reaction IDs. The compound IDs file is the same for all three models, “Model_builds/Input_files/minimal_compounds_2021-10-07.csv”. The reaction IDs are different for each build and are described below. To build a model, edit the input and output file paths in “Code/build_model.py”, then run this file. It will output the model object in “Model_builds/Models”.

HOW TO USE THE MODEL

There are several Python scripts provided in “Code/“ for running the simulations described in our manuscript:

biomass_yield.py - Edit the top of the script with your desired model and carbon substrate. The substrate must be a compound ID from the model with a transport reaction encoded. This will report the predicted biomass yield in mgDW/mmol substrate. This script can also handle multiple substrates.

run_standard_dFBA.py - This performs dynamic flux balance analysis (dFBA) with no fancy changes to show how substrates are consumed and biomass is generated over time. It’s used in the manuscript to investigate co-metabolism of multiple aromatic substrates. Edit the top of the script with your desired model, carbon substrates, number of timepoints, and intervals between timepoints.

run_PDC_dFBA.py - This script is specific to the PDC producing strain, although it can be used with any base model. This script is configured to take command line arguments for two substrates and their concentrations in mmol/L, so you would run it like:

python ./run_PDC_dFBA.py exVA 2.0 exC00031 3.0

You can also use the following syntax with a file of substrate inputs to run multiple versions of the model at once:

while read line; do python ./run_PDC_dFBA.py $line; done < substrate_inputs.txt

Make sure you edit the timepoint interval and number of timepoints in the script. This script outputs the rate of PDC production in g/L/hr.

iNovo_figures.R - This script takes output from the previous scripts and generates the figures in the manuscript. It is an R script and depends on the packages tidyverse, cowplot, and reshape2.

CONTENTS

-Code

	-biomass_yield.py	#Predicts biomass yield for given substrates
	
	-build_model.py	#Takes a file of compound IDs and a file of reaction IDs and builds a model object
	
	-run_standard_dFBA.py	#Runs dFBA on a given model and set of substrates with no modifications
	
	-run_PDC_dFBA.py	#Run a version of dFBA that simulates PDC production
	
	-iNovo_figures.R	#Generates the figures in the manuscript

-Model_builds

	-Input_files
	
		-hypothetical_demethylation_minimal_reactions_2021-12-09.txt	#Hypothetical demethylation with no energy gain
		
		-minimal_compounds_2021-10-07.csv	#Shared compound IDs file for all models
		
		-minimal_reactions_2021-12-04.txt	#Reaction IDs in the base iNovo model
		
		-vanAB_minimal_reactions_2021-12-09.txt	#Novo’s demethylation replaced with P. putida’s 
	-Models
		-hypothetical_demethylation_iNovo.json/xml	#Hypothetical demethylation with no energy gain
		
		-iNovo.json/xml	#Base model
		
		-vanAB_iNovo.json/xml	#Novo’s demethylation replaced with P. putida’s 
		
-Model_results

	-glu_SA_VA_pHBA_dFBA_results.csv	#Co-metabolism of these substrates from run_standard_dFBA.py
	
	-SA_VA_pHBA_dFBA_results.csv	#See above
	
	-aromatic_glucose_ratios.csv	#Aggregated results from run_PDC_dFBA.py
	
-Plots_and_Tables

	-input plots and tables from iNovo_figures.R and Excel. Used in the manuscript.
	
-Supplemental_Files

	-Supplemental files included with the manuscript
