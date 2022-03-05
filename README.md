*This GitHub repo accompanies the manuscript “iNovo479: metabolic modeling provides a roadmap to optimize bioproduct yield from deconstructed lignin aromatics by Novosphingobium aromaticivorans”
Alexandra M. Linz,Yanjun Ma, Samuel Scholz, Daniel R. Noguera, and Timothy J. Donohue*

*iNovo479 is a genome-scale metabolic model intended to assist in designing pipelines for converting depolymerized lignin to bioproducts. This site contains built models, the code and input files to build models, code to run the published simulations on iNovo479, and results from the model runs.*

HOW TO BUILD A MODEL

Go to the folder Model_builds/. Inside you'll find a Python script, build_iNovo.py. You'll need Python 3 installed as well as the packages cobra, ast, and pandas - another package, logging, is optional but helps with error reporting. This script takes three arguments in this order:
1. the path to the file of compound IDs
2. the path to the file of reaction IDs
3. the output path of the resulting model (do not specify a file extension)

For example:
> python build_iNovo.py Input_files/minimal_compounds_2022-03-03.csv Input_files/minimal_reactions_2022-02-11.txt Model_builds/iNovo_base_2022

This will write model in both XML (SBML) format and JSON format. The XML file is used in our scripts and is the more common format. The JSON format is helpful for plotting fluxes in Escher.

The input files and models used in this paper are included in Model_builds/. If you would like to modify iNovo479, you can edit the provided input files and use build_iNovo.py to make a new model version.

HOW TO USE THE MODEL

There are several Python scripts provided in Code/ for running the simulations described in our manuscript. Inside this folder, you'll find more detailed directions for running the scripts.

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
