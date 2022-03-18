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

	-Run_instructions.MD	#Instructions for running the code in this directory
	
	-PDC_dFBA.py		#Use dFBA to test aromatic:glucose ratios for maximizing PDC yield
	
	-bioproduct_dFBA.py	#Use dFBA to assess various bioproduct yields
	
	-calculate_biomass_yield.py	#Determine biomass yield from various substrates
	
	-cometabolism_dFBA.py	#Use dFBA to determine how multiple substrates are consumed
	
	-iNovo_figures.R	#Generate figures from the manuscript
	
	-iNovo.xml files	#Included copies of the models in Model_builds/Models/ for ease of use in these scripts

-Model_builds

	-build_iNovo.py		#Script for building models from input files

	-Input_files
	
		-hypothetical_demethylation_minimal_reactions_2022-03-02.txt	#Hypothetical demethylation with no energy gain
		
		-minimal_compounds_2022-03-03.csv		#Shared compound IDs file for all models
		
		-minimal_reactions_2022-02-11.txt		#Reaction IDs in the base iNovo model
		
		-vanAB_minimal_reactions_2022-03-03.txt		#Novo’s demethylation replaced with P. putida’s 
		
		-engineered_minimal_reactions_2022-03-02.txt	#Additional reactions added to allow new bioproducts
	-Models
		-iNovo_hypo_demeth_2022.json/xml	#Hypothetical demethylation with no energy gain
		
		-iNovo_base_2022.json/xml	#Base model
		
		-iNovo_vanAB_2022.json/xml	#Novo’s demethylation replaced with P. putida’s 
		
		-iNovo_engineered_2022.json/xml	#Base model + additional reactions to allow new bioproducts
		
-Model_results

	-SA_VA_pHBA_dFBA_results-2022.csv	#Cometabolism of G, S, and H aromatics
	
	-glu_SA_VA_pHBA_dFBA_results-2022.csv	#Cometabolism of G, S, and H aromatics and glucose
	
	-hypothetical_SA_VA_pHBA_dFBA_results-2022.csv	#Cometabolism of G, S, and H aromatics with hypothetical demethylation reaction
	
	-vanAB_SA_VA_pHBA_dFBA_results-2022.csv	#Cometabolism of G, S, and H aromatics with P. putida demethylation
	
	-aromatic_glucose_ratios-2022.csv	#PDC producion rate results of varying aromatic:glucose ratios
	
	-biomass_yields.csv			#Biomass yield predictions and experimental yields for aromatic substrates
