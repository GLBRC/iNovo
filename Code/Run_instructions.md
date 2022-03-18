CODE FOR RUNNING iNOVO

The scripts in this directory run the analyses described in the iNovo479 paper. They take command line arguments and, for most common use cases, will not require any editing of the scripts. Copies of all model versions have been provided in this directory for ease of use. Scripts will either print results to the command line, or write a file to this directory.

Required Python packages include cobra, pandas, and copy. The package logging is optional, but helpful for reporting errors.

BIOMASS YIELD

The script "calculate_biomass_yield.py" provides the maximum predicted biomass yield in mg dry weight biomass/mmol of substrate. It uses the following arguments in this order:
1. The model to use, including the xml extension
2. The compound ID of the substrate to test

For example:
> python calculate_biomass_yield.py iNovo_base_2022.xml exVA

This script only takes a single substrate. It writes a file called Model_fluxes.csv with fluxes through all reactions in the model in the optimal solution.

The following substrates are currently available for testing in iNovo479:

exC00031 (glucose), expHBA (p-hydroxybenzoic acid), exSA (syringic acid), exS (syringaldehyde), exVA (vanillic acid), exPCA (protocatechuic acid), exV (vanillin), exFA (ferulic acid), exGDK (G-diketone), exSDK (S-diketone), and exSSGGE or exSRGGE or exRSGGE or exRRGGE (GGE stereoisomers)

Note, there are some commented-out lines for removing HPV metabolism or guaiacol metabolism. These will not impact regular runs unless uncommented. If you want to test a substrate with these reactions removed, delete the # symbols from those lines.

CO-METABOLISM

We used dynamic flux balance analysis (dFBA) to investigate co-metabolism of multiple aromatic substrates. In dFBA, each time iNovo479 is run, we adjust the starting concentrations of media components and produced compounds accordingly, then use the new values as the starting concentrations for the next model run. This allows us to model growth curves and metabolism over time.

The script "cometabolism_dFBA.py" takes a model and a list of substrates. It can accomodate any number of substrates so long as they are encoded in the model (see calculate_biomass_yield.py instructions for list). Each substrate is provided at 1 mmol/L.

The arguments, in this order:
1. The model to use, including the xml extension
2. A list of compound IDs to test

Some examples:
> python cometabolism_dFBA.py iNovo_base_2022.xml exVA
> 
> python cometabolism_dFBA.py iNovo_base_2022.xml exC00031 exVA exSA expHBA
> 
> python cometabolism_dFBA.py iNovo_vanAB_2022.xml exVA exSA expHBA

It outputs a file called "dFBA_results.csv" with the concentrations of all tracked compounds at each timepoint. Currently, the default time interval is 30 minutes and this is encoded in the script. The script automatically stops when the total concentration of the provided substrates is zero.
 
PDC PRODUCTION BY AROMATIC:GLUCOSE RATIOS
 
When we grow Novo for PDC production, we often add glucose because the majority of carbon in the aromatic substrate goes to PDC instead of biomass. This means that we see very slow growth of our PDC-producing strain on aromatics only and poor overall PDC production rates when no glucose is added. However, glucose is an additional cost in the process and minimizing its use would be ideal. The script "PDC_dFBA.py" tests PDC yield in g/L/hr on two substrates with user-provided concentrations in mmol/L.

The arguments, in this order:
1. The model with the xml extension
2. Compound ID of substrate 1 (see biomass_yield.py for allowed substrates)
3. Concentration of substrate 1 in mmol/L
4. Compound ID of substrate 2
5. Concentration of substrate 2

For example:
> python PDC_dFBA.py exVA 3.0 exC00031 2.0

This script will print out both the rate of PDC production in g/L/hr and the maximum concentration of PDC achieved. It will also output a file called "PDC_dFBA_results.csv" with the concentrations over time of all tracked metabolites. It will stop automatically when all carbon is gone and will report several conditions during its run. It will tell you if a metabolite uptake rate is limiting, if a metabolite has been consumed, and if biomass production is no longer feasible. Since we know that Novo continues to convert aromatic substrates to PDC even in stationary phase, if biomass production has ceased, it will continue consuming aromatics at the last calculated rate until all aromatics are consumed. Biomass production will cease if the model is unable to meet the constraints of PDC production, the required energy cost of non-growth associated maintenance (NGAM reaction), and biomass production.

BIOPRODUCT DYNAMIC FLUX BALANCE

This script is similar to the previous dynamic flux balance analysis scripts in that it runs iNovo479 iteratively to simulate growth over time, but this script takes a model with "engineered" reactions, a desired bioproduct, and an "overexpression" amount. The model iNovo_engineered_2022.xml has specific additional reactions derived from the KEGG database but not expected to be found N. aromaticivorans, similar to a genetic cloning approach. Since production of a bioproduct generally detracts from biomass yield, we model production of the bioproduct by specifying an amount to produce and meeting that demand before optimizing biomass. We refer to this amount as the overexpression amount, as it is conceptually similar to overexpressing an enzyme to boost turnover of its substrates in vivo. For now, this script only takes 5 mmol/L vanillic acid as its substrate as in our paper, but that can be changed by editing the script.

The arguments, in this order:
1. The model with the xml extension (use the engineered version)
2. The desired bioproduct to test
3. The overexpression amount in mmol product per mmol substrate

For example:
> python bioproduct_dFBA.py iNovo_engineered_2022.xml C00033 0.6

You will likely want to test a lot of combinations. To do this, make a file containing arguments 2 and 3. For example:

C00033 0.5

C00033 0.6

C00033 0.7

C00158 0.3

C00158 0.4


Then use "while read line" in a Unix environment to run everything in the file. For example, if your file was named bioproducts.txt:
> while read line; do python bioproduct_dFBA.py iNovo_engineered_2022.xml $line; done < bioproducts.txt

This script will output several items - the time stopped (it will stop when all substrate is consumed), the production rate in mmol/L/hr of bioproduct, and the rate in g/L/hr of bioproduct. Both the overexpression amount and the g/L/hr are useful for comparing bioproducts. This script will also output a file called "bioproduct_dFBA_results.csv" with the amounts of substrates and products over time in the simulation.

The currently allowed bioproducts to test are C00489 (glutarate), C06098 (zeaxanthin), C02480 (cis-cis muconic acid), C00158 (citrate), C00163 (propanoate), C00084 (acetaldehyde), C00116 (glycerol), C00246 (butanoate), C00823 (1-hexadecanol), C00146 (phenol), C00086 (urea), C00033 (acetate), and C00189 (ethanolamine).

INOVO_FIGURES.R

This R script is not necessarily intended to be run by other users, but it provides all of the R code used to generate the figures in our manuscript. It takes input files from Model_results/ and outputs plots to a directory called Plots_and_Tables/ (not included here). It is intended to be a resource for making figures.
