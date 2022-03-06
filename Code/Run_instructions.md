CODE FOR RUNNING iNOVO

The scripts in this directory run the analyses described in the iNovo479 paper. They take command line arguments and, for most cases, will not require any editing of the scripts. Copies of all model versions have been provided in this directory for ease of use. Scripts will either print results to the command line, or write a file to this directory.

Required Python packages include cobra, pandas, and copy. The package logging is optional, but helpful for reporting errors.

BIOMASS YIELD

The script "calculate_biomass_yield.py" provides the maximum predicted biomass yield in mg dry weight biomass/mmol of substrate. It uses the following arguments:
1. The model to use, including the xml extension
2. The compound ID of the substrate to test

For example:
> python calculate_biomass_yield.py iNovo_base_2022.xml exVA

This script only takes a single substrate. It writes a file called Model_fluxes.csv with fluxes through all reactions in the model in the optimal solution.

The following substrates are currently available for testing in iNovo479:

exC00031 (glucose), expHBA (p-hydroxybenzoic acid), exSA (syringic acid), exS (syringaldehyde), exVA (vanillic acid), exPCA (protocatechuic acid), exV (vanillin), exFA (ferulic acid), exGDK (G-diketone), exSDK (S-diketone), exSSGGE or exSRGGE or exRSGGE or exRRGGE (GGE stereoisomers)

CO-METABOLISM

We used dynamic flux balance analysis (dFBA) to investigate co-metabolism of multiple aromatic substrates. In dFBA, each time iNovo479 is run, we adjust the starting concentrations of media components and produced compounds accordingly, then use the new values as the starting concentrations for the next model run. This allows us to model growth curves and metabolism over time.

The script "cometabolism_dFBA.py" takes a model and a list of substrates. It can accomodate any number of substrates so long as they are encoded in the model (see biomass_yield.py instructions for list). Each substrate is provided at 1 mmol/L.

The arguments:
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
 
When we grow Novo for PDC production, we often add glucose in a 1:1 ratio to the aromatic substrate because the majority of carbon in the aromatic substrate goes to PDC instead of biomass. This means that we see very slow growth of our PDC-producing strain on aromatics only and poor overall PDC production rates. However, glucose is an additional cost in the process and minimizing its use would be ideal. The script "PDC_dFBA.py" tests PDC yield in g/L/hr on two substrates with user-provided concentrations in mmol/L.

The arguments:
1. The model with the xml extension
2. Compound ID of substrate #1 (see biomass_yield.py for allowed substrates)
3. Concentration of substrate #1 in mmol/L
4. Compound ID of substrate #2
5. Concentration of substrate #2

For example:
> python PDC_dFBA.py exVA 3.0 exC00031 2.0

This script will output both the rate of PDC production in g/L/hr and the maximum concenration of PDC achieved. It will also output a file called "PDC_dFBA_results.csv" with the concentrations over time of all tracked metabolites. It will stop automatically when all carbon is gone and will report several conditions during its run. It will tell you if a metabolite uptake rate is limiting, if a metabolite has been consumed, and if biomass production is no longer feasible. Since we know that Novo continues to convert aromatic substrates to PDC even in stationary phase, if biomass production has ceased, it will consuming aromatics at the last calculate rate until all aromatics are consumed.
