CODE FOR RUNNING iNOVO

The scripts in this directory run the analyses described in the iNovo479 paper. They take command line arguments and, for most cases, will not require any editing of the scripts. Copies of all model versions have been provided in this directory for ease of use. Scripts will either print results to the command line, or write a file to this directory.

Required Python packagese include cobra, pandas, and copy. The package logging is optional, but helpful for reporting errors.

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

We used dynamic flux balance analysis to investigate co-metabolism of multiple aromatic substratesl
