###################
# calculate_biomass_yield.py
# Copyright 2022, Alexandra Linz, Daniel Noguera, and Timothy Donohue
#
# This script optimizes model flux for biomass from one or more carbon substrates
# Its output is a dataframe that can be plotted with a separate script
###################

# Import packages
import sys
import cobra
from cobra.flux_analysis.loopless import loopless_solution
import logging
import pandas
import copy
logging.basicConfig()

# There's a small (1e-7) error rate in the standard optimization method
# So since some of my fluxes are smaller than that, I needed to use the slower but more accurate glpk_exact
cobra_config = cobra.Configuration()
cobra_config.solver = "glpk_exact"

# Set paths
model_path = sys.argv[1]
output_path = 'Model_fluxes.csv' # Must be a csv file

# Note: there's a warning "Solver status infeasible" that may appear when running this script.
# This means that a model constraint cannot be met
# I have a constraint set on the S-type aromatic pathway, so running the model on G-type aromatics only will result in this error. No action needs to be taken in these cases.
# It may also indicate that no loopless solution could be found. You should take a closer look at these cases.

###############
substrates = {"exC00031": [0.0], "expHBA": [0.0], "exSA": [0.0], "exS": [0.0], "exVA": [0.0], "exPCA": [0.0], "exV": [0.0], "exFA": [0.0], "exGDK": [0.0], "exSDK": [0.0], "exSSGGE": [0.0], "exSRGGE": [0.0], "exRSGGE": [0.0], "exRRGGE": [0.0]} 

# Update with desired substrate
# Leave substrate concentration at 1 unless you want to do some conversions with the final biomass value
substrates[sys.argv[2]] = [1.0] 	

media_components = {"exC14818": [45.54], "exC00014": [1.0], "exC00009": [26.1], "exC00059": [8.]}
enviro = {"C00282": [10], "exC00001": [10], "exC00007": [10]}
outfluxes = {"C00067": [0], "C00058": [0], "C00033": [0], "C00010": [1], "C00162": [1], "C00010": [1], "C00132": [0], "C00054": [0], "C00011": [0], "C05198": [0], "C04425": [0], "C00266": [0], "C00153": [0]}

##############

# Kinetic parameters in mmol/L per min - these are estimates from related bacteria in the literature and not experimentally verified
nonaromatic_substrates = ["exC00031", "exC00243", "exC00033", "exC00022", "exC00095", "exC00181", "exC00208", "exC00185"]
def get_rate(cpd_ID):
    global Vm
    global Ks
    # glucose
    if cpd_ID in nonaromatic_substrates:
        Vm = 0.5
        Ks = 0.139
        Ki = 0.139
    # aromatics
    # S-type
    elif cpd_ID == "exSA" or cpd_ID == "exS" or cpd_ID == "exSDK":
        Vm = 0.582
        Ks = 0.05
        Ki = 0.05
    # H-type
    elif cpd_ID == "expHBA" or cpd_ID == "exPCA" or cpd_ID == "exC00633" or cpd_ID == "exC00180"or cpd_ID == "exC00156":
        Vm = 0.902
        Ks = 0.05
        Ki = 0.05
    # G-type
    elif cpd_ID == "exVA" or cpd_ID == "exV" or cpd_ID == "exFA" or cpd_ID == "exGDK" or cpd_ID == "exSRGGE" or cpd_ID == "exSSGGE" or cpd_ID == "exRRGGE" or cpd_ID == "exRSGGE":
        Vm = 0.569
        Ks = 0.1
        Ki = 0.1
    # ammonia
    elif cpd_ID == "exC00014":
        Vm = 0.5
        Ks = 0.1
        Ki = 0.1
    # phosphate
    elif cpd_ID == "exC00009":
        Vm = 0.060
        Ks = 0.002
        Ki = 0.1
    # sulfate
    elif cpd_ID == "exC00059":
        Vm = 0.0017
        Ks = 0.003
        Ki = 0.1
    # iron
    elif cpd_ID == "exC14818":
        Vm = 0.0017
        Ks = 0.003
        Ki = 0.1
    else:
        Vm = 0
        Ks = 0
        Ki = 0.1
    return Vm, Ks, Ki



##############
# Load and set up the model
Novo_model = cobra.io.read_sbml_model(model_path)
media_components.update(substrates) # Add carbon sources to the basic medium recipe

# Add the SA constraint
SA_flux = Novo_model.problem.Constraint(
    Novo_model.reactions.A031.flux_expression - Novo_model.reactions.A015.flux_expression * 0.15,
    lb=0,
    ub=0, name = 'SA_flux')
Novo_model.add_cons_vars(SA_flux)


# Combine all tracked items into the tracking dictionary
tracking = copy.deepcopy(media_components)
tracking.update(enviro)
tracking.update(outfluxes)


#############
# Run the model

track_rates = copy.deepcopy(substrates)
track_rates.update(media_components)

for metabolite in tracking:

    if metabolite in substrates or metabolite in media_components:
        if metabolite in substrates:
            # Set the rate based on the kinetic parameters above
            r = get_rate(metabolite)[0] * (tracking[metabolite][0] / ((tracking[metabolite][0] + get_rate(metabolite)[1]) * ( 1 + tracking[metabolite][0] / get_rate(metabolite)[2])))/1000
        else:
            r = get_rate(metabolite)[0] * (tracking[metabolite][0] / ((tracking[metabolite][0] + get_rate(metabolite)[1]) ))/1000

        exec ("Novo_model.reactions.EX_" + metabolite + ".lower_bound = -1 * r")
        exec ("Novo_model.reactions.EX_" + metabolite + ".upper_bound = 1 * r")
        exec ("track_rates[metabolite] = Novo_model.reactions.EX_" + metabolite + ".upper_bound")

# Run the optimization

opt = Novo_model.optimize()
solution = loopless_solution(Novo_model)
out = solution.fluxes

solution.fluxes.to_csv(output_path)

# We don't need to open the output file to get the biomass yield, though - print that out right here
# Biomass reaction take input in mmol and output in g, so just need to convert g to mg in the final value

carbon = list(substrates.keys())

carbon_flux = 0
for k in range(0,len(carbon)):
    exec("carbon_flux = carbon_flux + out[\"EX_" + carbon[k] + "\"]")

print("Biomass yield, mgDW/mmol: ")
print(out["biomass"] / (carbon_flux) * -1000)

