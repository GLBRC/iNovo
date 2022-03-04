###################
# build_iNovo.py
# Copyright 2022, Alexandra Linz, Daniel Noguera, and Timothy Donohue
#
# This script takes files containing compound and reaction info and builds them into a metabolic model
###################

# Load necessary packages
import sys
import logging
logging.basicConfig() # COBRApy uses this to report errors
import pandas
import ast # Used for interpreting literal strings
import cobra

# Set file paths
cpd_path = sys.argv[1]	# 1st argument is the path to the file of compound IDs
rxn_path = sys.argv[2]	# 2nd argument is the path to the file of reaction IDs
output_path = sys.argv[3]	# 3rd argument is the path and name of the model to output

# Set up model object
model = cobra.Model('Novo')

# Read in the file of compound IDs
# Includes compound ID, chemical formula, name of compound, and compartment
# Compartment options are c0 (cytosol), p0 (periplasm), and e (extracellular space)
compounds = pandas.read_csv(cpd_path, sep = ",", header = None)
compounds.columns = ["cpdID", "Formula", "Name", "compartment"]
compounds["cpdID"] = compounds["cpdID"].str.strip() # Removes leading and trailing white spaces in compound IDs
# Where do white spaces come from? I copy/pasted some of this data into Excel and may have accidentally brought some white space with it.

# Read in file of reaction IDs
# Includes reaction ID, reversibility of reaction (TRUE/FALSE), compounds in reaction (negative indicates consumed, positive indicates produced),
# name of reaction, genes encoding enzymes for the reaction (can be AND or OR), and SBO annotation

reactions = pandas.read_csv(rxn_path, sep = "\t", header = None)
reactions.columns = ["rxnID", "reversibility", "cpds", "names", "genes", "sbo"]
reactions["rxnID"] = reactions["rxnID"].str.strip() # Similar to above, str.strip() removes leading and trailing white spaces
reactions["names"] = reactions["names"].str.strip()
reactions["cpds"] = reactions["cpds"].str.strip()

# Add compounds to the model
# If adding a compound fails, report which line

for index, row in compounds.iterrows(): # iterrows() is a function from pandas specifically for looping through dataframes
    try:
        cpd_to_add = cobra.Metabolite(row["cpdID"], name = row["Name"], formula = row["Formula"], compartment = row["compartment"])
        model.add_metabolites([cpd_to_add])
        model.metabolites.get_by_id(row["cpdID"]).annotation["kegg.compound"] = row["cpdID"]
    except ValueError:
        print("Problem line:", index)

for m in model.metabolites:
    m.annotation["sbo"] = "SBO:0000247"

# Add reactions to the model in a similar method to compounds
# But accounting for reversibility as well. Although I have that encoded as TRUE/FALSE, COBRApy encodes it as a negative lower bound (reversible) or a zero or positive lower bound (irreversible)
# The bound limit of 1000 is arbitrary - I'm not incorporating any data right now to constrain those bounds

for index, row in reactions.iterrows():
    if row['reversibility'] == True:
        lb = -1000
    else:
        lb = 0.0
    try:
        rxn_to_add = cobra.Reaction(row["rxnID"], lower_bound = 0.0, upper_bound=1000., name=row["names"])
        rxn_to_add.gene_reaction_rule = '( ' + str(row["genes"]) + ' )' # Convert gene data to COBRApy format and add with specific function
        model.add_reactions([rxn_to_add])
        model.reactions.get_by_id(row["rxnID"]).add_metabolites(ast.literal_eval("{" + row["cpds"] + "}"))
        model.reactions.get_by_id(row["rxnID"]).annotation["kegg.reaction"] = row["rxnID"]
        model.reactions.get_by_id(row["rxnID"]).annotation["sbo"] = row["sbo"]
        if row['reversibility'] == True:
            model.reactions.get_by_id(row["rxnID"]).lower_bound = -1000.
        else:
            model.reactions.get_by_id(row["rxnID"]).lower_bound = 0.0

    except ValueError:
        print("Problem line:", index)

for g in model.genes:
    g.annotation["sbo"] = "SBO:0000243"
    g.annotation["refseq"] = g.id

model.objective = model.reactions.get_by_id("biomass")

# Optional: use these statements to check the mass balance of specific reactions
#print("engR00449: ", model.reactions.get_by_id("engR00449").check_mass_balance())
#print("engR02273: ", model.reactions.get_by_id("engR02273").check_mass_balance())

# Add exchanges for the following substrates, but set them to zero for now.

substrates = {"exC00031": [1.0], "expHBA": [1.0], "exSA": [1.0], "exS": [1.0], "exVA": [1.0], "exPCA": [1.0], "exV": [1.0], "exFA": [1.0], "exGDK": [1.0], "exSDK": [1.0], "exSSGGE": [1.0], "exSRGGE": [1.0], "exRSGGE": [1.0], "exRRGGE": [1.0]} 
media_components = {"exC14818": [45.54], "exC00014": [1.0], "exC00009": [26.1], "exC00059": [8.]}
# These are exchange reactions that are far in excess of others due to diffusion
enviro = {"C00282": [10], "exC00001": [10], "exC00007": [10]}
# These are things the model needs to be allowed to output or it will break - your run may not need all of these
outfluxes = {"C00067": [0], "C00058": [0], "C00033": [0], "C00010": [1], "C00162": [1], "C00010": [1], "C00132": [0], "C00054": [0], "C00011": [0], "C05198": [0], "C04425": [0], "C00266": [0], "C00153": [0]}


# Add media components, cofactors, virtually unlimited influxes, and outfluxes to the model as exchange reactions
for item in media_components.keys():

    new_exchange = "EX_" + item
    model.add_boundary(model.metabolites.get_by_id(item), type = "exchange",  ub = 1000., reaction_id=new_exchange)

for item in substrates.keys():

    new_exchange = "EX_" + item
    model.add_boundary(model.metabolites.get_by_id(item), type = "exchange",  ub = 0.0, reaction_id=new_exchange)

for item in enviro:

    new_enviro = "EX_" + item
    model.add_boundary(model.metabolites.get_by_id(item), ub = 1000., type = "exchange", reaction_id=new_enviro)

for item in outfluxes:
    new_outflux = "DM_" + item
    model.add_boundary(model.metabolites.get_by_id(item), ub = 1000., type = "demand", reaction_id=new_outflux)

# Set the model objective
objective = model.reactions.get_by_id("biomass")
model.objective = objective

# Set the non-growth associated maintenance requirement
model.reactions.get_by_id("NGAM").upper_bound = 0.00004
model.reactions.get_by_id("NGAM").lower_bound = 0.00004


# Write output - the following scripts use SBML format, but I also use JSON for plotting flux in Escher
output_path_json = output_path + '.json'
output_path_xml = output_path + '.xml'
cobra.io.write_sbml_model(model, output_path_xml)
cobra.io.save_json_model(model, output_path_json)
