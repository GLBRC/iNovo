###################
# build_model.py
# Copyright 2021, Alexandra Linz, Daniel Noguera, and Timothy Donohue
#
# This script takes files containing compound and reaction info and builds them into a metabolic model
###################

# Load necessary packages
import logging
logging.basicConfig() # COBRApy uses this to report errors
import pandas
import ast # Used for interpreting literal strings
import cobra

# Set file paths
cpd_path = '/Users/Alex/Desktop/iNovo/Model_builds/Input_files/minimal_compounds_2021-10-07.csv'
rxn_path = '/Users/Alex/Desktop/iNovo/Model_builds/Input_files/minimal_reactions_2021-12-04.txt'
output_path = '/Users/Alex/Desktop/iNovo/Model_builds/Models/iNovo'

# Set up model object
model = cobra.Model('Novo')

# Read in the file of compound IDs
# Includes compound ID, chemical formula, name of compound, and compartment
# Compartment options are c0 (cytosol), p0 (periplasm), and e (extracellular space)
compounds = pandas.read_csv(cpd_path, sep = ",", header = None)
compounds.columns = ["cpdID", "Formula", "Name", "compartment"]
compounds["cpdID"] = compounds["cpdID"].str.strip() # Removes leading and trailing white spaces in compound IDs
# Where do white spaces come from? I copy/pasted some of this data into Excel and may have accidentally
# brought some white space with it. Also, Excel does weird things to data

# Read in file of reaction IDs
# Includes reaction ID, reversibility of reaction (TRUE/FALSE), compounds in reaction (negative indicates consumed, positive indicates produced),
# name of reaction, genes encoding enzymes for the reaction (can be AND or OR), and where in the cell it takes place
# The reaction compartment is more for my own notes than the model - compartments are encoded in compound IDs
# Also, some of my reaction "compartments" are the membranes between compartments
reactions = pandas.read_csv(rxn_path, sep = "\t", header = None)
reactions.columns = ["rxnID", "reversibility", "cpds", "names", "genes", "compartment"]
reactions["rxnID"] = reactions["rxnID"].str.strip() # Similar to above, str.strip() removes leading and trailing white spaces
reactions["names"] = reactions["names"].str.strip()
reactions["cpds"] = reactions["cpds"].str.strip()

# Add compounds to the model
# If adding a compound fails, report which line

for index, row in compounds.iterrows(): # iterrows() is a function from pandas specifically for looping through dataframes
    try:
        cpd_to_add = cobra.Metabolite(row["cpdID"], name = row["Name"], formula = row["Formula"], compartment = row["compartment"])
        model.add_metabolites([cpd_to_add])
    except ValueError:
        print("Problem line:", index)

# Add reactions to the model in a similar method to compounds
# But accounting for reversibility as well. Although I have that encoded as TRUE/FALSE, COBRApy encodes it as
# a negative lower bound (reversible) or a zero or positive lower bound (irreversibile)
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
        if row['reversibility'] == True:
            model.reactions.get_by_id(row["rxnID"]).lower_bound = -1000
        else:
            model.reactions.get_by_id(row["rxnID"]).lower_bound = 0.0

    except ValueError:
        print("Problem line:", index)

model.objective = model.reactions.get_by_id("biomass")

# Optional: use these statements to check the mass balance of specific reactions
#print("ppA007: ", model.reactions.get_by_id("ppA007").check_mass_balance())
#print("altA018: ", model.reactions.get_by_id("altA018").check_mass_balance())

# Write output - the following scripts use SBML format, but I also use JSON for plotting flux in Escher
output_path_json = output_path + '.json'
output_path_xml = output_path + '.xml'
cobra.io.write_sbml_model(model, output_path_xml)
cobra.io.save_json_model(model, output_path_json)
