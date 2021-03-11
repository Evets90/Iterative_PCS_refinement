import os
import sys
import math
from paramagpy import protein, fit, dataparse, metal
import warnings
from Bio import BiopythonWarning
from Bio.PDB import *
warnings.simplefilter('ignore', BiopythonWarning)

# argv and paths
script, prot_file_name, npc_1_f, npc_2_f, npc_3_f = sys.argv
prot_f = os.getcwd() + "/" + prot_file_name
npc_1 = os.getcwd() + "/" + npc_1_f
npc_2 = os.getcwd() + "/" + npc_2_f
npc_3 = os.getcwd() + "/" + npc_3_f

# Universal dictionaries
res_dictionary = {1: 221, 2: 137, 3: 50}
Tag_dictionary = {1: 300, 2: 330, 3: 360}
Ori_dictionary = {1: [215, 216, 217, 218, 219, 220], 2: [129, 130, 131, 132, 133, 134], 3: [88, 89, 90, 91, 92, 93]}
q_fact_dict = {}

# Populate the dictionary with the models
for model in protein.load_pdb(prot_f):
    q_fact_dict[model.id] = ""

def get_qfac(protein_file, npc, res_num, q_fac_dict):
    """Fit and individual tensor to each model in the bundle and save each Q-factor into an universal dictionary"""
    # Load the protein, load the npc
    prot = protein.load_pdb(protein_file)
    rawData = dataparse.read_pcs(npc)
    qfactor_sep = {}

    # Initialize metal instance for search and set the initial position
    mStart = metal.Metal()
    mStart.position = prot[0]['A'][res_num]['CA'].position

    # Loop: for every model fit an individual tensor and store Q-factors and tensor components into a dict
    for model in prot:
        parsedData = prot.parse(rawData, models=model.id)
        [mGuess], _ = fit.svd_gridsearch_fit_metal_from_pcs([mStart], [parsedData], radius=10, points=10)
        [mFit], [data] = fit.nlr_fit_metal_from_pcs([mGuess], [parsedData])
        qfactor_sep[model.id] = fit.qfactor(data)
        # Save in universal dictionary
        if type(q_fact_dict[model.id]) == list:
            q_fac_dict[model.id].append(fit.qfactor(data))
        else:
            q_fac_dict[model.id] = [fit.qfactor(data)]

    minModel, minQfac = sorted(qfactor_sep.items(), key=lambda x: x[1])[0]

def manipulate_dict(dictionary):
    """Creates an average of the values of each key and find the best one (min average). Also writes down everything in
    a file."""
    # Initialize
    ave = {}
    name_q = os.getcwd() + "/Q-factors.txt"
    a = open(name_q, 'a+')
    a.write("Cycle " + str(prot_file_name.split(".")[0]) + "\nMod\tAverage\n")

    # Calculate average and write down individual Q-factors
    for key in dictionary:
        average = sum(dictionary[key]) / len(dictionary[key])
        ave[key] = average
        a.write(str(key+1) + "\t" + str(round(average, 5)) + "\n")

    # Identify best model
    min_model = min(ave.keys(), key=(lambda k: ave[k]))
    a.write("Best model: \t" + str(min_model+1) + "\n\n")

    return min_model

def get_tensor(protein_file, model, npc, res_num, tag_number):
    """Fit a tensor only to the specified model of the bundle pdb and writes down the corresponding components and
    ORI"""
    # Load the protein, load the npc
    prot = protein.load_pdb(protein_file)
    rawData = dataparse.read_pcs(npc)

    # Initialize metal instance for search and set the initial position
    mStart = metal.Metal()
    mStart.position = prot[model]['A'][res_num]['CA'].position

    # Get tensor on single structure
    for mod in prot:
        if mod.id == model:
            parsedData = prot.parse(rawData, models=mod.id)
            [mGuess], [data] = fit.svd_gridsearch_fit_metal_from_pcs([mStart], [parsedData], radius=10, points=10)
            [mFit], [data] = fit.nlr_fit_metal_from_pcs([mGuess], [parsedData])

    Axial = round(mFit.ax * 1E32, 3)
    Rhombicity = round((mFit.rh/mFit.ax), 3)

    # Generate .pcs metal center file
    name = os.path.splitext(protein_file)[0] + "_metal_centers.pcs"
    f = open(name, 'a+')
    f.write(5*" " + str(tag_number) + (14-len(str(Axial)))*" " + str(Axial) + "E+04      " + str(Rhombicity) + 5*" " + str(Tag_dictionary[tag_number]) + "\n")

    # Extract information about the metal center position relative to the protein backbone and write them
    ori = mFit.position
    name_ori_upl = os.path.splitext(protein_file)[0] + "_ORI_UPL.upl"
    name_ori_lol = os.path.splitext(protein_file)[0] + "_ORI_LOL.lol"
    u = open(name_ori_upl, 'a+')
    l = open(name_ori_lol, 'a+')
    for res_num in Ori_dictionary[tag_number]:
        res = prot[model]['A'][res_num]['CA'].position
        res_type = prot[model]['A'][res_num].get_resname()
        d = math.sqrt(((ori[0] - res[0]) ** 2) + ((ori[1] - res[1]) ** 2) + ((ori[2] - res[2]) ** 2))
        d_upl = round((d * 1E10) + 0.5, 2)
        d_lol = round((d * 1E10) - 0.5, 2)
        u.write((3 - len(str(res_num))) * " " + str(res_num) + " " + res_type + "  CA    " + str(Tag_dictionary[tag_number]) + " ORI  A0" + (10 - len(str(d_upl))) * " " + str(d_upl) + "\n")
        l.write((3 - len(str(res_num))) * " " + str(res_num) + " " + res_type + "  CA    " + str(Tag_dictionary[tag_number]) + " ORI  A0" + (10 - len(str(d_lol))) * " " + str(d_lol) + "\n")

def write_best_model(protein_file, model, suffix='_best.pdb'):
    """Just writes down the best model in the same folder with '_best.pdb'"""
    io = PDBIO()
    prot = protein.load_pdb(protein_file)
    io.set_structure(prot[model])
    save_name = os.path.splitext(protein_file)[0] + suffix
    io.save(save_name)

# Set log
log_name = os.getcwd() + "/" + os.path.splitext(prot_file_name)[0] + "_tensor_log.txt"
log = open(log_name, 'w+')
sys.stdout = log

# Start
# One get_qfac for each tag you have
get_qfac(prot_f, npc_1, res_dictionary[1], q_fact_dict)
get_qfac(prot_f, npc_2, res_dictionary[2], q_fact_dict)
get_qfac(prot_f, npc_3, res_dictionary[3], q_fact_dict)

# Leave this unchanged
min_model = manipulate_dict(q_fact_dict)
write_best_model(prot_f, min_model)

# One get_tensor for each tag you have
get_tensor(prot_f, min_model, npc_1, res_dictionary[1], 1)
get_tensor(prot_f, min_model, npc_2, res_dictionary[2], 2)
get_tensor(prot_f, min_model, npc_3, res_dictionary[3], 3)

# Close log
log.close()
