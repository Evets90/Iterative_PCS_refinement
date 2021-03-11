# INTRODUCTION

The following protocol capture is related to the paper "An automated iterative approach for protein structure refinement using pseudocontact shifts" by Stefano Cucuzza et al. Before attempting to replicate it, it is strongly encouraged to read the paper.

---

# DEPENDENCIES

Before you start, make sure you have installed the following dependencies:

* cyana - suggested version: 3.98.12 - obtain a valid license from [here](https://www.las.jp/english/products/cyana.html)
* python3 - suggested version: 3.7 - download from [here](https://www.python.org/downloads/)
* biopython (python 3 module) - download from [here](https://biopython.org/wiki/Download)
* paramagpy (python 3 module) - download from [here](https://pypi.org/project/paramagpy/)

It is also strongly encouraged to run this protocol on a cluster or high multi-core processor since it will take quite some time. Parallelization is easy to establish simply changing the variable number of processor `nproc` (see below).

---

# PROCESS OVERVIEW

The final goal of this process is to perform an iterative structure refinement on a *model.pdb* protein using experimental PCS through two main scripts, *process.py* and *CALC.cya*. In general the process will follow these steps:

1. The *CALC.cya* cyana macro is be called directly inside cyana and acts as the main script by calling all the relevant functions.
2. The initial *model.pdb* is loaded, renamed *1.pdb* (for incremental cycles) and *process.py* is used on it.
	1. For each available tag, a tensor is fitted on each individual model inside the bundle *model.pdb* in order to find the corresponding Q-factor. (Note: you can also start from a single structure, like in the example below)
	2. In each individual model of the bundle, Q-factors for each tag are averaged to find the minimum averaged Q-factor, representing the best structure by Q-factors, and this particular structure is saved as a *_best.pdb*.
	3. For each tag, a tensor is calculated on the best structure and the corresponding axial component and rhombicity are saved in a *_metal_centers.pcs* file while distances to 6 nearby CA atoms with a set tolerance are saved in *_ORI_LOL.lol* and *_ORI_UPL.upl* files.
3. The procedure becomes iterative for a specified number of cycles. The scaffold restrains are extracted using a *regularize.cya*-based approach and saved in *_regula.upl* and *regula_new.aco* files.
4. The main structure calculation starts, using all previously extracted restrains. The resulting top 30 (default) structures ranked by cyana's target function are saved in a new bundle with incremental number.
5. *process.py* is called again on this new bundle to determine the best individual model.
---

# INPUT FILES

While most of the restrains are created automatically during the procedure, the scripts requires two kind of initial input files: sequence and PCS.

**Sequence**: the aa sequence of the protein has to be saved in a *SEQ.seq* file in which each line starts with the three letter code aa followed by white spaces and an incremental residue number starting at 1 (see example file). At the C-terminus of the protein extra pseudoatoms to represent the paramagnetic metal centers have to be added. A single pseudoatoms *PL* after the last protein residue links the protein to the following linker, made by a variable number of *LL5* pseudoatoms representing a 5 A flexible linker. They are then followed by an *ORI* pseudoatom representing a single tag (metal center). Add an *ORI* for each tag available and keep track of the corresponding residue number which will be used in the *process.py* (see section below).

**PCS**: PCS needs to be provided in two formats, .npc for *process.py* and .pcs for *CALC.cya*
	.npc: one file for each tag, named as *NPC_tag_x.npc* (replacing x with the tag number). One line for each shift with the following fields tab or white spaces separated: residue number, atom, shift (ppm), tolerance. (see example file)
	.pcs: one single file with all tags, named as *PCS.pcs* One line for each shift with the following fields tab or white spaces separated: residue number, aa (three letter code), atom, shift (ppm), tolerance, weight, sample (tag). (see example file).

---

# IN-DEPTH PROCESS.py UNDERSTANDING AND VARIABLES

This python 3 script essentially applies experimental PCS to a pdb bundle (either the initial model or the result of each refinement cycle) first to determine which model inside the bundle is the best (according to average Q-factor) and then to extract tensor restrains on the best structure. These restrains take the form of two tensor components (axial component and rhombicity) in a *_metal_centers.pcs* file and 6 upper and lower distance limits of the tag to 6 nearby CA atoms (to provide cyana with info about the metal position) in *_ORI_LOL.lol* and *_ORI_UPL.upl* files.

This is accomplished by four custom functions:
* `get_qfac(protein_file, npc, res_num, q_fac_dict)` (lines 27-51): using paramagpy, loads a `protein_file`, applies PCS in the form of the `npc` and perform a tensor search first through grid search and then through non-linear regression starting the search from the CA atom of `res_num`. Results are appended to a `q_fac_dict` for further processing.
* `manipulate_dict(dictionary)` (lines 52-72): looks at every key (model in the bundle) of the `dictionary`, create an average of its values (q-factors) and then extract the lowest averaged Q-factor. Appends all the results to a *Q-factors.txt* file for manual consultation and can be modified to list all individual Q-factors if needed. Returns the index of the best module in the bundle.
* `write_best_model(protein_file, model, suffix='_best.pdb')` (lines 114-121): simply extract the model with index `model` from the bundle `protein_file` and writes it as a pdb with suffix `suffix`. Used both for further processing and easy access to the best structure at every cycle for manual analysis.
* `get_tensor(protein_file, model, npc, res_num, tag_number)` (lines 73-113): again using paramagpy, loads a `protein_file`, applies PCS in the form of the `npc` and perform a tensor search first through grid search and then through non-linear regression starting the search from the CA atom of `res_num`. When the tensor is found, its axial component and rhombicity are saved in a file using cyana .pcs format. At last, the distance (A) between the found position and 6 nearby CA atoms is saved with a set tolerance into an upper and lower distance limit files.

## Variables

The reported example is set for a situation in which there are three tags in specific locations. Since your situations might vary, it's extremely important that the following variables are correct before starting a calculation:

* Number of tags: Lines 11, 13-15 state how any and which files contains the .npc. Modify line 11 `script, prot_file_name, npc_1_f, npc_2_f, npc_3_f = sys.argv` to reflect how many available tags you have. For example, if you have another tag add a `, npc_4_f` after `npc_3_f` or if you have only two, delete tag 3 and so on. For each available tag there should be a statement below that telling the system where to find this file `npc_x = os.getcwd() + "/" + npc_x_f` where x is the number of the tag.
Additionally, lines 129-131 and 138-140 also have to reflect the available number of tags. Each tag should have a `get_qfac` statement such as `get_qfac(prot_f, npc_x, res_dictionary[x], q_fact_dict)` and a `get_tensor` statement such as `get_tensor(prot_f, min_model, npc_x, res_dictionary[x], x)` where x is the tag number.

* Attachment site dictionary: line 18 states to which residue is attached each tag. It is used to start the tensor search. Modify the dictionary to reflect your situation and number of tags. For example if you have only two tags, the first attached to residue 200 and the first to residue 50, the code will be `res_dictionary = {1: 200, 2:50}`

* Tag dictionary: line 19 states which tag is which *ORI* residue in the *SEQ.seq* file. Modify the dictionary to reflect your situation and number of tags.

* Nearby residues dictionary: line 20 states which 6 residues for each tag should be used to determine *ORI* distances in *get_tensor()*. Modify the dictionary to reflect your situation and number of tags. For example if you have only two tags and you want distances for the first to residues 1-6 and for the second 200-206 the code will be `Ori_dictionary = {1: [1, 2, 3, 4, 5, 6], 2: [200, 201, 202, 203, 204, 205, 206]}`

---

# IN-DEPTH CALC.cya UNDERSTANDING AND VARIABLES

This cyana macro runs the entire process from start to end. It starts by copying the input model structure *model.pdb* as *1.pdb* for incremental numbering and calling *process.py* on it (see section above).
Next, it starts the main iteration. The first step is to extract "scaffold" restrains in the form of UPL and ACO according to a *regularize.cya* macro (see paper for more details). This is accomplished in lines 17-31.
After that, all previously generated restrains are read and the main structure calculation is performed in lines 36-57. Finally, the resulting bundle is saved and *process.py* called on it to setup the next cycle.

## Variables

You will have to check and verify the following variables before starting your own calculation:

* Number of cycles: line 1 states how many cycles the refinement will go through, `sintax cycles=@i=x` where x is the total number of cycles.

* Number of tags: lines 9 and 66 call *process.py* on a certain structure with available tags. You will have to modify the part regarding tags to reflect your own situation, as already stated above.

* Scaffold restrains generation: lines 17-31 handle the restrains extraction and can be modified to reflect your needs. In the reported example, each cap/module is restrained a part in line 21.

* Weights and anneal weights: lines 46-53 handle weights of individual type of restrains. In this example, the UPL weight is set to 1.0 (default) while the PCS weight is increased to 30 to reflect the best conditions determined in the paper. Your situation may vary. I suggest to experiment with few values for the UPL/PCS weight ration while not changing ACO and VdW. Also it is possible but not recommended to change the weight of individual restrains only in specific anneal phases. For example if you want to have a weight value of 30 for PCS only in phase 1, 2 and 4, while you want half the value in phase 3, the following two lines of code should do the work `weight_pcs=30.0` `anneal_weight_pcs := 1.0, 1.0, 0.5, 1.0`

* Seed: line 54 specify which seed should be used for reproducibility

* Number of processor: line 55 specify the number of processors (`nproc`) when using parallelization.

* Number of annealed structures and MD steps: line 56 states how many random structures (`structures`) should be annealed and in how many MD steps (`steps`).

* Resulting bundle size and RMSD: line 57 specify how many structures ranked by best target function will compose the final bundle (`structures`) and the residue range to be used for RMSD calculations (`range`)

---

# EXAMPLE

This folder contains all the input files already setup to run an example refinement. As explained in the paper, the example is set to refine a model structure of the dArmRP YM4A (*model.pdb*) using three tags in 2 cycles annealing 50 structures in 25000 MD using 16 processors, to provide a fast calculation (for a proper run I suggest at least 5 cycles with 500-1000 structures, depending on the PCS quality).
Be sure to have installed the required dependencies, adjust the available number of processors and clone the folder. In a terminal, navigate to the folder, starts cyana (i.e. `/applications/cyana_3.98.12/cyana`) and then simply call the main cyana macro with `CALC`.
For debugging reasons, the output is quite verbose and should be printed in the terminal and several other log files.

---

# TROUBLESHOOT

If you encounter any technical problem using this protocol, feel free to contact the authors: stefano.cucuzza@uzh.ch, oliver.zerbe@uzh.ch

---

# TL;DR

Go back and read the all thing. This is not reddit, sorry.