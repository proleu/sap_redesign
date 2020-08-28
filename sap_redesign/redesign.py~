#!/usr/bin/env python
# This program accepts arguments like this:
#./remove_superfluous_trp.py pdb1.pdb pdb2.pdb pdb3.pdb
# or
#./remove_superfluous_trp.py -in:file:silent my.silent
# TODO fix max SASA calculation to read in a file, instead of calculating it everytime this damn script is called
# TODO remove all unused utility functions
# TODO figure out why voxel_array sometimes crashes?
# TODO better documentation, delta SAP? Delta SAP is somewhat nontrivial as it requires you to capture stdout while calling a function.
# TODO more options: 
# design recursively on small numbers since the offending ones are usually far away from each other
# i also need to make it dump an individual score file 
# so that there isn't any scorefile corruption if it is run in batches

# python libraries
from __future__ import division
__author__ = "Brian Coventry, Philip Leung"
__copyright__ = None
__credits__ = ["Brian Coventry", "Philip Leung", "Rosettacommons"]
__license__ = "MIT"
__version__ = "0.4.0"
__maintainer__ = "Philip Leung"
__email__ = "pleung@cs.washington.edu"
__status__ = "Prototype"
import argparse
from collections import defaultdict
import itertools
import math
import os
import subprocess
import sys
import time
# external libraries
import numpy as np
from pyrosetta import *
from pyrosetta.rosetta import *
from pyrosetta.rosetta.core.io.silent import (SilentFileData, 
        SilentFileOptions)
from pyrosetta.rosetta.core.chemical import aa_from_oneletter_code
from pyrosetta.rosetta.core.pose import Pose, pose_residue_is_terminal
from pyrosetta.rosetta.protocols.simple_moves import (MutateResidue, 
                                                      SimpleThreadingMover)
from pyrosetta.rosetta.core.select import get_residues_from_subset
from pyrosetta.rosetta.core.select import residue_selector
from pyrosetta.rosetta.core.select.residue_selector import (
    AndResidueSelector, NeighborhoodResidueSelector, NotResidueSelector,
    OrResidueSelector, PrimarySequenceNeighborhoodSelector,
    ResidueIndexSelector)
from pyrosetta.rosetta.core.scoring import ScoreFunction, ScoreType
from pyrosetta.rosetta.core.scoring.methods import EnergyMethodOptions
from pyrosetta.rosetta.core.pack.task import operation
from pyrosetta.rosetta.protocols.task_operations import (LinkResidues,
    LimitAromaChi2Operation, PruneBuriedUnsatsOperation)
from pyrosetta.rosetta.protocols.aa_composition import (
    AddCompositionConstraintMover)
from pyrosetta.rosetta.protocols.denovo_design.movers import FastDesign
from pyrosetta.rosetta.protocols.protein_interface_design import (
    FavorNativeResidue)
from pyrosetta.rosetta.protocols.relax import FastRelax
from pyrosetta.rosetta.protocols.rosetta_scripts import XmlObjects
# Bcov libraries
# TODO remove this dependency if possible?
sys.path.append("/mnt//home/bcov/sc/random/npose")
import voxel_array
import npose_util
import npose_util_pyrosetta as nup
# TODO remove holes?
flags = """
-corrections::beta_nov16
-holes:dalphaball 
/home/bcov/dev_rosetta/main/source/external/DAlpahBall/DAlphaBall.gcc
-ignore_unrecognized_res 1
-in:file:silent_struct_type binary 
-keep_input_scores false
-mute core.select.residue_selector.SecondaryStructureSelector
-mute core.select.residue_selector.PrimarySequenceNeighborhoodSelector
-mute protocols.DsspMover
"""
pyrosetta.init(' '.join(flags.replace('\n\t', ' ').split()))
# TODO add more defense here
parser = argparse.ArgumentParser()
parser.add_argument("-in:file:silent", type=str, default='')
parser.add_argument("--pdbs", type=str, nargs='*')
parser.add_argument("-zero_adjust", type=float, default=0)
parser.add_argument("-worst_n", type=int, default=25)
parser.add_argument("-radius", type=int, default=5)
parser.add_argument("-flexbb", dest='flexbb', action='store_true')
parser.add_argument("-use_sc_neighbors", dest='use_sc_neighbors',
        action='store_true')
parser.add_argument("-lock_resis", type=int, nargs='*', default=[])
parser.add_argument("-cutoffs", type=float, nargs='*', default=[20,40])
parser.add_argument("-relax_script", type=str, default='MonomerDesign2019')
parser.add_argument("-up_ele", dest='up_ele', action='store_true')
parser.add_argument("-no_prescore", dest='prescore', action='store_false')
parser.add_argument("-no_rescore", dest='rescore', action='store_false')

args = parser.parse_args(sys.argv[1:])

pdbs = args.pdbs
silent = args.__getattribute__("in:file:silent")
worst_n = args.worst_n
zero_adjust = args.zero_adjust
radius = args.radius
flexbb = args.flexbb
use_sc_neighbors = args.use_sc_neighbors
lock_resis = args.lock_resis
cutoffs = tuple(args.cutoffs)
relax_script = args.relax_script
up_ele = args.up_ele
prescore = args.prescore
rescore = args.rescore
# TODO
print(args)
# TODO put this info into a file and just load the file, it should be faster
alpha = "ACDEFGHIKLMNPQRSTVWY"
seq = ''
for letter in alpha:
    seq += "AA%sAA"%letter

sasa_pose = pose_from_sequence(seq)
scorefxn = get_fa_scorefxn()
all_sub = core.select.residue_selector.TrueResidueSelector().apply(sasa_pose)
protocols.toolbox.pose_manipulation.repack_these_residues(all_sub, sasa_pose, scorefxn)

def get_per_atom_sasa(pose, probe_size=1.1):
    atoms = core.id.AtomID_Map_bool_t()
    atoms.resize(pose.size())
    for i in range(1, pose.size()+1):
        atoms.resize( i, pose.residue(i).natoms(), True)
    surf_vol = core.scoring.packing.get_surf_vol( pose, atoms, probe_size)
    # print(surf_vol.tot_surf)
    # print(surf_vol.surf(2, 1))  # this is per atom sasa (residue 2, atom 1)
    return surf_vol

surf_vol = get_per_atom_sasa(sasa_pose)

max_sasa = {}
for i in range(len(alpha)):
    resnum = i*5+3
    letter = alpha[i]

    sasa = 0
    res = sasa_pose.residue(resnum)
    assert(res.name1() == letter)
    for atno in range(1, res.natoms()+1):
        if ( res.atom_is_backbone(atno) ):
            continue
        sasa += surf_vol.surf(resnum, atno)

    max_sasa[letter] = sasa

# Development of hydrophobicity parameters to analyze proteins which bear post- or cotranslational modifications
# then you subtract 0.5 from scaled
hydrophobicity = {
    'A': 0.116,
    'C': 0.18,
    'D': -0.472,
    'E': -0.457,
    'F': 0.5,
    'G': 0.001,
    'H': -0.335,
    'I': 0.443,
    'K': -0.217,
    'L': 0.443,
    'M': 0.238,
    'N': -0.264,
    'P': 0.211,
    'Q': -0.249,
    'R': -0.5,
    'S': -0.141,
    'T': -0.05,
    'V': 0.325,
    'W': 0.378,
    'Y': 0.38,
}
# TODO remove this?
scorefxn = core.scoring.ScoreFunctionFactory.create_score_function("beta_nov16")

def fix_scorefxn(sfxn, allow_double_bb=False):
    opts = sfxn.energy_method_options()
    opts.hbond_options().decompose_bb_hb_into_pair_energies(True)
    opts.hbond_options().bb_donor_acceptor_check(not allow_double_bb)
    opts.hbond_options().use_hb_env_dep(True)
    #opts.elec_context_dependent(True)

    sfxn.set_energy_method_options(opts)

fix_scorefxn(scorefxn)

def my_rstrip(string, strip):
    if (string.endswith(strip)):
        return string[:-len(strip)]
    return string

# Developability index: a rapid in silico tool for the screening of antibody aggregation propensity.
def sap_score(pose, radius, name_no_suffix, out_score_map, out_string_map,
        suffix):
    # R from the paper
    R = radius
    pose = pose.split_by_chain()[1]
    surf_vol = get_per_atom_sasa(pose)
    # get the per res base stats
    res_max_sasa = [None]
    res_hydrophobicity = [None]

    for resnum in range(1, pose.size()+1):
        letter = pose.residue(resnum).name1()
        res_max_sasa.append(max_sasa[letter])
        res_hydrophobicity.append(hydrophobicity[letter] + zero_adjust)
    # make the things required to find 5A neighbors 
    idx_to_atom = []
    xyzs = []
    atom_sasa = []
    for resnum in range(1, pose.size()+1):
        res = pose.residue(resnum)
        for at in range(1, res.natoms()+1):
            if ( res.atom_is_backbone(at) ):
                continue
            xyzs.append(nup.from_vector(res.xyz(at)))
            idx_to_atom.append([resnum, at])
            atom_sasa.append(surf_vol.surf(resnum, at))
    
    atom_sasa = np.array(atom_sasa)
    idx_to_atom = np.array(idx_to_atom)
    xyzs = np.array(xyzs)

    resl = 1

    low = np.min(xyzs, axis=0) - R*2 - resl*2
    high = np.max(xyzs, axis=0) + R*2 + resl*2

    print("Making neighbor grid")
    clashgrid = voxel_array.VoxelArray(low, high, np.array([resl]*3), object)
    for idx, _ in enumerate(clashgrid.arr.flat):
        clashgrid.arr.flat[idx] = []

    for ixyz, xyz in enumerate(xyzs):
        indices = clashgrid.indices_within_x_of(R+resl, xyz)
        for index in indices:
            clashgrid.arr[tuple(index)].append(ixyz)

    atom_grid_indices = clashgrid.floats_to_indices(xyzs)

    sap_scores = []

    pdb_info = pose.pdb_info()

    for iatom in range(len(xyzs)):
        xyz = xyzs[iatom]
        resnum, at = idx_to_atom[iatom]
        grid_index = atom_grid_indices[iatom]

        grid_list = np.array(list(clashgrid.arr[tuple(grid_index)]))

        distances = np.linalg.norm( xyzs[grid_list] - xyz, axis=-1)

        idx_within_R = grid_list[distances <= R]

        atoms_within_R = idx_to_atom[idx_within_R]
        resnums = np.unique(atoms_within_R[:,0])

        atom_score = 0
        for ot_resnum in resnums:
            ats_idx = idx_within_R[atoms_within_R[:,0] == ot_resnum]
            res_sasa = np.sum(atom_sasa[ats_idx])

            res_score = res_sasa / res_max_sasa[ot_resnum] * res_hydrophobicity[ot_resnum]
            if ( res_score > 1000 ):
                print(ot_resnum, pose.residue(ot_resnum).name1(),
                        res_sasa, res_max_sasa[ot_resnum], 
                        res_hydrophobicity[ot_resnum])

            atom_score += res_score

        pdb_info.bfactor(resnum, at, atom_score)
        sap_scores.append(atom_score)

    sap_scores = np.array(sap_scores)

    sap_score = np.sum( sap_scores[sap_scores > 0])
    print("sap score: %.1f"%sap_score)
    out_score_map['sap_score'] = sap_score
    return pose

def sfxn_hard_maker(const_bb=True, up_ele=False) -> ScoreFunction:
    """Sets up Bcov's reweighted score function that penalizes buried
    unsatisfied polars more highly so that Rosetta doesn't make as many
    mistakes. Also unsets lk_ball because lk_ball is slow. 
    
    Args:
        const_bb (bool): Set this to False if you don't know where to expect
        PRO. Sets approximate_buried_unsat_penalty_assume_const_backbone.
        up_ele (bool): Increase the bonus for electrostatics, good for making
        Rosetta design salt bridges and hydrogen bonds.

    Returns:
        sfxn_hard (ScoreFunction): The modified score function.
    
    TODO: 
        as of 08172020 there is not much evidence that this reweighting is 
        necessary, but it doesn't seem to harm anything. 
        approximate_buried_unsat_penalty_hbond_energy_threshold may depend, 
        -0.5 is default but -0.2/-0.25 might work too. Might be good to check
        if it is really a good idea to unset lk_ball. Current best practices 
        are to set beta16_nostab.wts as the weights. Hopefully beta_nov20 
        will be much better. For sequence conservation, currently constraints
        (res_type_constraint) is set externally and this seems to work fine.
        as of 08262020 I am not sure if up_ele is considered best practices, 
        but since this implementation is primarily intended for resurfacing I
        think it should be okay to leave in for now.
    """
    sfxn_hard = pyrosetta.create_score_function("beta_nov16.wts")
    sfxn_hard.set_weight(ScoreType.aa_composition, 1.0)
    sfxn_hard.set_weight(ScoreType.approximate_buried_unsat_penalty, 5.0)
    emo = sfxn_hard.energy_method_options()
    emo.approximate_buried_unsat_penalty_burial_atomic_depth(4.0)
    emo.approximate_buried_unsat_penalty_hbond_energy_threshold(-0.2)	
    if const_bb:
        emo.approximate_buried_unsat_penalty_assume_const_backbone(1)
    else:
        emo.approximate_buried_unsat_penalty_assume_const_backbone(0)
    if up_ele:
        sfxn_hard.set_weight(ScoreType.fa_elec, 1.4)
        sfxn_hard.set_weight(ScoreType.hbond_sc, 2.0)
    else:
        pass
    sfxn_hard.set_energy_method_options(emo)
    sfxn_hard.set_weight(ScoreType.lk_ball, 0)
    sfxn_hard.set_weight(ScoreType.lk_ball_iso, 0)
    sfxn_hard.set_weight(ScoreType.lk_ball_bridge, 0)
    sfxn_hard.set_weight(ScoreType.lk_ball_bridge_uncpl, 0)
    return sfxn_hard

def generic_layer_dict_maker() -> dict:
    """Just a function that puts all of the standard layer definitions and 
    their corresponding allowed amino acids into a convenient dictionary.

    Args:
        None

    Returns:
        layer_dict (dict): The dict mapping standard layer definitions to 
        their allowed amino acids.
    """
    layer_dict = {"core AND helix_start": 'AFILVWYNQSTHP',
                  "core AND helix": 'AFILVWYNQHM',
                  "core AND loop": 'AFGILPVWYDENQSTHM',
                  "core AND sheet": 'FILVWYDENQSTH',
                  "boundary AND helix_start": 'ADEHIKLNPQRSTVY',
                  "boundary AND helix": 'ADEHIKLNQRSTVYM',
                  "boundary AND loop": 'ADEFGHIKLNPQRSTVY',
                  "boundary AND sheet": 'DEFHIKLNQRSTVY',
                  "surface AND helix_start": 'DEHKPQR',
                  "surface AND helix": 'EHKQR',
                  "surface AND loop": 'DEGHKNPQRST',
                  "surface AND sheet": 'EHKNQRST',
                  "helix_cap": 'DNSTP'}
    return layer_dict

def layer_design_maker(cutoffs:tuple, use_dssp:bool, use_sc_neighbors:bool
                      ) -> operation.DesignRestrictions:
    """Given options, returns a layer design task operation.

    Args:
        cutoffs (tuple): The tuple to set the number of sidechain neighbors or
        SASA accessibility. cutoffs[0] is for core, cuttoffs[1] is for surface
        and everything in between is considered boundary. 
        use_dssp (bool): Whether to use DSSP to determine secondary structure. 
        This is probably a good idea since sometimes poses don't have info on 
        secondary structure.
        use_sc_neighbors: Whether to use the number of sidechain neighbors to 
        determine burial. If false, will use SASA to determine burial. 
        Sidechain neighbors seems a little less noisy.
        
    Returns:
        layer_design (operation.DesignRestrictions): The task operation set by
        the options, ready to be passed to a task factory. 
    """
    # setup residue selectors: find surfaces, boundaries, cores, ss elements 
    surface = residue_selector.LayerSelector()
    surface.set_layers(0, 0, 1), surface.set_cutoffs(*cutoffs)
    surface.set_use_sc_neighbors(int(use_sc_neighbors))
    boundary = residue_selector.LayerSelector()
    boundary.set_layers(0, 1, 0), boundary.set_cutoffs(*cutoffs)
    boundary.set_use_sc_neighbors(int(use_sc_neighbors))
    core = residue_selector.LayerSelector()
    core.set_layers(1, 0, 0), core.set_cutoffs(*cutoffs)
    core.set_use_sc_neighbors(int(use_sc_neighbors))
    sheet = residue_selector.SecondaryStructureSelector("E")
    sheet.set_overlap(0)
    sheet.set_minH(3), sheet.set_minE(3)
    sheet.set_include_terminal_loops(0)
    sheet.set_use_dssp(int(use_dssp))
    entire_loop = residue_selector.SecondaryStructureSelector("L")
    entire_loop.set_overlap(0)
    entire_loop.set_minH(3), entire_loop.set_minE(3)
    entire_loop.set_include_terminal_loops(1)
    entire_loop.set_use_dssp(int(use_dssp))
    entire_helix = residue_selector.SecondaryStructureSelector("H")
    entire_helix.set_overlap(0)
    entire_helix.set_minH(3), entire_helix.set_minE(3)
    entire_helix.set_include_terminal_loops(0)
    entire_helix.set_use_dssp(int(use_dssp))
    lower_helix = PrimarySequenceNeighborhoodSelector(1, 0, entire_helix)
    helix_cap = AndResidueSelector(entire_loop, lower_helix)
    upper_helix_cap = PrimarySequenceNeighborhoodSelector(0, 1, helix_cap)
    helix_start = AndResidueSelector(entire_helix, upper_helix_cap)
    not_helix_start = NotResidueSelector(helix_start)
    helix = AndResidueSelector(entire_helix, not_helix_start)
    not_helix_cap = NotResidueSelector(helix_cap)
    loop = AndResidueSelector(entire_loop, not_helix_cap)
    # setup layer design
    layer_dict = generic_layer_dict_maker()
    layer_design = operation.DesignRestrictions()
    for selector_logic, aas in layer_dict.items():
        rlto = operation.RestrictAbsentCanonicalAASRLT()
        rlto.aas_to_keep(aas)
        selector_objs = []
        if 'AND' in selector_logic:
            selectors = selector_logic.split(" AND ")
            for selector_str in selectors:
                selector_objs.append(eval(selector_str))
            selector = AndResidueSelector(*tuple(selector_objs))
        elif 'OR' in selector_logic:
            selectors = selector_logic.split(" OR ")
            for selector_str in selectors:
                selector_objs.append(eval(selector_str))
            selector = OrResidueSelector(*tuple(selector_objs))
        else:
            selector = eval(selector_logic)
        layer_design.add_selector_rlto_pair(selector, rlto)
    return layer_design

def design_pack_lock_maker(design_resis:list) -> tuple:
    """Given options, returns a design, pack, and lock task operations.

    Args:
        design_resis (list): A list of residues to design. Neighbors will be 
        allowed to repack, everything else will be locked. Cysteine is not 
        allowed as a residue for design, everything else is, so should be used
        in combination with other task operations. 
        
    Returns:
        design, pack, lock (tuple): The task operations set by the options,
        ready to be passed to a task factory. 
    """
    designable = ResidueIndexSelector()
    designable.set_index(','.join([str(x) for x in design_resis]))
    not_designable = NotResidueSelector()
    not_designable.set_residue_selector(designable)
    packable = NeighborhoodResidueSelector()
    packable.set_focus_selector(designable)
    not_packable = NotResidueSelector()
    not_packable.set_residue_selector(packable)
    no_cys = operation.RestrictAbsentCanonicalAASRLT()
    no_cys.aas_to_keep('ADEFGHIKLMNPQRSTVWY')
    no_design = operation.RestrictToRepackingRLT()
    no_repack = operation.PreventRepackingRLT()
    design = operation.OperateOnResidueSubset(no_cys, designable)
    pack = operation.OperateOnResidueSubset(no_design, not_designable)
    lock = operation.OperateOnResidueSubset(no_repack, not_packable)
    return design, pack, lock

def favor_native_residue_maker(sfxn: ScoreFunction, restraint: float
                              ) -> FavorNativeResidue:
    """Given options, returns a design, pack, and lock task operations.

    Args:
        sfxn (ScoreFunction): A Rosetta ScoreFunction. It will have the weight
        for the res_type_constraint set to 1.
        restraint (float): What bonus to pass to the FavorNativeResidue mover. 
        
    Returns:
        favor_native_residue (FavorNativeResidue): The the FavorNativeResidue
        mover set by the options, ready to be applied to a pose. 
    """
    sfxn.set_weight(ScoreType.res_type_constraint, 1.0)
    xml_string = """
    <MOVERS>
        <FavorNativeResidue name="favor_native_residue" bonus="{}"/>
    </MOVERS>
    """.format(restraint)
    xml_obj = XmlObjects.create_from_string(xml_string)
    favor_native_residue = xml_obj.get_mover('favournative')
    return favor_native_residue

def relax_script_maker(relax_script:str
                      ) -> pyrosetta.rosetta.std.vector_std_string():
    """Given an absolute or local path or a database relax script name, sets
    up a relax script for rosetta to read in after reading it in line by line.

    Args:
        relax_script (str): Somewhat flexibly implemented and can be an 
        absolute or local path or a database relax script
        
    Returns:
        script (pyrosetta.rosetta.std.vector_std_string): The relax script,
        ready to read into a mover. 
    """
    script = pyrosetta.rosetta.std.vector_std_string()
    path = "/software/rosetta/latest/database/sampling/relax_scripts/"
    # assumes database script if only the base name of a script is given
    if '/' not in relax_script and ".txt" not in relax_script:
        absolute_path = path + relax_script + ".txt"
    # assumes database script if the name of a script is given without a path
    elif '/' not in relax_script:
        absolute_path = path + relax_script
    # if there is a full name and path assumes a custom script
    else: 
        absolute_path = relax_script
    with open(absolute_path) as f:
        lines = f.readlines()
        for line in lines:
            script.append(' '.join(line.split()))
    return script

# TODO documentation
def fast_design_with_options(pose:Pose, to_design=[], cutoffs=(20,40), 
        flexbb=True, relax_script="MonomerDesign2019", restraint=0,
        up_ele=False, use_dssp=True, use_sc_neighbors=False) -> Pose:
    """"""
    sfxn_hard = sfxn_hard_maker(up_ele=up_ele)
    # determine which residues are designable
    true_sel = residue_selector.TrueResidueSelector()
    if len(to_design) == 0:
        design_resis = list(get_residues_from_subset(true_sel.apply(pose)))
    else:
        design_resis = to_design.copy()
    # setup task operations
    task_factory = pyrosetta.rosetta.core.pack.task.TaskFactory()
    design, pack, lock = design_pack_lock_maker(design_resis)
    layer_design = layer_design_maker(cutoffs, use_dssp, use_sc_neighbors)
    ic = operation.IncludeCurrent()
    arochi = LimitAromaChi2Operation()
    arochi.include_trp(True)
    ex1_ex2 = operation.ExtraRotamersGeneric()
    ex1_ex2.ex1(True), ex1_ex2.ex2(True)
    prune = PruneBuriedUnsatsOperation()
    for op in [design, pack, lock, layer_design, ic, arochi, ex1_ex2, prune]:
        task_factory.push_back(op)
    # setup movemap
    mm = pyrosetta.rosetta.core.kinematics.MoveMap()
    mm.set_bb(flexbb), mm.set_chi(True), mm.set_jump(False)
    # optionally enable FNR
    if restraint != 0:
        favor_native_residue = favor_native_residue_maker(sfxn=sfxn_hard,
                                                          restraint=restraint)
        favor_native_residue.apply(pose)
    else:
        pass
    # setup fast design
    fast_design = FastDesign(scorefxn_in=sfxn_hard, standard_repeats=1)
    script = relax_script_maker(relax_script)
    fast_design.set_script_from_lines(script)
    fast_design.cartesian(False)
    fast_design.set_task_factory(task_factory)
    fast_design.set_movemap(mm)
    fast_design.constrain_relax_to_start_coords(True)
    fast_design.minimize_bond_angles(False)
    fast_design.minimize_bond_lengths(False)
    fast_design.min_type("lbfgs_armijo_nonmonotone")
    fast_design.ramp_down_constraints(False)
    fast_design.apply(pose)
    return pose

def residue_sap_list_maker(pose:Pose) -> list:
    residue_sap_list = []
    for resi in range(1, pose.size()+1):
        residue = pose.residue(resi)
        residue_sap = 0
        for atom in range(1, residue.natoms()+1):
            if residue.atom_is_backbone(atom):
                continue
            else:
                residue_sap += pose.pdb_info().bfactor(resi, atom)
        residue_sap_list.append((resi, residue_sap))
        
    return residue_sap_list

# TODO implement actual main?
############### BEGIN MAIN FUNCTION ###########################


if silent != '':
    sfd_in = SilentFileData(SilentFileOptions())
    sfd_in.read_file(silent)
    pdbs = list(sfd_in.tags())
    sfd_out = SilentFileData("out.silent", False, False, "binary",
            SilentFileOptions())
# TODO is num used?
num = -1
for pdb in pdbs:
    t0 = time.time()
    print("Attempting pose: " + pdb)

    for k in [1]:
        if ( silent == '' ):
            pose = pose_from_file(pdb)
        else:
            pose = Pose()
            sfd_in.get_structure(pdb).fill_pose(pose)

        name_no_suffix = my_rstrip(my_rstrip(os.path.basename(pdb), ".gz"),
                ".pdb")
        sfd = core.io.raw_data.ScoreFileData("score.sc")
        score_map = std.map_std_string_double()
        string_map = std.map_std_string_std_string()
        if prescore:
            # get SAP score for the pose
            print("prescoring SAP:")
            pre_pose = sap_score(pose, radius, name_no_suffix, score_map,
                    string_map, '')
            core.io.raw_data.ScoreMap.add_arbitrary_score_data_from_pose(pose,
                    score_map)
            core.io.raw_data.ScoreMap.add_arbitrary_string_data_from_pose(
                    pose, string_map)
        else:
            # if prescore is set to false, assumes pose already has SAP info
            pre_pose = pose.clone()
        # use per residue SAP to make a list of the worst offenders
        residue_sap_list = residue_sap_list_maker(pre_pose)
        sorted_residue_sap_list = sorted(residue_sap_list, key=lambda x: x[1],
                                         reverse=True)
        # check to see if each worst resi is allowed to be designed
        if len(lock_resis) == 0:
            worst_resis = [x[0] for x in sorted_residue_sap_list[:worst_n]]
        else:
            worst_resis = []
            print("The residues that will not be designed:",
                    ' '.join(str(x) for x in lock_resis))
            for residue, sap in sorted_residue_sap_list:
                if len(worst_resis) >= worst_n:
                    break
                elif residue in lock_resis:
                    pass
                else:
                    worst_resis.append(residue)
        print("Worst residues by SAP that are allowed to be designed:",
                ' '.join(str(x) for x in worst_resis))
        # redesign a new pose targeting the worst residues
        new_pose = fast_design_with_options(pre_pose, to_design=worst_resis,
                cutoffs=cutoffs, flexbb=flexbb, relax_script=relax_script,
                restraint=0, up_ele=up_ele, use_dssp=True,
                use_sc_neighbors=use_sc_neighbors)
        if rescore:
            # rescore the designed pose
            print("rescoring SAP:")
            name_no_suffix += '_resurf'
            post_pose = sap_score(new_pose, radius, name_no_suffix, score_map,
                    string_map, '')
            sfd.write_pose(post_pose, score_map, name_no_suffix, string_map)
        else:
            post_pose = new_pose.clone()
        if (pre_pose != None):
            if ( silent == '' ):
                post_pose.dump_pdb(name_no_suffix + ".pdb")
            else:
                struct = sfd_out.create_SilentStructOP()
                struct.fill_struct(post_pose, name_no_suffix)
                sfd_out.add_structure(struct)

        seconds = int(time.time() - t0)
        print("protocols.jd2.JobDistributor: {0} reported success in {1} \
                seconds".format(name_no_suffix, seconds))

if ( silent != '' ):
    sfd_out.write_all("out.silent", False)

