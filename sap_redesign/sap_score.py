#!/usr/bin/env python
from __future__ import division

# This program accepts arguments like this:

#./remove_superfluous_trp.py pdb1.pdb pdb2.pdb pdb3.pdb
# or
#./remove_superfluous_trp.py -in:file:silent my.silent

import os
import sys
import math

from pyrosetta import *
from pyrosetta.rosetta import *

import numpy as np
from collections import defaultdict
import time
import argparse
import itertools
import subprocess
import time

sys.path.append("/home/bcov/sc/random/npose")
import voxel_array
import npose_util
import npose_util_pyrosetta as nup


init("-corrections:beta_nov16  -in:file:silent_struct_type binary -keep_input_scores false"
    " -holes:dalphaball /home/bcov/dev_rosetta/main/source/external/DAlpahBall/DAlphaBall.gcc")


parser = argparse.ArgumentParser()
parser.add_argument("-in:file:silent", type=str, default="")
parser.add_argument("pdbs", type=str, nargs="*")
parser.add_argument("-zero_adjust", type=float, default=0)


args = parser.parse_args(sys.argv[1:])

pdbs = args.pdbs
silent = args.__getattribute__("in:file:silent")




alpha = "ACDEFGHIKLMNPQRSTVWY"

seq = ""
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

def get_per_atom_sasa2(pose, probe_size=1.1):
    sasas = core.id.AtomID_Map_double_t()
    rsd_sasa = utility.vector1_double()
    core.scoring.calc_per_atom_sasa(pose, sasas, rsd_sasa, probe_size, False)
    return sasas


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
    "A": 0.116,
    "C": 0.18,
    "D": -0.472,
    "E": -0.457,
    "F": 0.5,
    "G": 0.001,
    "H": -0.335,
    "I": 0.443,
    "K": -0.217,
    "L": 0.443,
    "M": 0.238,
    "N": -0.264,
    "P": 0.211,
    "Q": -0.249,
    "R": -0.5,
    "S": -0.141,
    "T": -0.05,
    "V": 0.325,
    "W": 0.378,
    "Y": 0.38,
}



# script_dir = os.path.dirname(os.path.realpath(__file__))
# xml = script_dir + "/py_xml/remove_superfluous_nonpolar.xml"


# objs = protocols.rosetta_scripts.XmlObjects.create_from_file(xml)


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


def add_to_score_map(og_map, to_add, prefix, suffix=""):
    for name, score in list(to_add.items()):    # this iterator is broken. use list()
        og_map[prefix + name + suffix] = score

def move_chainA_far_away(pose):
    pose = pose.clone()
    sel = core.select.residue_selector.ChainSelector("A")
    subset = sel.apply(pose)

    x_unit = numeric.xyzVector_double_t(1, 0, 0)
    far_away = numeric.xyzVector_double_t(100, 0, 0)

    protocols.toolbox.pose_manipulation.rigid_body_move(x_unit, 0, far_away, pose, subset)

    return pose


def which_chain(pose, resi):
    for i in range(1, pose.num_chains()+1):
        if ( pose.conformation().chain_begin(i) <= resi and pose.conformation().chain_end(i) >= resi):
            return i
    assert(False)

def get_monomer_score(pose, scorefxn):
    pose = pose.split_by_chain()[1]
    return scorefxn(pose)


def get_filter_by_name(filtername):
    the_filter = objs.get_filter(filtername)

    # Get rid of stochastic filter
    if ( isinstance(the_filter, pyrosetta.rosetta.protocols.filters.StochasticFilter) ):
        the_filter = the_filter.subfilter()

    return the_filter

def add_filter_to_results(pose, filtername, out_score_map):
    filter = get_filter_by_name(filtername)
    print("protocols.rosetta_scripts.ParsedProtocol.REPORT: ============Begin report for " + filtername + "=======================")
    if (isinstance(filter, protocols.simple_filters.ShapeComplementarityFilter)):
        value = filter.compute(pose)
        out_score_map[filtername] = value.sc
        out_score_map[filtername+"_median_dist"] = value.distance
    else:
        value = filter.report_sm(pose)
        out_score_map[filtername] = value
    print("============End report for " + filtername + "=======================")

def score_with_these_filters(pose, filters, out_score_map):
    for filtername in filters:
        add_filter_to_results(pose, filtername, out_score_map)

def cmd(command, wait=True):
    # print ""
    # print command
    the_command = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    if (not wait):
        return
    the_stuff = the_command.communicate()
    return str(the_stuff[0]) + str(the_stuff[1])

def atid(resnum, atno):
    return core.id.AtomID( atno, resnum )

abego_man = core.sequence.ABEGOManager()
def get_abego(pose, seqpos):
    return abego_man.index2symbol(abego_man.torsion2index_level1( pose.phi(seqpos), pose.psi(seqpos), pose.omega(seqpos)))

def dump_region(pose, start, end, name):
    residues = utility.vector1_unsigned_long()
    for i in range(start, end ):
        residues.append(i)

    to_dump = core.pose.Pose()
    core.pose.pdbslice(to_dump, pose, residues)
    pdbinfo = core.pose.PDBInfo( to_dump )
    to_dump.pdb_info(pdbinfo)
    to_dump.dump_pdb(name)



def get_consensus(letters):
    counts = defaultdict(lambda : 0, {})
    for letter in letters:
        counts[letter] += 1

    maxx_letter = 0
    maxx = 0
    for key in counts:
        if ( counts[key] > maxx ):
            maxx = counts[key]
            maxx_letter = key
    return maxx_letter

# this is 1 indexed with the start and end with loops converted to nearby dssp
# and HHHHHH turns identified
def better_dssp2(pose, length=-1):
    if ( length < 0 ):
        length = pose.size()

    dssp = core.scoring.dssp.Dssp(pose)
    dssp.dssp_reduced()
    the_dssp = "x" + dssp.get_dssp_secstruct()
    the_dssp = list(the_dssp)

    n_consensus = get_consensus(the_dssp[3:6])

    the_dssp[1] = n_consensus
    the_dssp[2] = n_consensus
    the_dssp[3] = n_consensus
    the_dssp[4] = n_consensus
    the_dssp[5] = n_consensus

    c_consensus = get_consensus(the_dssp[-5:-2])

    the_dssp[-1] = c_consensus
    the_dssp[-2] = c_consensus
    the_dssp[-3] = c_consensus
    the_dssp[-4] = c_consensus
    the_dssp[-5] = c_consensus

    the_dssp = "".join(the_dssp)

    my_dssp = "x"

    for seqpos in range(1, length+1):
        abego = get_abego(pose, seqpos)
        this_dssp = the_dssp[seqpos]
        if ( the_dssp[seqpos] == "H" and abego != "A" ):
            # print("!!!!!!!!!! Dssp - abego mismatch: %i %s %s !!!!!!!!!!!!!!!"%(seqpos, the_dssp[seqpos], abego))

            # This is the Helix-turn-helix HHHH case. See the test_scaffs folder
            if ( abego == "B" or abego == "E" and seqpos > 5 and seqpos < len(the_dssp)-5 ):
                this_dssp = "L"

        my_dssp += this_dssp

    return my_dssp



# assumes dssp starts with X and removes it
def get_ss_elements2(dssp):
    assert(dssp[0] == "x")
    ss_elements = []

    offset = 0
    ilabel = -1
    for label, group in itertools.groupby(dssp):
        ilabel += 1
        this_len = sum(1 for _ in group)
        next_offset = offset + this_len

        ss_elements.append( (label, offset, next_offset-1))

        offset = next_offset
    return ss_elements[1:]


def delete_residues_smart(pose, start, end):

    ft = pose.fold_tree()
    if ( start == 1 ):
        ft.reorder(pose.size())
    else:
        ft.reorder(1)
    pose.fold_tree(ft)

    pose.delete_residue_range_slow(start, end)
    pose.conformation().detect_disulfides()

def ss_elem_subset(elem, pose):
    subset = utility.vector1_bool(pose.size())

    label, start, end = elem

    for i in range(start, end+1):
        subset[i] = True

    return subset


the_locals = None


# Developability index: a rapid in silico tool for the screening of antibody aggregation propensity.
def sap_score(pose, name_no_suffix, out_score_map, out_string_map, suffix):


    # R from the paper
    R = 5

    pose = pose.split_by_chain()[1]
    surf_vol = get_per_atom_sasa(pose)


    # get the per res base stats

    res_max_sasa = [None]
    res_hydrophobicity = [None]

    for resnum in range(1, pose.size()+1):
        letter = pose.residue(resnum).name1()
        res_max_sasa.append(max_sasa[letter])
        res_hydrophobicity.append(hydrophobicity[letter] + args.zero_adjust)


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
            # if ( atom_sasa[-1] > 1000 ):
            #     print("AARF", atom_sasa[-1], resnum, at, res.natoms(), len(atom_sasa), surf_vol.vol(resnum, at))
            # print("AAAA", atom_sasa[-1], resnum, at, res.natoms(), len(atom_sasa), surf_vol.vol(resnum, at))

    # print("TOTAL: ", surf_vol.tot_surf)

    # return None
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
            # print(clashgrid.arr[tuple(index)])
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
                print(ot_resnum, pose.residue(ot_resnum).name1(), res_sasa, res_max_sasa[ot_resnum], res_hydrophobicity[ot_resnum])

            atom_score += res_score

        pdb_info.bfactor(resnum, at, atom_score)
        sap_scores.append(atom_score)

    sap_scores = np.array(sap_scores)

    sap_score = np.sum( sap_scores[sap_scores > 0])
    print("sap score: %.1f"%sap_score)



    out_score_map['sap_score'] = sap_score




    return pose



############### BEGIN MAIN FUNCTION ###########################

if ( silent != "" ):
    sfd_in = rosetta.core.io.silent.SilentFileData(rosetta.core.io.silent.SilentFileOptions())
    sfd_in.read_file(silent)

    pdbs = list(sfd_in.tags())

    sfd_out = core.io.silent.SilentFileData( "out.silent", False, False, "binary", core.io.silent.SilentFileOptions())


num = -1
for pdb in pdbs:
    t0 = time.time()
    print("Attempting pose: " + pdb)

    # try:
    for k in [1]:
        if ( silent == "" ):
            pose = pose_from_file(pdb)
        else:
            pose = Pose()
            sfd_in.get_structure(pdb).fill_pose(pose)

        name_no_suffix = my_rstrip(my_rstrip(os.path.basename(pdb), ".gz"), ".pdb")

        sfd = core.io.raw_data.ScoreFileData("score.sc")

        score_map = std.map_std_string_double()
        string_map = std.map_std_string_std_string()


        out_pose = sap_score(pose, name_no_suffix, score_map, string_map, "")


        core.io.raw_data.ScoreMap.add_arbitrary_score_data_from_pose( pose, score_map)
        core.io.raw_data.ScoreMap.add_arbitrary_string_data_from_pose( pose, string_map)

        sfd.write_pose( pose, score_map, name_no_suffix, string_map)
        if (out_pose != None):


            # pdb_info = core.pose.PDBInfo(pose)
            # pose.pdb_info(pdb_info)
            if ( silent == "" ):
                out_pose.dump_pdb(name_no_suffix + ".pdb")
            else:
                struct = sfd_out.create_SilentStructOP()
                struct.fill_struct(out_pose, name_no_suffix)
                sfd_out.add_structure(struct)


        seconds = int(time.time() - t0)

        print("protocols.jd2.JobDistributor: " + name_no_suffix + " reported success in %i seconds"%seconds)

    # except Exception as e:
    #     print("Error!!!")
    #     print(e)



if ( silent != "" ):
    sfd_out.write_all("out.silent", False)




















