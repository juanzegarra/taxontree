from mimetypes import common_types
from pickletools import dis

from Bio import Phylo
from Bio import SeqIO
from Bio import AlignIO
import re
import sys
from Bio.Align import substitution_matrices
from Bio.Align import PairwiseAligner
from Bio.Phylo.TreeConstruction import DistanceCalculator

# 1. Initialize the aligner
aligner = PairwiseAligner()
aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")

entry1 = sys.argv[1]
entry2 = sys.argv[2]
format1 = entry1.split(".")[1]
format2 = entry2.split(".")[1]
distance_mode = sys.argv[3]
alignment_file = sys.argv[4]


format_dict = {"nwk" : "newick", "tree" : "newick", "nex" : "nexus", "dnd" : "newick"}

query_tree = Phylo.read(entry1, format_dict[format1])

subject_tree = Phylo.read(entry2, format_dict[format2])

query_leaves = set(query_tree.get_terminals())
subject_leaves = set(subject_tree.get_terminals())

nodes1 = list(query_tree.get_nonterminals())
nodes2 = list(subject_tree.get_nonterminals())

common_leaves = query_leaves.intersection(subject_leaves)

depth_node_query = query_tree.depths(unit_branch_lengths=True)
depth_node_subject = subject_tree.depths(unit_branch_lengths=True)

if distance_mode == "genetic":
    sequences = list(SeqIO.parse(alignment_file, "fasta"))
    calculator = DistanceCalculator("blosum62")

# Walk through most recent nodes, and find the range of each node by finding the most two most distant leaves within the branch
# Implement both genetic distance and branch length, maybe transform branch into string representation for blast alignment ()
def find_clades_ranges_query(query_leaves, subject_leaves, common_leaves, nodes, distance_mode):
    unique_leaves = query_leaves - subject_leaves
    clades = []
    comparisons_made = {}
    for node in reversed(nodes):
        # Dont add the root node to the clades
        if nodes.index(node) == 0:
            break
        if len(node.get_terminals()) <= 2:
            clades.append(tuple(node.get_terminals()))
        else:
            leaves = [leaf for leaf in node.get_terminals() if leaf in common_leaves]
            max_distance = (None, None, 0)
            if len(leaves) < 2:
                continue
            if distance_mode == 'branch':
                for i, leaf in enumerate(leaves):
                    for j, other_leaf in enumerate(leaves):
                        if i != j and (i, j) not in comparisons_made:
                            distance = query_tree.distance(leaf, other_leaf)
                            if distance > max_distance[2]:
                                max_distance = (leaf, other_leaf, distance)
                            comparisons_made[(leaf, other_leaf)] = distance
                else:
                    if max_distance[2] > 0:
                        clades.append(tuple(max_distance[:1]))
            # Genetic distance measures the exact distance between all possible pair of leaves inside the clade
            # The range of the clade is delimited by the maximum distance between any pair of leaves
            elif distance_mode == 'genetic':
                leaves_ids = [leaf.id for leaf in leaves]
                max_distance = (None, None, 0)
                aligned_seqs = [seq for seq in sequences if seq.id == leaf.id]
                reconstrured_msa = aligner.align(aligned_seqs[0], aligned_seqs[1])
                matrix = calculator.get_distance(reconstrured_msa)
                matrix_names = matrix.names
                for i in range(len(matrix_names)):
                    for j in range(i + 1, len(matrix_names)):
                        seq_a = matrix_names[i]
                        seq_b = matrix_names[j]
                        distance = matrix[i, j]
                        if distance > max_distance[2]:
                            max_distance = (seq_a, seq_b, distance)
                        comparisons_made[(seq_a, seq_b)] = distance
                else:
                    if max_distance[2] > 0:
                        clades.append(tuple(max_distance[:1]))
            elif distance_mode == 'blast':
                branch_string = "".join([leaf.id for leaf in leaves])
                clades.append(branch_string)
    return clades, comparisons_made


def fast_compare_trees(query, subject, distance_mode):
    query_clades, query_comparisons = find_clades_ranges_query(query_leaves, subject_leaves, common_leaves, nodes1, distance_mode)
    subject_clades, subject_comparisons = find_clades_ranges_query(subject_leaves, query_leaves, common_leaves, nodes2, distance_mode)

    if distance_mode == 'blast':
        aligner = PairwiseAligner()
        aligner.mode = "global"
        aligner.match_score = 2
        aligner.mismatch_score = -1
        aligner.open_gap_score = -2
        aligner.extend_gap_score = -1

        for query_clade in query_clades:
            max_score = 0
            for subject_clade in subject_clades:
                score = aligner.score(query_clade, subject_clade)
                if score > max_score:
                    max_score = score
