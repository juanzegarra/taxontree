# Author: Juan Zegarra
# Date 10/11/26
#
# This is a script that parses a newick tree file and uses taxallnomy phylogenetic grouping for determining LCA (Lowest Common ancestor)
# The algorithm searches progressively the inner nodes of a tree and determines pairwise LCA of all leafs in the group

#######################################################################################

from Bio import Phylo
from Bio import SeqIO
import re
import sys
from pathlib import Path
from collections import defaultdict

help_str = '''Help: Correct use of the program is as follows
        python3 find_lca.py <tree.nwk> <Path-to-taxallnomy_db> <query_id>

    tree_file = Tree input needs to be input in the .nwk file format
    Path-to-taxallnomy_db = System path for taxallnomy database that needs to be downloaded previously
    query_id = String for the fasta ID of the query used for finding LCA
'''
# Open taxallnomy file and build two dictionaries for fast search, indexes are stored as ints
# child_dict maps every node to its children (organized in a list)
# parent_dict maps every node to its parent node, the last node has parent 0
def load_taxonomy_tree(taxallnomy_file):
    child_dict = defaultdict(list)
    parent_dict = defaultdict(float)

    if not taxallnomy_file.exists() or not taxallnomy_file.isfile():
        print("Failed locating taxallnomy file, check if file exists or if the path is correct.")
        return
    with open(taxallnomy_file, "r") as f:

        for line in f:
            taxa = line.strip().split("\t")
            if len(taxa) == 2:
                parent_dict[taxa[0]] = float(taxa[1])
                child_dict[taxa[1]].append(float(taxa[0)])
            elif len(taxa) == 1 :
               parent_dict[float(taxa[0])] = 0.0

    return child_dict, parent_dict

# Parses a tree structure to find the taxids that are inside the header of uniprot sequences (OX=11111111)
def parse_leaf_names(tree):
    leaf_taxids = []
    taxa_names = tree.get_terminals()
    pattern = r"^.*\bOX=([0-9]+)\b.*"
    for taxa in taxa_names:
        match = re.search(pattern, taxa.name)
        if match:
            leaf_taxids.append(float(match[1]))
        else:
            print(f"OX field with taxon id missing from entry {taxa.name}")
    return leaf_taxids

# Backtraces the taxonomy of each of the leafs to return its entire taxonomic line
def build_taxonomy_for_leafs(parent_dict, leaf_taxids):

    ancestry_dict = defaultdict(list)

    for taxid in leaf_taxids:
        taxonomy = []
        parent_node = parent_dict[int(taxid)]
        while parent_node != 0.0:
            parent_node = parent_dict[parent_node]
            taxonomy.append(parent_node)
        else:
            ancestry_dict[float(taxid)] = taxonomy

    return ancestry_dict
# For two given nodes, finds the LCA between them
def get_LCA_from_node(node1, node2, ancestry_dict=defaultdict(list)):
    A_ancestry = ancestry_dict[node1]
    B_ancestry = ancestry_dict[node2]
    for i in range(len(A_ancestry)):
        while A_ancestry[i] != B_ancestry[i]:
            i += 1
        else:


def find_tree_paralogs(tree, ancestry_dict)

def main():
    tree_file = Path(sys.argv[1])
    taxallnomy_file = Path(sys.argv[2])
    query = sys.argv[3]

    try:
        if ".nwk" in tree_file.name:
            tree = Phylo.read(tree_file.name, "nwk")
        elif not tree_file.is_file() or tree_file.exists():
            print("Tree file not found, check name or path")
            sys.exit()
        else:
            print("Non-supported tree format, please input nwk file")
            sys.exit()
    except Exception as e:
        print(f"Error parsing tree file with Biopython, please check tree file\n Error type {e}")


    # Extract taxids from tree leafs (leafs are named in the uniprot entry format and contain OX=taxid field)
    child_dict, parent_dict = load_taxonomy_tree(taxallnomy_file.name)
    leaf_taxids = parse_leaf_names(tree)
    ancestry_dict = build_taxonomy_for_leafs(parent_dict, leaf_taxids)

    internal_nodes = tree.get_non_terminals()

    for node in internal_nodes[:-1]:
