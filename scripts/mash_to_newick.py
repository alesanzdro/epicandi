#!/usr/bin/env python3
"""
Convert Mash triangle distance matrix to Newick phylogenetic tree format.

Usage:
    python mash_to_newick.py <input_mash_distances> <output_newick>
"""

import sys
import numpy as np
from scipy.cluster.hierarchy import linkage, to_tree
from scipy.spatial.distance import squareform
import argparse


def parse_mash_triangle(file_path):
    """
    Parse the output of 'mash triangle' and return sample names and square distance matrix.

    Args:
        file_path (str): Path to mash triangle output file

    Returns:
        tuple: (sample_names, distance_matrix)
    """
    with open(file_path, 'r') as f:
        lines = [line.strip() for line in f if line.strip()]

    # Skip the first line if it's just the number of samples
    try:
        int(lines[0])
        proc_lines = lines[1:]
    except ValueError:
        proc_lines = lines

    names = [proc_lines[0].split('\t')[0]]
    distances = []

    for line in proc_lines[1:]:
        parts = line.split('\t')
        names.append(parts[0])
        row_distances = [float(d) for d in parts[1:]]
        distances.append(row_distances)

    # Build symmetric distance matrix
    num_items = len(names)
    dist_matrix = np.zeros((num_items, num_items))

    row_idx = 0
    for i in range(1, num_items):
        for j in range(i):
            dist_matrix[i, j] = distances[row_idx][j]
            dist_matrix[j, i] = distances[row_idx][j]
        row_idx += 1

    return names, dist_matrix


def scipy_tree_to_newick(node, leaf_names):
    """
    Convert a SciPy tree node to Newick format recursively.

    Args:
        node: SciPy tree node
        leaf_names (list): Names of leaf nodes

    Returns:
        str: Newick format string for this node
    """
    if node.is_leaf():
        return f"{leaf_names[node.id]}"
    else:
        left = scipy_tree_to_newick(node.get_left(), leaf_names)
        right = scipy_tree_to_newick(node.get_right(), leaf_names)

        left_len = node.dist - node.get_left().dist
        right_len = node.dist - node.get_right().dist

        return f"({left}:{left_len:.6f},{right}:{right_len:.6f})"


def main():
    parser = argparse.ArgumentParser(
        description="Convert Mash triangle distance matrix to Newick tree format"
    )
    parser.add_argument("input", help="Input Mash triangle distance file")
    parser.add_argument("output", help="Output Newick tree file")
    parser.add_argument("--method", default="average", choices=["average", "single", "complete", "ward"],
                        help="Linkage method for hierarchical clustering (default: average)")

    args = parser.parse_args()

    # Parse mash distances
    names, dist_matrix = parse_mash_triangle(args.input)

    # Clean sample names (remove path and extension)
    clean_names = []
    for name in names:
        # Remove directory path and common extensions
        clean_name = name.split('/')[-1]  # Remove path
        clean_name = clean_name.replace('.fasta', '').replace('.fa', '').replace('.fas', '')
        clean_names.append(clean_name)

    # Convert to condensed format expected by linkage()
    condensed_dist_matrix = squareform(dist_matrix)

    # Perform hierarchical clustering
    linkage_matrix = linkage(condensed_dist_matrix, method=args.method)

    # Convert to SciPy tree
    scipy_tree = to_tree(linkage_matrix, rd=False)

    # Convert to Newick format
    newick_string = scipy_tree_to_newick(scipy_tree, clean_names)

    # Write output
    with open(args.output, 'w') as f:
        f.write(f"{newick_string};\n")

    print(f"Newick tree written to {args.output}")
    print(f"Tree includes {len(clean_names)} samples: {', '.join(clean_names)}")


if __name__ == "__main__":
    main()