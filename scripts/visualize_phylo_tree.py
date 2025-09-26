#!/usr/bin/env python3
"""
Visualize phylogenetic trees using ETE3 with professional styling.

Usage:
    python visualize_phylo_tree.py <input_tree> <output_prefix> [--tree-type {triage,phylogeny}]
"""

import sys
import os
import argparse
from pathlib import Path
from ete3 import Tree, TreeStyle, NodeStyle, TextFace, CircleFace


def setup_tree_style(tree_type="triage"):
    """
    Configure tree visualization style.

    Args:
        tree_type (str): Type of tree ("triage" or "phylogeny")

    Returns:
        TreeStyle: Configured tree style object
    """
    ts = TreeStyle()
    ts.show_leaf_name = True
    ts.show_branch_length = True
    ts.show_branch_support = False
    ts.scale = 100 if tree_type == "triage" else 150
    ts.branch_vertical_margin = 10
    ts.margin_top = 20
    ts.margin_bottom = 20
    ts.margin_left = 20
    ts.margin_right = 50

    # Add title
    if tree_type == "triage":
        title = "Triage Tree - Assembly-to-Assembly Distances"
    else:
        title = "Phylogenetic Tree - Orthologous Gene Analysis"

    ts.title.add_face(TextFace(title, fsize=14, bold=True), column=0)

    return ts


def get_strain_info():
    """
    Define strain information for coloring and annotation.

    Returns:
        dict: Strain information with colors and metadata
    """
    return {
        "EPI003361": {"color": "#E74C3C", "clade": "Unknown", "type": "Sample"},
        "EPI003102": {"color": "#3498DB", "clade": "Unknown", "type": "Sample"},
        "EPI003135": {"color": "#2ECC71", "clade": "Unknown", "type": "Sample"},
        "B8441_ref": {"color": "#F39C12", "clade": "Clade I", "type": "Reference"},
        "Chaemulonii_out": {"color": "#95A5A6", "clade": "Outgroup", "type": "Outgroup"},
        # Add more strains as needed
    }


def style_tree_nodes(tree, strain_info, tree_type="triage"):
    """
    Apply styling to tree nodes based on strain information.

    Args:
        tree (Tree): ETE3 tree object
        strain_info (dict): Strain information dictionary
        tree_type (str): Type of tree for styling adjustments
    """
    for node in tree.traverse():
        if node.is_leaf():
            # Extract strain name (remove common suffixes)
            strain = node.name.split('.')[0]
            strain = strain.replace('_clean', '').replace('_masked', '')

            if strain in strain_info:
                ns = NodeStyle()
                ns["bgcolor"] = strain_info[strain]["color"]
                ns["size"] = 15 if tree_type == "triage" else 12
                ns["fgcolor"] = "white"
                node.set_style(ns)

                # Add strain type annotation
                strain_type = strain_info[strain]["type"]
                clade = strain_info[strain]["clade"]

                if strain_type != "Sample":
                    face = TextFace(f" ({strain_type})", fsize=8, fgcolor="gray")
                    node.add_face(face, column=1, position="branch-right")

                if clade != "Unknown":
                    face = TextFace(f" [{clade}]", fsize=8, fgcolor="blue")
                    node.add_face(face, column=2, position="branch-right")

        else:
            # Internal nodes
            ns = NodeStyle()
            ns["size"] = 3
            ns["fgcolor"] = "black"
            node.set_style(ns)

            # Show bootstrap support if available
            if hasattr(node, 'support') and node.support > 0:
                if node.support >= 0.7:  # High support
                    support_face = TextFace(f"{node.support:.2f}", fsize=6, fgcolor="green")
                else:  # Low support
                    support_face = TextFace(f"{node.support:.2f}", fsize=6, fgcolor="red")
                node.add_face(support_face, column=0, position="branch-top")


def root_tree(tree):
    """
    Root the tree using outgroup if available.

    Args:
        tree (Tree): ETE3 tree object to root

    Returns:
        Tree: Rooted tree
    """
    # Try to root at outgroup
    outgroup_names = ["Chaemulonii", "haemulonii", "outgroup", "out"]

    for leaf in tree.get_leaves():
        for outgroup_pattern in outgroup_names:
            if outgroup_pattern.lower() in leaf.name.lower():
                try:
                    tree.set_outgroup(leaf)
                    print(f"Tree rooted at outgroup: {leaf.name}")
                    return tree
                except Exception as e:
                    print(f"Warning: Could not root at {leaf.name}: {e}")

    print("No outgroup found. Tree remains unrooted.")
    return tree


def main():
    parser = argparse.ArgumentParser(
        description="Visualize phylogenetic trees with professional styling"
    )
    parser.add_argument("input_tree", help="Input tree file (Newick format)")
    parser.add_argument("output_prefix", help="Output file prefix")
    parser.add_argument("--tree-type", choices=["triage", "phylogeny"], default="triage",
                        help="Type of tree for appropriate styling")
    parser.add_argument("--width", type=int, default=1200, help="Image width in pixels")
    parser.add_argument("--height", type=int, default=800, help="Image height in pixels")
    parser.add_argument("--dpi", type=int, default=300, help="Image resolution (DPI)")

    args = parser.parse_args()

    # Check input file exists
    if not Path(args.input_tree).exists():
        print(f"Error: Input tree file '{args.input_tree}' not found.")
        sys.exit(1)

    # Set QT backend for headless rendering
    os.environ["QT_QPA_PLATFORM"] = "offscreen"

    try:
        # Load tree
        tree = Tree(args.input_tree, format=1)
        print(f"Loaded tree with {len(tree.get_leaves())} leaves")

        # Root tree if possible
        tree = root_tree(tree)

        # Get strain information and style tree
        strain_info = get_strain_info()
        style_tree_nodes(tree, strain_info, args.tree_type)

        # Configure tree style
        ts = setup_tree_style(args.tree_type)

        # Render tree in multiple formats
        output_png = f"{args.output_prefix}.png"
        output_svg = f"{args.output_prefix}.svg"
        output_pdf = f"{args.output_prefix}.pdf"

        tree.render(output_png, w=args.width, h=args.height, tree_style=ts, dpi=args.dpi)
        tree.render(output_svg, tree_style=ts)
        tree.render(output_pdf, tree_style=ts)

        print(f"Tree visualizations saved:")
        print(f"  PNG: {output_png}")
        print(f"  SVG: {output_svg}")
        print(f"  PDF: {output_pdf}")

        # Print tree statistics
        print(f"\n=== TREE STATISTICS ===")
        print(f"Number of leaves: {len(tree.get_leaves())}")
        print(f"Samples: {[leaf.name for leaf in tree.get_leaves()]}")

        # Print ASCII tree for quick inspection
        print(f"\n=== TREE TOPOLOGY (ASCII) ===")
        print(tree.get_ascii(show_internal=True, compact=True))

    except Exception as e:
        print(f"Error processing tree: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()