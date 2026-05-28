#!/usr/bin/env python3
"""
upgma_network.py · Red UPGMA epidemiológica desde matriz snp-dists (v2)

v2: muestra el tier (short_read/long_read) en el título y leyenda.
    Si tier == long_read (Nanopore-only), añade aviso de menor resolución.
"""

import argparse
import sys
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import squareform


def load_matrix(matrix_path):
    # snp-dists emits TSV by default and CSV when invoked with -c. Auto-detect
    # by counting separators on the header line so callers do not need to align.
    with open(matrix_path) as fh:
        first = fh.readline()
    sep = ',' if first.count(',') > first.count('\t') else '\t'
    return pd.read_csv(matrix_path, sep=sep, index_col=0)


def load_samplesheet(samplesheet_path):
    df = pd.read_csv(samplesheet_path)
    id_col = 'id' if 'id' in df.columns else 'sample_id'
    return df.set_index(id_col)


def build_upgma_tree(matrix_df):
    condensed = squareform(matrix_df.values, checks=False)
    Z = linkage(condensed, method='average')
    return Z, matrix_df.index.tolist()


def tree_to_network(Z, labels):
    G = nx.Graph()
    n = len(labels)
    for label in labels:
        G.add_node(label, kind='leaf')
    for i, (a, b, dist, _) in enumerate(Z):
        internal_id = f"node_{n + i}"
        G.add_node(internal_id, kind='internal')
        a_label = labels[int(a)] if int(a) < n else f"node_{int(a)}"
        b_label = labels[int(b)] if int(b) < n else f"node_{int(b)}"
        G.add_edge(a_label, internal_id, weight=dist / 2)
        G.add_edge(b_label, internal_id, weight=dist / 2)
    return G


def annotate_metadata(G, samples_df):
    for node in G.nodes():
        if G.nodes[node].get('kind') == 'leaf' and node in samples_df.index:
            for col in ['batch', 'clade', 'reference_tag']:
                if col in samples_df.columns:
                    G.nodes[node][col] = str(samples_df.loc[node, col])


def color_by_batch(G):
    batches = sorted({
        G.nodes[n].get('batch', 'na')
        for n in G.nodes() if G.nodes[n].get('kind') == 'leaf'
    })
    cmap = plt.cm.tab20
    palette = {b: cmap(i % 20) for i, b in enumerate(batches)}
    colors = []
    for n in G.nodes():
        if G.nodes[n].get('kind') == 'leaf':
            colors.append(palette[G.nodes[n].get('batch', 'na')])
        else:
            colors.append((0.7, 0.7, 0.7, 0.5))
    return colors, palette


def draw_network(G, out_svg, cohort_id, threshold, tier):
    colors, palette = color_by_batch(G)
    sizes = [350 if G.nodes[n].get('kind') == 'leaf' else 80 for n in G.nodes()]
    pos = nx.kamada_kawai_layout(G, weight='weight')

    fig, ax = plt.subplots(figsize=(14, 12))
    nx.draw_networkx_edges(G, pos, alpha=0.6, ax=ax)
    nx.draw_networkx_nodes(G, pos, node_color=colors, node_size=sizes, ax=ax)
    leaf_labels = {n: n for n in G.nodes() if G.nodes[n].get('kind') == 'leaf'}
    nx.draw_networkx_labels(G, pos, labels=leaf_labels, font_size=8, ax=ax)
    edge_labels = {
        (u, v): f"{int(round(d['weight'] * 2))}"
        for u, v, d in G.edges(data=True)
        if G.nodes[u].get('kind') == 'leaf' or G.nodes[v].get('kind') == 'leaf'
    }
    nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels, font_size=6, ax=ax)

    handles = [
        plt.Line2D([0], [0], marker='o', linestyle='', color=col,
                   label=str(b), markersize=10)
        for b, col in palette.items()
    ]
    ax.legend(handles=handles, title='Batch', loc='best', fontsize=8)

    title = f"{cohort_id} · UPGMA · transmission threshold = {threshold} SNPs · tier={tier}"
    if tier == 'long_read':
        title += "\n⚠ Nanopore-only: esperar ±2-3 SNPs de ruido vs duplicados técnicos"
    ax.set_title(title, fontsize=11)
    ax.axis('off')
    plt.tight_layout()
    plt.savefig(out_svg, format='svg', bbox_inches='tight')
    plt.close()


def find_transmission_pairs(matrix_df, samples_df, threshold, tier):
    pairs = []
    labels = matrix_df.index.tolist()
    for i in range(len(labels)):
        for j in range(i + 1, len(labels)):
            dist = matrix_df.iloc[i, j]
            if dist <= threshold:
                row = {
                    's1': labels[i],
                    's2': labels[j],
                    'snp_distance': int(dist),
                    'tier': tier,
                }
                for col in ['batch', 'clade']:
                    if col in samples_df.columns:
                        row[f's1_{col}'] = samples_df.loc[labels[i], col] if labels[i] in samples_df.index else 'na'
                        row[f's2_{col}'] = samples_df.loc[labels[j], col] if labels[j] in samples_df.index else 'na'
                pairs.append(row)
    return pd.DataFrame(pairs).sort_values('snp_distance') if pairs else pd.DataFrame()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--matrix',      required=True)
    ap.add_argument('--molten',      required=True)
    ap.add_argument('--samplesheet', required=True)
    ap.add_argument('--cohort',      required=True)
    ap.add_argument('--tier',        required=True, choices=['short_read', 'long_read'])
    ap.add_argument('--threshold',   type=int, default=6)
    ap.add_argument('--out-svg',     required=True)
    ap.add_argument('--out-gexf',    required=True)
    ap.add_argument('--out-pairs',   required=True)
    args = ap.parse_args()

    matrix_df = load_matrix(args.matrix)
    samples_df = load_samplesheet(args.samplesheet)

    Z, labels = build_upgma_tree(matrix_df)
    G = tree_to_network(Z, labels)
    annotate_metadata(G, samples_df)

    draw_network(G, args.out_svg, args.cohort, args.threshold, args.tier)
    nx.write_gexf(G, args.out_gexf)

    pairs_df = find_transmission_pairs(matrix_df, samples_df, args.threshold, args.tier)
    pairs_df.to_csv(args.out_pairs, sep='\t', index=False)
    print(f"[upgma_network] {args.cohort} (tier={args.tier}): {len(pairs_df)} pares ≤ {args.threshold} SNPs",
          file=sys.stderr)


if __name__ == '__main__':
    main()
