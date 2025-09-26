#!/usr/bin/env python3
"""
Convert Mash triangle output to square distance matrix.

Usage:
    python mash_to_matrix.py <mash_distances> <distance_matrix>
"""

import sys
import argparse
from pathlib import Path


def mash_to_matrix(mash_distances, distance_matrix):
    """
    Convert Mash triangle output to square distance matrix.
    """
    # Read Mash triangle output
    distances = {}
    samples = []

    with open(mash_distances, 'r') as f:
        # First line contains sample names
        first_line = next(f).strip()
        samples = [Path(sample).stem for sample in first_line.split('\t')]

        # Read distance matrix (lower triangle)
        for i, line in enumerate(f):
            values = line.strip().split('\t')
            sample1 = samples[i + 1]  # i+1 because we skip the diagonal

            for j, distance in enumerate(values):
                if distance:  # Skip empty cells
                    sample2 = samples[j]
                    distances[(sample1, sample2)] = float(distance)
                    distances[(sample2, sample1)] = float(distance)  # Symmetric

    # Write square matrix
    with open(distance_matrix, 'w') as f:
        # Write header
        f.write('\t' + '\t'.join(samples) + '\n')

        # Write matrix
        for sample1 in samples:
            row = [sample1]
            for sample2 in samples:
                if sample1 == sample2:
                    distance = 0.0
                else:
                    distance = distances.get((sample1, sample2), 1.0)

                row.append(f"{distance:.6f}")

            f.write('\t'.join(row) + '\n')


def main():
    parser = argparse.ArgumentParser(description="Convert Mash triangle output to square distance matrix")
    parser.add_argument("mash_distances", help="Mash triangle distances file")
    parser.add_argument("distance_matrix", help="Output square distance matrix file")

    args = parser.parse_args()

    mash_to_matrix(args.mash_distances, args.distance_matrix)
    print(f"Square distance matrix written to: {args.distance_matrix}")


if __name__ == "__main__":
    main()