#!/usr/bin/env python3
"""
Convert FastANI results to distance matrix.

Usage:
    python ani_to_distance.py <ani_results> <distance_matrix> [--min-ani 0.8]
"""

import sys
import argparse
from pathlib import Path


def ani_to_distance(ani_results, distance_matrix, min_ani=0.8):
    """
    Convert FastANI results to distance matrix.
    Distance = 1 - (ANI/100)
    """
    # Read ANI results
    ani_data = {}
    samples = set()

    with open(ani_results, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 3:
                query = Path(parts[0]).stem
                reference = Path(parts[1]).stem
                ani = float(parts[2])

                samples.add(query)
                samples.add(reference)

                # Only keep if ANI meets threshold
                if ani >= min_ani * 100:
                    ani_data[(query, reference)] = ani

    # Convert to distance matrix
    samples = sorted(samples)

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
                    ani = ani_data.get((sample1, sample2), 0.0)
                    if ani == 0.0:
                        ani = ani_data.get((sample2, sample1), 0.0)

                    distance = 1.0 - (ani / 100.0) if ani > 0 else 1.0

                row.append(f"{distance:.6f}")

            f.write('\t'.join(row) + '\n')


def main():
    parser = argparse.ArgumentParser(description="Convert FastANI results to distance matrix")
    parser.add_argument("ani_results", help="FastANI results file")
    parser.add_argument("distance_matrix", help="Output distance matrix file")
    parser.add_argument("--min-ani", type=float, default=0.8, help="Minimum ANI threshold")

    args = parser.parse_args()

    ani_to_distance(args.ani_results, args.distance_matrix, args.min_ani)
    print(f"Distance matrix written to: {args.distance_matrix}")


if __name__ == "__main__":
    main()