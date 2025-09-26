#!/usr/bin/env python3
"""
Script to update resource allocation patterns in Snakefile to use standardized functions.
"""

import re
import sys


def update_resource_patterns(file_path):
    """Update resource allocation patterns in a Snakefile."""

    with open(file_path, 'r') as f:
        content = f.read()

    # Pattern to match old resource allocation
    patterns = [
        # Pattern 1: Direct config access
        (
            r'threads: config\["resources"\]\["(\w+)"\]\["threads"\]',
            r'threads: get_resource(config, "\1", "threads")'
        ),
        (
            r'mem_mb=config\["resources"\]\["(\w+)"\]\["mem"\]',
            r'mem_mb=get_resource(config, "\1", "mem")'
        ),
        (
            r'walltime=config\["resources"\]\["(\w+)"\]\["walltime"\]',
            r'walltime=get_resource(config, "\1", "walltime")'
        ),

        # Pattern 2: Conditional access with fallbacks
        (
            r'threads: config\["resources"\]\["(\w+)"\]\["threads"\] if "(\w+)" in config\["resources"\] else (\d+)',
            r'threads: get_resource(config, "\1", "threads", \3)'
        ),
        (
            r'mem_mb=config\["resources"\]\["(\w+)"\]\["mem"\] if "(\w+)" in config\["resources"\] else (\d+)',
            r'mem_mb=get_resource(config, "\1", "mem", \3)'
        ),
        (
            r'walltime=config\["resources"\]\["(\w+)"\]\["walltime"\] if "(\w+)" in config\["resources"\] else "(\d+):(\d+):(\d+)"',
            r'walltime=get_resource(config, "\1", "walltime", \3)'
        ),
    ]

    original_content = content

    for pattern, replacement in patterns:
        content = re.sub(pattern, replacement, content)

    if content != original_content:
        with open(file_path, 'w') as f:
            f.write(content)
        print(f"Updated resource patterns in {file_path}")
        return True
    else:
        print(f"No changes needed in {file_path}")
        return False


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python update_resources.py <snakefile_path>")
        sys.exit(1)

    file_path = sys.argv[1]
    updated = update_resource_patterns(file_path)

    if updated:
        print("Resource patterns updated successfully!")
    else:
        print("No resource patterns needed updating.")