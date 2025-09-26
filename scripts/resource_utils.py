#!/usr/bin/env python3
"""
Resource utility functions for EpiCandi pipeline.
Standardizes resource allocation across all rules.
"""

from typing import Union, Dict, Any


def get_resource(config: Dict[str, Any], rule_name: str, resource_type: str,
                default_value: Union[int, str] = None) -> Union[int, str]:
    """
    Get resource allocation for a specific rule and resource type.

    Args:
        config: Configuration dictionary
        rule_name: Name of the rule (matches config.yaml resources section)
        resource_type: Type of resource ('threads', 'mem', 'walltime')
        default_value: Default value if not found in config

    Returns:
        Resource value or default
    """
    try:
        # Try to get specific rule resource
        rule_resources = config.get("resources", {}).get(rule_name, {})
        if rule_resources and resource_type in rule_resources:
            return rule_resources[resource_type]

        # Final fallback to provided default first
        if default_value is not None:
            return default_value

        # Fall back to default resources if no custom default provided
        default_resources = config.get("resources", {}).get("default", {})
        if default_resources and resource_type in default_resources:
            return default_resources[resource_type]

        # Sensible defaults by resource type
        defaults = {
            "threads": 1,
            "mem": 1024,  # 1 GB
            "walltime": 60  # 1 hour
        }

        return defaults.get(resource_type, 1)

    except Exception as e:
        print(f"Warning: Error getting resource {resource_type} for rule {rule_name}: {e}")
        return default_value if default_value is not None else 1


def get_threads(config: Dict[str, Any], rule_name: str, default: int = 1) -> int:
    """Get thread count for a rule."""
    return get_resource(config, rule_name, "threads", default)


def get_memory(config: Dict[str, Any], rule_name: str, default: int = 1024) -> int:
    """Get memory allocation (MB) for a rule."""
    return get_resource(config, rule_name, "mem", default)


def get_walltime(config: Dict[str, Any], rule_name: str, default: int = 60) -> int:
    """Get walltime (minutes) for a rule."""
    return get_resource(config, rule_name, "walltime", default)


def format_resource_block(config: Dict[str, Any], rule_name: str) -> str:
    """
    Generate a standardized resource block for a Snakemake rule.

    Args:
        config: Configuration dictionary
        rule_name: Name of the rule

    Returns:
        Formatted resource block string
    """
    threads = get_threads(config, rule_name)
    memory = get_memory(config, rule_name)
    walltime = get_walltime(config, rule_name)

    return f"""    threads: {threads}
    resources:
        mem_mb={memory},
        walltime={walltime}"""


if __name__ == "__main__":
    # Test the utility functions
    import yaml

    # Example config structure
    test_config = {
        "resources": {
            "default": {"threads": 1, "mem": 1024, "walltime": 60},
            "flye": {"threads": 16, "mem": 65536, "walltime": 1440},
            "medaka": {"threads": 8, "mem": 49152, "walltime": 480}
        }
    }

    print("Testing resource utilities:")
    print(f"Flye threads: {get_threads(test_config, 'flye')}")
    print(f"Medaka memory: {get_memory(test_config, 'medaka')}")
    print(f"Unknown rule threads: {get_threads(test_config, 'unknown_rule')}")

    print("\nResource block for flye:")
    print(format_resource_block(test_config, 'flye'))