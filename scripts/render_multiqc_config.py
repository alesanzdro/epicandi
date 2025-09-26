#!/usr/bin/env python3
"""Renderiza el template de multiqc_config reemplazando placeholders con valores del config.yaml.
Si alguna clave falta, se usa "UNKNOWN".
"""
import yaml
from pathlib import Path
import sys

BASE = Path(__file__).resolve().parent.parent
CONFIG_YAML = BASE / 'config.yaml'
TEMPLATE = BASE / 'conf' / 'multiqc' / 'multiqc_config.template.yaml'
OUTPUT = BASE / 'conf' / 'multiqc' / 'multiqc_config.rendered.yaml'

PLACEHOLDERS = {
    '__RUN_NAME__': ('run_metadata', 'run_name'),
    '__SEQUENCER_ID__': ('run_metadata', 'sequencer_id'),
    '__DORADO_VERSION__': ('run_metadata', 'dorado_version'),
    '__DORADO_MODEL__': ('run_metadata', 'dorado_model'),
    '__SEQUENCING_PLATFORM__': ('run_metadata', 'sequences_platform'),  # posible typo deliberado para fallback
}

# Corrección: aceptar ambas claves sequencing_platform / sequencing_platform
ALT_KEYS = {
    ('run_metadata', 'sequences_platform'): ('run_metadata', 'sequencing_platform'),
}

def get_value(cfg, path_tuple):
    # Intentar path directo
    cur = cfg
    try:
        for p in path_tuple:
            cur = cur[p]
        return cur
    except Exception:
        # Intentar alias
        if path_tuple in ALT_KEYS:
            return get_value(cfg, ALT_KEYS[path_tuple])
        return 'UNKNOWN'

def main():
    if not TEMPLATE.exists():
        print(f"ERROR: No existe template {TEMPLATE}", file=sys.stderr)
        return 1
    if not CONFIG_YAML.exists():
        print(f"ERROR: No existe config.yaml {CONFIG_YAML}", file=sys.stderr)
        return 1

    cfg = yaml.safe_load(CONFIG_YAML.read_text()) or {}
    txt = TEMPLATE.read_text()

    for placeholder, key_path in PLACEHOLDERS.items():
        value = get_value(cfg, key_path)
        txt = txt.replace(placeholder, str(value))

    OUTPUT.write_text(txt)
    print(f"Renderizado MultiQC config -> {OUTPUT}")
    return 0

if __name__ == '__main__':
    sys.exit(main())
