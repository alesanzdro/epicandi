#!/usr/bin/env python3
"""
Setup FunAnnotate database and environment variables.

Usage:
    python setup_funannotate.py <project_dir> [--genemark-dir GENEMARK_DIR]
"""

import os
import sys
import argparse
import subprocess
import urllib.request
import tarfile
import shutil
from pathlib import Path
import logging


def setup_logging():
    """Set up logging configuration."""
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.StreamHandler(sys.stdout),
            logging.FileHandler('setup_funannotate.log')
        ]
    )
    return logging.getLogger(__name__)


def run_command(cmd, check=True):
    """Run shell command and return result."""
    logger = logging.getLogger(__name__)
    logger.info(f"Running: {cmd}")

    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)

    if result.returncode != 0 and check:
        logger.error(f"Command failed: {cmd}")
        logger.error(f"Error: {result.stderr}")
        raise subprocess.CalledProcessError(result.returncode, cmd)

    return result


def setup_funannotate_db(project_dir):
    """
    Set up FunAnnotate database.

    Args:
        project_dir (Path): Project directory path

    Returns:
        Path: Database directory path
    """
    logger = logging.getLogger(__name__)

    db_dir = project_dir / "resources" / "funannotate_db"
    db_dir.mkdir(parents=True, exist_ok=True)

    logger.info(f"Setting up FunAnnotate database in: {db_dir}")

    # Set environment variable for funannotate setup
    env = os.environ.copy()
    env['FUNANNOTATE_DB'] = str(db_dir.absolute())

    # Run funannotate setup
    try:
        cmd = "funannotate setup --update --force"
        result = subprocess.run(cmd, shell=True, env=env, capture_output=True, text=True)

        if result.returncode != 0:
            logger.warning(f"FunAnnotate setup returned non-zero exit code: {result.returncode}")
            logger.warning(f"Stderr: {result.stderr}")
            logger.info("Continuing with setup despite warnings...")

        logger.info("FunAnnotate database setup completed")

    except Exception as e:
        logger.error(f"Error during funannotate setup: {e}")
        raise

    return db_dir


def download_busco_lineage(db_dir, lineage="saccharomycetes_odb10"):
    """
    Download BUSCO lineage dataset.

    Args:
        db_dir (Path): Database directory
        lineage (str): BUSCO lineage name
    """
    logger = logging.getLogger(__name__)

    lineage_url = f"https://busco-data.ezlab.org/v5/data/lineages/{lineage}.2024-01-08.tar.gz"
    lineage_file = db_dir / f"{lineage}.2024-01-08.tar.gz"

    logger.info(f"Downloading BUSCO lineage: {lineage}")

    try:
        # Download lineage dataset
        urllib.request.urlretrieve(lineage_url, lineage_file)
        logger.info(f"Downloaded: {lineage_file}")

        # Extract dataset
        with tarfile.open(lineage_file, 'r:gz') as tar:
            tar.extractall(db_dir)

        logger.info(f"Extracted BUSCO lineage dataset to: {db_dir}")

        # Clean up downloaded archive
        lineage_file.unlink()

    except Exception as e:
        logger.error(f"Error downloading BUSCO lineage: {e}")
        raise


def setup_conda_env_vars(db_dir, genemark_dir=None):
    """
    Set up conda environment variables for FunAnnotate.

    Args:
        db_dir (Path): Database directory path
        genemark_dir (Path, optional): GeneMark directory path
    """
    logger = logging.getLogger(__name__)

    conda_prefix = os.environ.get('CONDA_PREFIX')
    if not conda_prefix:
        logger.warning("CONDA_PREFIX not found. Skipping conda environment variable setup.")
        return

    activate_dir = Path(conda_prefix) / "etc" / "conda" / "activate.d"
    deactivate_dir = Path(conda_prefix) / "etc" / "conda" / "deactivate.d"

    activate_dir.mkdir(parents=True, exist_ok=True)
    deactivate_dir.mkdir(parents=True, exist_ok=True)

    activate_script = activate_dir / "funannotate_env_vars.sh"
    deactivate_script = deactivate_dir / "funannotate_env_vars.sh"

    # Create activation script
    with open(activate_script, 'w') as f:
        f.write(f"export FUNANNOTATE_DB={db_dir.absolute()}\n")

        if genemark_dir:
            f.write(f"export GENEMARK_PATH={genemark_dir.absolute()}\n")

    # Create deactivation script
    with open(deactivate_script, 'w') as f:
        f.write("unset FUNANNOTATE_DB\n")

        if genemark_dir:
            f.write("unset GENEMARK_PATH\n")

    # Make scripts executable
    activate_script.chmod(0o755)
    deactivate_script.chmod(0o755)

    logger.info(f"Conda environment variables configured:")
    logger.info(f"  Activation script: {activate_script}")
    logger.info(f"  Deactivation script: {deactivate_script}")

    if genemark_dir:
        logger.info(f"  GeneMark path: {genemark_dir}")


def fix_aux_scripts_permissions():
    """Fix permissions for FunAnnotate auxiliary scripts."""
    logger = logging.getLogger(__name__)

    conda_prefix = os.environ.get('CONDA_PREFIX')
    if not conda_prefix:
        logger.warning("CONDA_PREFIX not found. Skipping auxiliary scripts permission fix.")
        return

    # Find aux_scripts directory
    conda_lib = Path(conda_prefix) / "lib"
    aux_scripts_dirs = list(conda_lib.glob("**/funannotate/aux_scripts"))

    if not aux_scripts_dirs:
        logger.warning("FunAnnotate aux_scripts directory not found")
        return

    aux_scripts_dir = aux_scripts_dirs[0]
    logger.info(f"Found aux_scripts directory: {aux_scripts_dir}")

    # Make all scripts executable
    for script_file in aux_scripts_dir.glob("*"):
        if script_file.is_file():
            script_file.chmod(0o755)

    logger.info("Fixed permissions for FunAnnotate auxiliary scripts")


def main():
    parser = argparse.ArgumentParser(
        description="Setup FunAnnotate database and environment"
    )
    parser.add_argument("project_dir", help="Project directory path")
    parser.add_argument("--genemark-dir", help="GeneMark installation directory")
    parser.add_argument("--busco-lineage", default="saccharomycetes_odb10",
                        help="BUSCO lineage dataset to download")

    args = parser.parse_args()

    # Set up logging
    logger = setup_logging()

    try:
        project_dir = Path(args.project_dir).absolute()
        genemark_dir = Path(args.genemark_dir).absolute() if args.genemark_dir else None

        logger.info(f"Setting up FunAnnotate for project: {project_dir}")

        # Setup database
        db_dir = setup_funannotate_db(project_dir)

        # Download BUSCO lineage
        download_busco_lineage(db_dir, args.busco_lineage)

        # Setup conda environment variables
        setup_conda_env_vars(db_dir, genemark_dir)

        # Fix auxiliary scripts permissions
        fix_aux_scripts_permissions()

        logger.info("FunAnnotate setup completed successfully!")
        logger.info("Please reactivate your conda environment to load the new variables:")
        logger.info("  conda deactivate && conda activate epicandi_phylonew")

    except Exception as e:
        logger.error(f"Setup failed: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()