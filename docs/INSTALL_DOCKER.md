# Docker Installation Guide

Docker is required for running InterProScan functional annotation in the EpiCandi pipeline. This is an optional component used in the robust phylogenetic analysis workflow.

## Prerequisites

- Administrative/sudo access on your system
- Internet connection for downloading Docker images

## Installation by Operating System

### Ubuntu/Debian

```bash
# Update package index
sudo apt-get update

# Install required packages
sudo apt-get install -y \
    ca-certificates \
    curl \
    gnupg \
    lsb-release

# Add Docker's official GPG key
sudo mkdir -p /etc/apt/keyrings
curl -fsSL https://download.docker.com/linux/ubuntu/gpg | sudo gpg --dearmor -o /etc/apt/keyrings/docker.gpg

# Set up the repository
echo \
  "deb [arch=$(dpkg --print-architecture) signed-by=/etc/apt/keyrings/docker.gpg] https://download.docker.com/linux/ubuntu \
  $(lsb_release -cs) stable" | sudo tee /etc/apt/sources.list.d/docker.list > /dev/null

# Update package index again
sudo apt-get update

# Install Docker Engine
sudo apt-get install -y docker-ce docker-ce-cli containerd.io docker-compose-plugin

# Start Docker service
sudo systemctl start docker
sudo systemctl enable docker
```

### CentOS/RHEL/Fedora

```bash
# Install Docker using dnf/yum
sudo dnf install -y docker

# Or for older versions:
# sudo yum install -y docker

# Start Docker service
sudo systemctl start docker
sudo systemctl enable docker
```

### macOS

1. Download Docker Desktop from: https://www.docker.com/products/docker-desktop/
2. Install the .dmg file
3. Launch Docker Desktop from Applications

### Windows

1. Download Docker Desktop from: https://www.docker.com/products/docker-desktop/
2. Install the .exe file
3. Launch Docker Desktop

## Post-Installation Configuration

### Add User to Docker Group (Linux)

To run Docker commands without sudo:

```bash
# Add your user to the docker group
sudo usermod -aG docker $USER

# Log out and log back in, or run:
newgrp docker
```

### Verify Installation

```bash
# Check Docker version
docker --version

# Test Docker installation
docker run hello-world
```

## InterProScan Docker Image

The EpiCandi pipeline uses the `blaxterlab/interproscan:latest` Docker image. This will be automatically pulled when needed, but you can pre-download it:

```bash
# Pre-download InterProScan image (optional)
docker pull blaxterlab/interproscan:latest
```

## Memory and Storage Requirements

- **Memory**: InterProScan requires significant RAM (recommended: 16GB+)
- **Storage**: Docker images and containers require ~10GB of disk space
- **Network**: Initial downloads may be large (several GB)

## Troubleshooting

### Permission Denied Errors

If you get permission denied errors:

```bash
# Ensure Docker daemon is running
sudo systemctl status docker

# Check if user is in docker group
groups $USER
```

### Docker Daemon Not Running

```bash
# Start Docker daemon
sudo systemctl start docker

# Check status
sudo systemctl status docker
```

### Network Issues

If Docker cannot pull images:

```bash
# Test network connectivity
docker run --rm alpine ping -c 3 google.com

# Check Docker daemon logs
sudo journalctl -fu docker.service
```

## Alternative: Podman

If Docker is not available, you can use Podman as an alternative:

```bash
# Install Podman (Ubuntu/Debian)
sudo apt-get install -y podman

# Create alias for docker command
alias docker=podman
```

## Usage in EpiCandi

Once Docker is installed, InterProScan will be automatically used in the pipeline when running:

```bash
# Full pipeline with InterProScan functional annotation
snakemake -s Snakefile_robust --cores 16 --use-conda robust_phylogeny
```

## Support

For Docker-specific issues:
- Docker Documentation: https://docs.docker.com/
- Docker Community Forums: https://forums.docker.com/
- InterProScan Image: https://hub.docker.com/r/blaxterlab/interproscan