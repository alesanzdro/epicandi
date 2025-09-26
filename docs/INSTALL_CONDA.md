# Instalación de Conda y Configuración del Pipeline

Este documento describe cómo instalar Conda desde cero en tu directorio home y configurar el ambiente para ejecutar el pipeline de nanocheck.

## Configuración de Proxy y Variables de Entorno

```bash
# Aseguramos tema proxy
export http_proxy=http://proxy.san.gva.es:8080
export https_proxy=http://proxy.san.gva.es:8080
export HTTP_PROXY=$http_proxy
export HTTPS_PROXY=$https_proxy
```

## Instalación de Conda

### 1. Descarga el instalador
```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-py313_25.5.1-1-Linux-x86_64.sh -O miniconda.sh
```

### 2. Lanza la instalación en tu directorio home
```bash
bash miniconda.sh -b -p ${HOME}/miniconda3
```

### 3. Inicializa conda SOLO para tu usuario
```bash
${HOME}/miniconda3/bin/conda init bash
```

### 4. Configuración proxy y canales
```bash
cat > ${HOME}/miniconda3/.condarc <<'EOF'
proxy_servers:
  http: http://proxy.san.gva.es:8080
  https: http://proxy.san.gva.es:8080
channels:
  - conda-forge
  - bioconda
  - defaults
channel_priority: strict
EOF
```

### 5. Recarga bashrc
```bash
source ${HOME}/.bashrc
```

### 6. Aceptar términos de servicio
```bash
conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main
conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r
```

### 7. Actualizar conda
```bash
conda update -n base -c defaults conda -y
```

### 8. Instalar mamba
```bash
conda install -n base -c conda-forge mamba -y
mamba update -n base mamba -y
```

### 9. Instalar git y snakemake en base
```bash
mamba install -n base -c conda-forge -c bioconda git snakemake=9.11.1 -y
```

### 10. Configurar proxy para git
```bash
git config --global http.proxy http://proxy.san.gva.es:8080
git config --global https.proxy http://proxy.san.gva.es:8080
```

### 11. Eliminar instalador
```bash
rm miniconda.sh
```

## Configuración del Environment del Pipeline

### Crear el archivo de ambiente epinano
```bash
tee envs/epinano.yml <<'EOF'
name: epinano
channels:
  - conda-forge
  - bioconda
  - defaults
dependencies:
  - python==3.12.11
  - perl==5.32.1
  - biopython==1.85
  - nanoplot==1.46.1
  - filtlong==0.3.0
  - porechop==0.2.4
  - kraken2==2.1.6
  - blast==2.16.0
  - diamond==2.1.13
  - minimap2==2.30
  - bwa==0.7.19
  - bcftools==1.22
  - bedtools==2.31.1
  - samtools==1.22.1
  - entrez-direct==22.4
  - pandas==2.3.2
  - flye==2.9.6
  - quast==5.3.0
  - canu==2.3
  - hifiasm==0.25.0
  - medaka==2.1.1
  - multiqc==1.30
EOF
```

## Ejecución del Pipeline

Una vez completada la instalación:

1. **Cargar el ambiente base** (que contiene snakemake):
   ```bash
   source ~/.bashrc
   ```

2. **Ejecutar el pipeline**:
   ```bash
   snakemake --cores 32 --use-conda --configfile config.yaml
   ```

## Notas Importantes

- La instalación de conda se realiza en tu directorio home (`${HOME}/miniconda3`)
- El ambiente base incluye snakemake 9.11.1, git y mamba
- Snakemake creará automáticamente el ambiente `epinano` la primera vez que se ejecute el pipeline
- La configuración de proxy es específica para el entorno de trabajo de la GVA

## Verificación de la Instalación

Para verificar que todo está correctamente instalado:

```bash
# Verificar conda
conda --version

# Verificar mamba
mamba --version

# Verificar snakemake
snakemake --version

# Verificar git
git --version
```