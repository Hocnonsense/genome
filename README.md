<!--
 * @Date: 2022-10-10 15:01:33
 * @LastEditors: Hwrn
 * @LastEditTime: 2022-11-24 17:11:45
 * @FilePath: /genome/README.md
 * @Description:
-->
genome
===

- genome storage and analysis unit

---
## introduction
- seperated from any unpublished database

## installation
### create environment
```bash
git clone https://github.com/Hocnonsense/genome.git
cd genome/
mamba env create -f envs/genome.yaml
conda activate genome
python setup.py develop
```

### update
```bash
conda activate genome
cd $(dirname $(python -c "import genome; print(genome.__file__)"))/..
git pull
mamba env update -f envs/genome.yaml
```

### create other conda environment
```python
from genome.create_conda_env import create_conda_env_gene_clust

create_conda_env_gene_clust()

from genome.create_conda_env import list_envs, create_conda_env

create_conda_env(*list_envs())
```

## compositon
- snakemake
- python
- bash


## changelog
- 0.0.2:
    - move functions to new packages
        - changing:
            - `checkm` from `genome.binning` to `genome.bin_statistic_ext`
            - `contig2bin` from `genome.binning` to `genome.bin_statistic`
        - now these functions can still be loaded from old packages, but cannot in next minor release.
- 0.0.1: init
