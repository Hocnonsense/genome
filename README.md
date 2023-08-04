<!--
 * @Date: 2022-10-10 15:01:33
 * @LastEditors: Hwrn hwrn.aou@sjtu.edu.cn
 * @LastEditTime: 2023-08-04 11:07:34
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

- to set database for checkm, you should specific a directory (e.g., ~/.checkm):

```bash
export CHECKM_DATA_PATH=~/.checkm
wget https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz
tar -zxvf checkm_data_2015_01_16.tar.gz -C $CHECKM_DATA_PATH
checkm data setRoot $CHECKM_DATA_PATH
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
- 0.1.1:
    - add mantis rules for snakemake in `genome.pyrule`, only for snakemake>=7.31
    - remove `checkm` and `contig2bin` from `genome.binning`
- 0.1.0:
    - [Version match with meta-snakemake-minimal v0.1.0](http://202.120.45.162:12080/Metabolic_Modeling/genome/releases/tag/version-0.1.0)
- 0.0.3:
    - remove dependence of bcbio-gff (no other change)
- 0.0.2:
    - move functions to new packages
        - changing:
            - `checkm` from `genome.binning` to `genome.bin_statistic_ext`
            - `contig2bin` from `genome.binning` to `genome.bin_statistic`
        - now these functions can still be loaded from old packages, but cannot in next minor release.
- 0.0.1: init
