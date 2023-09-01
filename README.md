<!--
 * @Date: 2022-10-10 15:01:33
 * @LastEditors: Hwrn hwrn.aou@sjtu.edu.cn
 * @LastEditTime: 2023-08-07 15:36:58
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
- [changelog](changelog.md)
