<!--
 * @Date: 2022-10-10 15:01:33
 * @LastEditors: hwrn hwrn.aou@sjtu.edu.cn
 * @LastEditTime: 2024-04-28 20:08:38
 * @FilePath: /genome/README.md
 * @Description:
-->
genome
===

- genome storage and analysis unit

---
## introduction
- mantains some snakemake modules, and can also import python functions to handle the output.

## installation
### create environment
the main environment is defined [here](genome/pyrule/envs/genome.yaml)
```bash
mamba env create -p .snakemake/conda/snakemake-g -f https://raw.githubusercontent.com/Hocnonsense/genome/master/genome/pyrule/envs/genome.yaml
conda activate .snakemake/conda/snakemake-g
pip install git+https://github.com/Hocnonsense/genome.git
```

or with a specific version:
```bash
version="0.2.3"
mamba env create -p snakemake-g-$version -f https://raw.githubusercontent.com/Hocnonsense/genome/$version/genome/pyrule/envs/genome.yaml
conda activate snakemake-g-$version
pip install git+https://github.com/Hocnonsense/genome.git@$version
```

some post-deploy settings are defined [here](genome/pyrule/envs/genome.post-deploy.sh)

to set database for checkm, you should specific a directory (e.g., ~/.checkm):
```bash
export CHECKM_DATA_PATH=~/.checkm
wget https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz
tar -zxvf checkm_data_2015_01_16.tar.gz -C $CHECKM_DATA_PATH
checkm data setRoot $CHECKM_DATA_PATH
```

### update
it is recommanded to create a new conda environment

### create other conda environment
```python
from genome.create_conda_env import create_conda_env_gene_clust

create_conda_env_gene_clust()
```

```python
from genome.create_conda_env import list_envs, create_conda_env

create_conda_env(*list_envs())
```

## compositon
- snakemake
- python
- bash

## changelog
- [changelog](changelog.md)


## tests
```bash
# normal tests
python -m pytest -vv --cov --cov-report=html

# also run tests with computation-heavily functions e.g. checkm and binning
python -m pytest -vv --cov --cov-report=html -m "not limit_resource"
```
