<!--
 * @Date: 2022-10-10 15:01:33
 * @LastEditors: Hwrn
 * @LastEditTime: 2022-10-13 10:45:29
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
mamba env update -f envs/genome.yaml
cd $(dirname $(python -c "import genome; print(genome.__file__)"))/..
git pull
```

## compositon
- snakemake
- python
- bash
