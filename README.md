<!--
 * @Date: 2022-10-10 15:01:33
 * @LastEditors: Hwrn
 * @LastEditTime: 2022-10-13 16:23:04
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

## compositon
- snakemake
- python
- bash
