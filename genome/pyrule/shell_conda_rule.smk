"""
 * @Date: 2023-12-20 19:30:57
 * @LastEditors: Hwrn hwrn.aou@sjtu.edu.cn
 * @LastEditTime: 2023-12-20 20:19:54
 * @FilePath: /genome/genome/pyrule/shell_conda_rule.smk
 * @Description:
"""


rule _py:
    input:
        **config["inputs"],
    output:
        **config["outputs"],
    params:
        **config["params"],
    shadow:
        config["shadow"]
    conda:
        config["conda"]
    threads: config["threads"]
    shell:
        config["shell"]
