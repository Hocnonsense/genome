"""
 * @Date: 2023-12-20 19:30:57
 * @LastEditors: Hwrn hwrn.aou@sjtu.edu.cn
 * @LastEditTime: 2023-12-22 16:32:30
 * @FilePath: /genome/genome/pyrule/shell_conda_rule.smk
 * @Description:
"""


rule _py:
    input:
        **config.get("inputs", {}),
    output:
        **config.get("outputs", {}),
    params:
        **config.get("params", {}),
    shadow:
        config.get("shadow")
    conda:
        config.get("conda")
    threads: config.get("threads")
    shell:
        config["shell"]
