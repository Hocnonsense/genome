"""
 * @Date: 2022-10-10 15:30:31
 * @LastEditors: Hwrn hwrn.aou@sjtu.edu.cn
 * @LastEditTime: 2023-12-20 19:47:40
 * @FilePath: /genome/workflow/gene_clust.smk
 * @Description:
    use mmseq to cluster genes

useage:

"""

from genome.pyrule import mmseq_clust_95

mmseq_clust_95.register(workflow, name="gene_clust_workflow")(
    rules=["_py"], exclude_rules=[], name_modifier="mmseq_clust_95"
)
