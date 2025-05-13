# -*- coding: utf-8 -*-
"""
* @Date: 2025-05-13 22:03:08
* @LastEditors: hwrn hwrn.aou@sjtu.edu.cn
* @LastEditTime: 2025-05-13 22:22:27
* @FilePath: /genome/tests/genome/test_gene_annot.py
* @Description:
"""
# """

from genome import gene_annot

from tests import Path, temp_output, test_files, test_temp


class TestMantisAnnot:
    @staticmethod
    def make_file():
        lines = (
            "Query\tRef_Files\tRef_Hits\tConsensus_hits\tTotal_hits\t|\tLinks\n"
            "BCX55083.1\tNOGG_merged;Pfam-A\t688245.CtCNB1_4570;RHH_4\t2\t2\t|\tcog:COG4321\tdescription:Ribbon-helix-helix domain\teggnog:1MZP2\teggnog:2VU7X\teggnog:4AEFF\teggnog:COG4321\tpfam:PF13467\tpfam:RHH_4\n"
            "BCX53579.1\tNOGG_merged;Pfam-A;kofam_merged\t688245.CtCNB1_1751;tRNA_edit;K03976\t3\t4\t|\tcog:COG2606\tdescription:Aminoacyl-tRNA editing domain\tdescription:Cys-tRNA(Pro)/Cys-tRNA(Cys) deacylase\tdescription:aminoacyl-tRNA editing activity\teggnog:1RGX5\teggnog:2VQAY\teggnog:4ABKA\teggnog:COG2606\tenzyme_ec:3.1.1.-\tgo:0002161\tgo:0043907\tkegg_ko:K03976\tpfam:PF04073\tpfam:tRNA_edit\n"
            "BCX54275.1\tNOGG_merged;Pfam-A;kofam_merged\t688245.CtCNB1_1153;Aldedh;K21802\t3\t4\t|\tcog:COG1012\tdescription:Aldehyde dehydrogenase family\tdescription:oxidoreductase activity\tdescription:vanillin dehydrogenase\teggnog:1MU1V\teggnog:2VH71\teggnog:4ABR0\teggnog:COG1012\tenzyme_ec:1.2.1.67\tgo:0016491\tgo:0050608\tkegg_ko:K21802\tpfam:Aldedh\tpfam:PF00171\n"
            "BCX50661.1\tNOGG_merged;Pfam-A;Pfam-A;kofam_merged\t688245.CtCNB1_0252;HTH_1;LysR_substrate;K09681\t4\t5\t|\tcog:COG0583\tdescription:Bacterial regulatory helix-turn-helix protein, lysR family\tdescription:DNA-binding transcription factor activity\tdescription:LysR family transcriptional regulator, transcription activator of glutamate synthase operon\tdescription:LysR substrate binding domain\tdescription:regulation of transcription, DNA-templated\teggnog:1NUAB\teggnog:2VP1P\teggnog:4AF5A\teggnog:COG0583\tgo:0003700\tgo:0006355\tkegg_ko:K09681\tpfam:HTH_1\tpfam:LysR_substrate\tpfam:PF00126\tpfam:PF03466\n"
        )
        with open(test_temp / "mantis.txt", "w") as text:
            print(lines.strip(), file=text)
        return Path(text.name)

    def test_get_ref_annots(self):
        file = self.make_file()
        annots = gene_annot.MantisAnnot.get_ref_annots(file)
        assert annots.to_csv() == (
            ",NOGG_merged,Pfam-A,kofam_merged\n"
            "BCX55083.1,['688245.CtCNB1_4570'],['RHH_4'],\n"
            "BCX53579.1,['688245.CtCNB1_1751'],['tRNA_edit'],['K03976']\n"
            "BCX54275.1,['688245.CtCNB1_1153'],['Aldedh'],['K21802']\n"
            "BCX50661.1,['688245.CtCNB1_0252'],\"['HTH_1', 'LysR_substrate']\",['K09681']\n"
        )

    def test_get_link_annots(self):
        file = self.make_file()
        annots = gene_annot.MantisAnnot.get_link_annots(file)
        assert annots.to_csv() == (
            """,KO,COG,arCOG,PFAM,EC,description,eggnog,go\n"""
            """BCX55083.1,[],['COG4321'],[],"['PF13467', 'RHH_4']",[],['Ribbon-helix-helix domain'],"['1MZP2', '2VU7X', '4AEFF', 'COG4321']",\n"""
            """BCX53579.1,['K03976'],['COG2606'],[],"['PF04073', 'tRNA_edit']",['3.1.1.-'],"['Aminoacyl-tRNA editing domain', 'Cys-tRNA(Pro)/Cys-tRNA(Cys) deacylase', 'aminoacyl-tRNA editing activity']","['1RGX5', '2VQAY', '4ABKA', 'COG2606']","['0002161', '0043907']"\n"""
            """BCX54275.1,['K21802'],['COG1012'],[],"['Aldedh', 'PF00171']",['1.2.1.67'],"['Aldehyde dehydrogenase family', 'oxidoreductase activity', 'vanillin dehydrogenase']","['1MU1V', '2VH71', '4ABR0', 'COG1012']","['0016491', '0050608']"\n"""
            """BCX50661.1,['K09681'],['COG0583'],[],"['HTH_1', 'LysR_substrate', 'PF00126', 'PF03466']",[],"['Bacterial regulatory helix-turn-helix protein, lysR family', 'DNA-binding transcription factor activity', 'LysR family transcriptional regulator, transcription activator of glutamate synthase operon', 'LysR substrate binding domain', 'regulation of transcription, DNA-templated']","['1NUAB', '2VP1P', '4AF5A', 'COG0583']","['0003700', '0006355']"\n"""
        )
