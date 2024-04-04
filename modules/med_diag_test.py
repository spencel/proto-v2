
# med_diag_test is short for Medical Diagnostic Test


import logging
import os

from devtools import debug

import modules as m


# Configs




class MedDiagTest():

    @classmethod
    def proccess_foreign_data_immunoassay(cls,
        test_data=None,
        pg_cur=None
    ):
        data = {}
        # Collect fields
        for column_name in test_data["columns"]:
            this_column = test_data["columns"][column_name]

            if this_column["foreign_table"] == "immunoassay":
                data[this_column["foreign_column"]] = this_column["value"]
                # logging.debug(f"{scope_pref}: data: {data}")


        # Get reference genome ID
        reference_genome_id = m.Genome.get_id(genome = data["reference_genome_id"], pg_cur = pg_cur)
        data["reference_genome_id"] = reference_genome_id

        # Get gene ID
        ## Use v3 gene isoforms for now
        # target_gene_id = m.Gene.get_isoform_id_v3(gene=data["target_gene_id"], pg_cur = pg_cur)
        gene_name = data["target_gene_id"]
        # logging.debug(f"{scope_pref}: gene_name: {gene_name}")
        target_gene_id = m.GeneIsoformV3.objects.get(
            gene_name = gene_name
        ).id
        # logging.debug(f"{scope_pref}: target_gene_id: {target_gene_id}")
        data["target_gene_id"] = target_gene_id

        # Add the kit data to the immunoassay table
        immunoassay, created = m.Immunoassay.objects.get_or_create(
            target_aa_start = data["target_aa_start"],
            target_aa_end = data["target_aa_end"],
            reference_genome_id = reference_genome_id,
            target_gene_id = target_gene_id
        )

        return immunoassay
