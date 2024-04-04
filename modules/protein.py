

import json
import logging
import os
import uuid

from django.db import models


import modules as m
import util


# Configs



bio_config = json.load(open(os.path.join("data", "bio.json"), 'r'))


class Protein(models.Model):
    
    class Meta:
        db_table = "protein"

    id = models.UUIDField(
        primary_key=True,
        default=uuid.uuid4,
    )
    location = models.TextField()
    gene_symbol = models.CharField(max_length=255)
    locus_tag = models.CharField(max_length=255)
    product = models.TextField()
    note = models.TextField()
    ncbi_protein_id = models.CharField(max_length=255)
    ncbi_db_xref = models.CharField(max_length=255)
    genbank_id = models.CharField(max_length=255)
    gene_cds_id = models.UUIDField()

    # @classmethod
    # def insert_old(cls,
    #     data = None,
    #     source_type = None,
    #     gene_cds_id = None,
    #     pg_cur = None
    # ):
    #     pg_con = None
    #     close_con = False
    #     if not pg_cur:
    #         pg_con, pg_cur = m.Pg.get_cursor()
    #         close_con = True
    #
    #     if source_type == "GenBank":
    #         feature = data
    #
    #         location = ""
    #         if "location" in vars(feature):
    #             location = util.json_to_str(feature.location)
    #         genbank_id = ""
    #         if "id" in vars(feature):
    #             genbank_id = m.Pg.clean_text_value(feature.id)
    #
    #         qualifiers = feature.qualifiers
    #         gene_symbol = ""
    #         if "gene" in qualifiers:
    #             gene_symbol = m.Pg.clean_text_value(qualifiers["gene"])
    #         locus_tag = ""
    #         if "locus_tag" in qualifiers:
    #             locus_tag = m.Pg.clean_text_value(qualifiers["locus_tag"])
    #         product = ""
    #         if "product" in qualifiers:
    #             product = m.Pg.clean_text_value(qualifiers["product"])
    #         note = ""
    #         if "note" in qualifiers:
    #             note = m.Pg.clean_text_value(qualifiers["note"])
    #         ncbi_protein_id = ""
    #         if "protein_id" in qualifiers:
    #             ncbi_protein_id = m.Pg.clean_text_value(qualifiers["protein_id"])
    #         ncbi_db_xref = ""
    #         if "db_xref" in qualifiers:
    #             ncbi_db_xref = m.Pg.clean_text_value(qualifiers["db_xref"])
    #
    #         if not gene_cds_id:
    #             gene_cds_id = "NULL"
    #
    #
    #     pg_cur.execute(
    #         "INSERT INTO protein ("
    #             "location, "
    #             "gene_symbol, "
    #             "locus_tag, "
    #             "product, "
    #             "note, "
    #             "ncbi_protein_id, "
    #             "ncbi_db_xref, "
    #             "genbank_id, "
    #             "gene_cds_id"
    #         ") VALUES ("
    #             f"NULLIF('{location}', ''), "
    #             f"NULLIF('{gene_symbol}', ''), "
    #             f"NULLIF('{locus_tag}', ''), "
    #             f"NULLIF('{product}', ''), "
    #             f"NULLIF('{note}', ''), "
    #             f"NULLIF('{ncbi_protein_id}', ''), "
    #             f"NULLIF('{ncbi_db_xref}', ''), "
    #             f"NULLIF('{genbank_id}', ''), "
    #             f"{gene_cds_id}"
    #         ") RETURNING id;"
    #     )
    #     gene_cds_id = pg_cur.fetchone()[0]
    #
    #     # Close DB connection if it was opened
    #     if close_con == True:
    #         pg_cur.close()
    #         pg_con.close()
    #
    #     return gene_cds_id

    @classmethod
    def insert(cls,
        data = None,
        source_type = None,
        gene_cds_id = None
    ):

        if source_type == "GenBank":
            feature = data

            location = None
            if "location" in vars(feature):
                location = util.json_to_str(feature.location)
            genbank_id = None
            if "id" in vars(feature):
                genbank_id = m.Pg.clean_text_value(feature.id)

            qualifiers = feature.qualifiers
            gene_symbol = None
            if "gene" in qualifiers:
                gene_symbol = m.Pg.clean_text_value(qualifiers["gene"])
            locus_tag = None
            if "locus_tag" in qualifiers:
                locus_tag = m.Pg.clean_text_value(qualifiers["locus_tag"])
            product = None
            if "product" in qualifiers:
                product = m.Pg.clean_text_value(qualifiers["product"])
            note = None
            if "note" in qualifiers:
                note = m.Pg.clean_text_value(qualifiers["note"])
            ncbi_protein_id = None
            if "protein_id" in qualifiers:
                ncbi_protein_id = m.Pg.clean_text_value(qualifiers["protein_id"])
            ncbi_db_xref = None
            if "db_xref" in qualifiers:
                ncbi_db_xref = m.Pg.clean_text_value(qualifiers["db_xref"])

            if not gene_cds_id:
                gene_cds_id = None


        gene_cds, created = cls.objects.get_or_create(
            location = location, 
            gene_symbol = gene_symbol,
            locus_tag = locus_tag,
            product = product,
            note = note,
            ncbi_protein_id = ncbi_protein_id,
            ncbi_db_xref = ncbi_db_xref,
            genbank_id = genbank_id,
            gene_cds_id = gene_cds_id
        )

        return gene_cds.id

    @classmethod
    def get_id(cls, genome, pg_cur):
        scope_pref = f"{cls.__name__}.get_id"
        # Args
        #  genome: can be genome filepath or NCBI accession and version
        # Returns
        #  A single ID of the genome in this table

        # Try to get by NCBI accession and version
        pg_cur.execute(f"SELECT id FROM genome WHERE ncbi_version = '{genome}';")
        res = pg_cur.fetchone()
        if res:
            genome_id = res[0]
            return genome_id

        # Try to get by genome filepath
        pg_cur.execute(f"SELECT id FROM genome WHERE fasta_filename = '{genome}';")
        res = pg_cur.fetchone()
        if res:
            genome_id = res[0]
            return genome_id

    @classmethod
    def get_row(cls, id, pg_cur):
        scope_pref = f"{cls.__name__}.get_id"
        pg_cur.execute(f"SELECT * FROM genome WHERE id = {id}")
        res = pg_cur.fetchone()
        # logging.debug(f"{scope_pref}: res: {res}")
        if res:
            results = dict()
            for i, column in enumerate(cls.columns):
                results[column] = res[i]
            results = util.dotdict(results)
            return results

    @classmethod
    def reconstruct_sequence_alignment(cls,
        reference_sequence = None, # The reference nucleic acid sequence
        reference_offset_origin = None, # Where the reference sequence starts if it's a subsection of the whole sequence
        reference_segments = None, # A list of these
        variant_segments   = None, # A list of these
        variant_segment_start_loci_on_reference_sequence = None, # A list of these
        mutation_types = None # A list of these
    ):
        scope_pref = __name__

        # Case 1
        if reference_sequence and variant_segments and variant_segment_start_loci_on_reference_sequence and mutation_types:

            # Handle arguments
            na_sequence = reference_sequence
            tech_seq_start = reference_offset_origin
            ref_sequences = reference_segments
            var_sequences = variant_segments
            ref_starts = variant_segment_start_loci_on_reference_sequence
            mut_types = mutation_types

            # Loop through variant segments and create reference and variant technical sequence alignment
            # Set the offset to account for insertions
            ref_offset = 0
            var_offset = 0
            # Initialize nucleotide sequences
            ref_tech_sequence = na_sequence
            var_tech_sequence = na_sequence
            for i in range(len(ref_starts)):
                # logging.debug(f"{scope_pref}: i: {i}")
                # Alias variant start relative to technical sequence, not start on variant genome
                ref_start = ref_starts[i] - tech_seq_start + ref_offset
                var_start = ref_starts[i] - tech_seq_start + var_offset
                ref_sequence = ref_sequences[i]
                var_sequence = var_sequences[i]
                mut_type = mut_types[i]
                # logging.debug(f"{scope_pref}: ref_sequence: {ref_sequence}")
                # logging.debug(f"{scope_pref}: var_sequence: {var_sequence}")

                # Handle insertion
                if mut_type == bio_config["mutation_types"]["insertion"]:
                    # logging.debug(f"{scope_pref}: insertion")
                    ref_tech_sequence = ref_tech_sequence[:ref_start] + ref_sequence + ref_tech_sequence[ref_start:]
                    # logging.debug(f"{scope_pref}: ref_tech_sequence: {ref_tech_sequence}")
                    var_tech_sequence = var_tech_sequence[:var_start] + var_sequence + var_tech_sequence[var_start:]
                    # logging.debug(f"{scope_pref}: var_tech_sequence: {var_tech_sequence}")
                    ref_offset += len(ref_sequence)

                # Handle deletion
                if mut_type == bio_config["mutation_types"]["deletion"]:
                    # logging.debug(f"{scope_pref}: deletion")
                    var_tech_sequence = var_tech_sequence[:ref_start] + var_sequence + var_tech_sequence[ref_start + len(var_sequence):]
                    # logging.debug(f"{scope_pref}: var_tech_sequence: {var_tech_sequence}")

                # Handle substitution
                elif mut_type == bio_config["mutation_types"]["substitution"]:
                    # logging.debug(f"{scope_pref}: substitution")
                    var_tech_sequence = \
                          var_tech_sequence[:ref_start] \
                        + var_sequence \
                        + var_tech_sequence[ref_start + len(var_sequence):]
                    # logging.debug(f"{scope_pref}: var_tech_sequence: {var_tech_sequence}")

            return ref_tech_sequence, var_tech_sequence