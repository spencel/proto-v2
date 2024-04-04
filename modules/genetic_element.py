

import json
import logging
import os
import uuid

from devtools import debug
from django.db import models


import modules as m
import util


# Configs



bio_config = json.load(open(os.path.join("data", "bio.json"), 'r'))


class GeneticElement(models.Model):

    class Meta:
        db_table = "genetic_element"

    id = models.UUIDField(
        primary_key=True,
        default=uuid.uuid4,
    )
    name = models.CharField(max_length=255)
    type = models.CharField(max_length=255)
    length_bp = models.IntegerField()


    @classmethod
    def insert(cls,
        name = None,
        type = None,
        length_bp = None,
        genome_id = None,
        pg_cur = None
    ):
        pg_con = None
        close_con = False
        if not pg_cur:
            pg_con, pg_cur = m.Pg.get_cursor()
            close_con = True

        if not name:
            name = ""
        if not type:
            type = ""
        if not length_bp:
            length_bp = "NULL"
        if not genome_id:
            genome_id = "NULL"
            
        pg_cur.execute(
            "INSERT INTO genetic_element ("
                "name, "
                "type, "
                "length_bp"
            ") VALUES ("
                f"NULLIF('{name}', ''), "
                f"NULLIF('{type}', ''), "
                f"{length_bp}"
            ") RETURNING id;"
        )
        genetic_element_id = pg_cur.fetchone()[0]

        # Relate genome to genetic element
        pg_cur.execute(
            "INSERT INTO genome_to_genetic_element ("
                "genome_id, "
                "genetic_element_id"
            ") VALUES ("
                f"{genome_id}, "
                f"{genetic_element_id}"
            ")"
        )

        # Close DB connection if it was opened
        if close_con == True:
            pg_cur.close()
            pg_con.close()

        return genetic_element_id


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
    
    @classmethod
    def add(cls,
        name = None,
        type = None,
        length_bp = None,
        genome_id = None
    ):

        genetic_element, created_genetic_element = cls.objects.get_or_create(
            name = name,
            type = type,
            length_bp = length_bp
        )

        # Relate genome to genetic element
        # debug(genome_id)
        # debug(genetic_element.id)
        m.GenomeToGeneticElement.objects.get_or_create(
            genome_id = genome_id,
            genetic_element_id = genetic_element.id
        )

        return genetic_element.id


class GeneticElementToGene(models.Model):

    class Meta:
        db_table = "genetic_element_to_gene"

    id = models.UUIDField(
        primary_key=True,
        default=uuid.uuid4,
    )
    # genetic_element_id = models.UUIDField()
    genetic_element = models.ForeignKey("GeneticElement", on_delete=models.CASCADE)
    # gene_id = models.UUIDField()
    gene = models.ForeignKey("Gene", on_delete=models.CASCADE)