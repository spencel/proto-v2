
import json
import logging
import os
import uuid

from devtools import debug
from django.db import models


import modules as m


# Config



bio_config = json.load(open(os.path.join("data", "bio.json"), 'r'))


class VariantSegment(models.Model):

    class Meta:
        db_table = "variant_segment"

    id = models.UUIDField(
        primary_key=True,
        default=uuid.uuid4,
    )
    original_start = models.IntegerField()
    original_stop = models.IntegerField()
    original_sequence = models.TextField()
    variant_start = models.IntegerField()
    variant_end = models.IntegerField()
    variant_sequence = models.TextField()
    variant_code = models.TextField()
    mutation_type = models.IntegerField()
    mutation_type_values = {
        "deletion": 0,
        0: "deletion",
        "insertion": 1,
        1: "insertion",
        "substitution": 2,
        2: "substitution"
    }
    annotation_status = models.CharField(max_length=32)
    annotation_status_values = {
        "unannotated": 0,
        "annotated": 1
    }
    genome_severity = models.IntegerField()
    is_intergenic = models.BooleanField()
    is_coding_change = models.BooleanField()
    is_stop_codon = models.BooleanField()
    aa_position = models.IntegerField()
    aa_position_stop = models.IntegerField()
    original_aa_sequence = models.TextField()
    variant_aa_sequence = models.TextField()
    variant_aa_code = models.TextField()
    protein_severity = models.IntegerField()
    gene_id = models.IntegerField()
    gene_symbol = models.CharField(max_length=32)
    updated_at = models.DateTimeField()
    is_frameshift = models.BooleanField()
    reference_genome_id = models.UUIDField()
    pathogen_id = models.IntegerField()

    # Replace upsert() method with get or create
    @classmethod
    def add_from_lastz_output(cls,
        lastz_output_filepath = None,
        organism_sample_id = None,
        reference_genome_id = None
    ):
    
        with open(lastz_output_filepath) as lastz_output_file:
            for line in lastz_output_file:
                # Convert line to an list
                line = line.split('\t')
                # Initialize variant or mutation type
                mutation_type = "substitution"
                original_sequence = line[10]
                if original_sequence[0] == "-":
                    mutation_type = "insertion"
                original_start = int(line[1]) + 1
                variant_sequence = line[11]
                if variant_sequence[0] == "-":
                    mutation_type = "deletion"
                variant_code = f"{line[10]}{int(line[1])+1}{line[11]}"
                if line[3] == '+':
                    line[3] = 1
                else:
                    line[3] = 0
                if line[8] == '+':
                    line[8] = 1
                else:
                    line[8] = 0

                variant_segment, created = cls.objects.get_or_create(
                    variant_code = variant_code,
                    defaults = {
                        "mutation_type": cls.mutation_type_values[mutation_type],
                        "original_start": original_start,
                        "original_stop": original_start+len(original_sequence)-1,
                        "original_sequence": original_sequence,
                        "variant_sequence": variant_sequence,
                        "annotation_status": cls.annotation_status_values["unannotated"],
                        "reference_genome_id": reference_genome_id
                    }
                )

                # Insert relationship into variant_segment_to_organism_sample table
                m.VariantSegmentToOrganismSample.objects.get_or_create(
                    variant_segment_id = variant_segment.id,
                    organism_sample_id = organism_sample_id
                )

    # @classmethod
    # def upsert_old(cls,
    #     lastz_output_filepath = None,
    #     pg_cur = None,
    #     organism_sample_id = None,
    #     reference_genome_id = None
    # ):
    #     pg_con = None
    #     close_con = False
    #     if not pg_cur:
    #         close_con = True
    #         pg_con, pg_cur = m.Pg.get_cursor()
    #
    #     # Upsert variant segments to database
    #     with open(lastz_output_filepath) as lastz_output_file:
    #         for line in lastz_output_file:
    #             # Convert line to an list
    #             line = line.split('\t')
    #             # Initialize variant or mutation type
    #             mutation_type = "substitution"
    #             original_sequence = line[10]
    #             if original_sequence[0] == "-":
    #                 mutation_type = "insertion"
    #             original_start = int(line[1]) + 1
    #             variant_sequence = line[11]
    #             if variant_sequence[0] == "-":
    #                 mutation_type = "deletion"
    #             variant_code = f"{line[10]}{int(line[1])+1}{line[11]}"
    #             if line[3] == '+':
    #                 line[3] = 1
    #             else:
    #                 line[3] = 0
    #             if line[8] == '+':
    #                 line[8] = 1
    #             else:
    #                 line[8] = 0
    #             pg_cur.execute(
    #                 f"INSERT INTO variant_segment ("
    #                     f"reference_genome_id,"
    #                     f"variant_code,"
    #                     f"mutation_type, "
    #                     f"original_start,"
    #                     f"original_stop,"
    #                     # f"original_strand,"
    #                     # f"variant_start,"
    #                     # f"variant_end,"
    #                     # f"variant_strand,"
    #                     f"original_sequence,"
    #                     f"variant_sequence,"
    #                     "annotation_status"
    #                 f") VALUES ("
    #                     f"'{reference_genome_id}',"
    #                     f"'{variant_code}',"
    #                     f"'{mutation_type}',"
    #                     f"{original_start},"
    #                     f"{original_start+len(original_sequence)-1},"
    #                     # f"B\"{line[3]}\","
    #                     # f"{int(line[6])+1},"
    #                     # f"{line[7]},"
    #                     # f"B\"{line[8]}\","
    #                     f"'{original_sequence}',"
    #                     f"'{variant_sequence}',"
    #                     "'unannotated'"
    #                 f") ON CONFLICT (variant_code) DO UPDATE SET "
    #                     f"reference_genome_id = EXCLUDED.reference_genome_id, "
    #                     f"original_start = EXCLUDED.original_start, "
    #                     f"original_stop = EXCLUDED.original_stop, "
    #                     # f"original_strand = EXCLUDED.original_strand, "
    #                     # f"variant_start = EXCLUDED.variant_start, "
    #                     # f"variant_end = EXCLUDED.variant_end, "
    #                     # f"variant_strand = EXCLUDED.variant_strand, "
    #                     f"original_sequence = EXCLUDED.original_sequence, "
    #                     f"variant_sequence = EXCLUDED.variant_sequence, "
    #                     "annotation_status = EXCLUDED.annotation_status "
    #                 "RETURNING id;"
    #             )
    #             variant_segment_id = pg_cur.fetchone()[0]
    #
    #             # Insert relationship into variant_segment_to_organism_sample table
    #             pg_cur.execute(
    #                 f"INSERT INTO variant_segment_to_organism_sample ("
    #                     "variant_segment_id, "
    #                     "organism_sample_id"
    #                 ") VALUES ("
    #                     f"{variant_segment_id}, "
    #                     f"'{organism_sample_id}'"
    #                 ")"
    #             )
    #
    #     if close_con:
    #         pg_cur.close()
    #         pg_con.close()

class VariantSegmentToOrganismSample(models.Model):

    class Meta:
        db_table = "variant_segment_to_organism_sample"

    id = models.UUIDField(
        primary_key=True,
        default=uuid.uuid4,
    )
    variant_segment_id = models.UUIDField()
    organism_sample_id = models.UUIDField()