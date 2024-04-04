
from datetime import datetime
import hashlib
import json
import logging
import lzma
import os
import sys
import uuid

from devtools import debug
from django.db import models


import modules as m

# Configs
import util



bio_config = json.load(open(os.path.join("data", "bio.json"), 'r'))
root_dir = os.getcwd()
data_dir = os.path.join(root_dir, "data")
genomes_dir = os.path.join(data_dir, "genomes")
gisaid_dir = os.path.join(data_dir, "gisaid")
gisaid_samples_dir = os.path.join(gisaid_dir, "samples")


class OrganismSample(models.Model):

    class Meta:
        db_table = "organism_sample"

    id = models.UUIDField(
        primary_key=True,
        default=uuid.uuid4,
    )
    organism_sample_genome_id = models.UUIDField()
    organism_sample_source_id = models.UUIDField()

    # Returns ID of record add to organism_sample table
    @classmethod
    def insert(cls,
        data,
        source,
        reference_genome_id,
        pg_cur
    ):
        scope_pref = f"{cls.__name__}.insert"

        # logging.debug(f"{scope_pref}: source: {source}")

        if source == "gisaid":

            json_line = data

            # Don't add genetic sequence
            json_line.pop("sequence", None)

            pg_cur.execute(
                "INSERT INTO organism_sample ("
                    "meta, "
                    "reference_genome_id"
                ") VALUES ("
                    f"'{json.dumps(json_line)}', "
                    f"{reference_genome_id}"
                ") RETURNING id;"
            )
            organism_sample_id = pg_cur.fetchone()[0]
            # logging.debug(f"{scope_pref}: organism_sample_id: {organism_sample_id}")

        return organism_sample_id



    # Reconstructing sample sequence doesn't work quite right because the ends are not accounted for during variant calling.
    # A:                                                     TTGTAGATCTGTTCTCTAAACGAACTTTAAAATCTGTGTGGCTGTCACTCGGCTGCATGCTTAGTGCACTCACGCAGTAT
    # B:                                                     TTGTAGATCTGTTCTCTAAACGAACTTTAAAATCTGTGTGGCTGTCACTCGGCTGCATGCTTAGTGCACTCACGCAGTAT
    # C:                                                         AGATCTGTTCTCTAAACGAACTTTAAAATCTGTGTGGCTGTCACTCGGCTGCATGCTTAGTGCACTCACGCAGTAT
    # Ref: ATTAAAGGTTTATACCTTCCCAGGTAACAAACCAACCAACTTTCGATCTCTTGTAGATCTGTTCTCTAAACGAACTTTAAAATCTGTGTGGCTGTCACTCGGCTGCATGCTTAGTGCACTCACGCAGTATAATTAATAACTAATTACTGTCGTTGACAGGACACGAGTAACTCG
    @classmethod
    def reconstruct_na_seq(cls,
        organism_sample_id = None,
        pg_cur = None
    ):

        pg_con = None
        close_con = False
        if not pg_cur:
            pg_con, pg_cur = m.Pg.get_cursor()
            close_con = True

        # If organism sample IDs are not provided, then loop through all
        pg_cur.execute(
            "SELECT id, reference_genome_id FROM organism_sample;"
        )
        organism_samples = pg_cur.fetchall()

        for organism_sample in organism_samples:

            organism_sample_id = organism_sample[0]
            # logging.debug(f"{scope_pref}: sample_id: {sample_id}")
            reference_genome_id = organism_sample[1]
            # logging.debug(f"{scope_pref}: reference_genome_id: {reference_genome_id}")

            # Get the reference genome sequence
            ncbi_genome_version = m.Genome.get_ncbi_version(reference_genome_id)
            naseq_filepath = os.path.join(genomes_dir, f"{ncbi_genome_version}.naseq")
            f = open(naseq_filepath, 'r')

            # Get the variants in order ascending order
            pg_cur.execute(
                "SELECT "
                    "variant_segment.original_start, "
                    "variant_segment.variant_sequence "
                "FROM variant_segment_to_organism_sample "
                "JOIN variant_segment "
                f"ON variant_segment_to_organism_sample.variant_segment_id = variant_segment.id "
                f"WHERE variant_segment_to_organism_sample.organism_sample_id = {organism_sample_id} "
                "ORDER BY original_start ASC;"
            )
            res = pg_cur.fetchall()
            # logging.debug(f"{scope_pref}: res: {res}")

            # Todo: finish this method

            f.close()

        if close_con:
            pg_cur.close()
            pg_con.close()


class OrganismSampleGenome(models.Model):

    class Meta:
        db_table = "organism_sample_genome"

    id = models.UUIDField(
        primary_key=True,
        default=uuid.uuid4,
    )
    na_sequence_hash = models.TextField()


class OrganismSampleSource(models.Model):

    class Meta:
        db_table = "organism_sample_source"

    id = models.UUIDField(
        primary_key=True,
        default=uuid.uuid4,
    )
    name = models.CharField(max_length=32)


class OrganismSampleGisaid(models.Model):

    class Meta:
        db_table = "organism_sample_gisaid"

    # 11,692,139 records of GISAID data were checked in order to:
    ## Determine the data field names and or check for unexpected field names.
    ### There were no unexpected field names.
    ## Determine the max string length of each field

    # Project specific data
    id = models.UUIDField(
        primary_key=True,
        default=uuid.uuid4,
    )
    covv_accession_id = models.CharField(24)
    covv_add_location = models.TextField()
    covv_clade = models.TextField()
    covv_collection_date = models.DateField()
    covv_lineage = models.TextField()
    pangolin_lineages_version = models.TextField()
    covv_location = models.TextField()
    covv_specimen = models.TextField()
    covv_subm_date = models.DateField()
    covv_type = models.TextField()
    covv_variant = models.TextField()
    covv_virus_name = models.TextField()
    organism_sample_id = models.UUIDField()
    data_file_line_num = models.PositiveIntegerField()
    provision_file_id = models.UUIDField()
    file_position = models.BigIntegerField()
    
    # GISAID specific data
    # Total quantity of records: 11692139
    covv_accession_id = models.CharField(max_length=24)  # max character length: 18
    # covv_add_host_info =  # max character length: 284
    # covv_sampling_strategy =  # max character length: 311
    # covv_add_location =  # max character length: 193
    # covv_assembly_method =  # max character length: 477
    # covv_clade =  # max character length: 3
    # covv_collection_date =  # max character length: 10
    # covsurver_prot_mutations =  # max character length: 68884
    # gc_content =  # max character length: 0
    # covv_gender =  # max character length: 29
    # covv_host =  # max character length: 35
    # is_high_coverage =  # max character length: 0
    # is_reference =  # max character length: 0
    # is_complete =  # max character length: 0
    # covv_lineage =  # max character length: 10
    # pangolin_lineages_version =  # max character length: 15
    # covv_location =  # max character length: 434
    # n_content =  # max character length: 0
    # covv_outbreak =  # max character length: 128
    # covv_passage =  # max character length: 123
    # covv_patient_age =  # max character length: 85
    # covv_patient_status =  # max character length: 179
    # covsurver_existingmutlist =  # max character length: 28765
    # sequence =  # max character length: 35205
    # sequence_length =  # max character length: 0
    # covv_seq_technology =  # max character length: 227
    # covv_specimen =  # max character length: 231
    # covv_subm_date =  # max character length: 10
    # covv_type =  # max character length: 15
    # covsurver_uniquemutlist =  # max character length: 42742
    # covv_variant =  # max character length: 82
    # covv_virus_name =  # max character length: 92


        
    @classmethod
    def add(cls, data, line_number, provision_filepath, file_position):

        # Get the filepath ID
        file_system, file_created = m.Os.objects.get_or_create(
            filepath = provision_filepath
        )



        covv_collection_date = util.str_to_datetime(data["covv_collection_date"], "-")
        covv_subm_date = util.str_to_datetime(data["covv_subm_date"], "-")

        gisaid_sample, sample_created = m.OrganismSampleGisaid.objects.get_or_create(
            covv_accession_id         = data["covv_accession_id"],
            covv_add_location         = data["covv_add_location"],
            covv_clade                = data["covv_clade"],
            covv_collection_date      = covv_collection_date,
            covv_lineage              = data["covv_lineage"],
            pangolin_lineages_version = data["pangolin_lineages_version"],
            covv_location             = data["covv_location"],
            covv_specimen             = data["covv_specimen"],
            covv_subm_date            = covv_subm_date,
            covv_type                 = data["covv_type"],
            covv_variant              = data["covv_variant"],
            covv_virus_name           = data["covv_virus_name"],
            provision_file_id         = file_system.id,
            file_position             = file_position,
            defaults                  = { "data_file_line_num": line_number}
        )

        # # Create json.xz for this GISAID sample
        # filepath = os.path.join(
        #     gisaid_samples_dir,
        #     f"{gisaid_sample.id}.json.xz"
        # )
        # with lzma.open(filepath, "wb") as f:
        #     f.write(json.dumps(data, indent=2).encode())

        if sample_created:

            na_sequence_hash = hashlib.md5(data["sequence"].encode("utf-8")).hexdigest()

            sample_genome, c = m.OrganismSampleGenome.objects.get_or_create(
                na_sequence_hash = na_sequence_hash
            )


            sample_source = m.OrganismSampleSource.objects.get(
                name = "GISAID"
            )

            sample = m.OrganismSample.objects.create(
                organism_sample_genome_id = sample_genome.id,
                organism_sample_source_id = sample_source.id
            )

            gisaid_sample.organism_sample_id = sample.id
            gisaid_sample.save()
            
    @classmethod
    def add_samples_from_provision(cls,
        filepath = None,
        start_at_line = 0,  # 0 is the first line and default start
        stop_at_line = -2,  # -2 or else stop_at_line + 1 won't break
        interval = 1        # 1 indicates all records will be added
    ):

        f = lzma.open(filepath, 'r')

        file_position = 0
        i_interval = 1
        for i, line in enumerate(f):

            # logging.debug(f"{scope_pref}: file_position: {file_position}")

            if i < start_at_line:
                file_position = f.tell()
                continue
            elif i == stop_at_line + 1:
                break
            elif i_interval != interval:
                file_position = f.tell()
                i_interval += 1
                continue

            i_interval = 1

            sys.stdout.write(f"\r  Processed {i+1-start_at_line} lines of GISAID data.")

            # Convert line to json
            line = json.loads(line)

            m.OrganismSampleGisaid.add(
                data = line,
                line_number = i,
                provision_filepath = filepath,
                file_position = file_position
            )

            file_position = f.tell()

        sys.stdout.write("\n")


# DROP VIEW IF EXISTS view_organism_sample;
# CREATE VIEW view_organism_sample AS
# SELECT
#     organism_sample.id AS sample_id,
#     organism_sample_genome.na_sequence_hash,
#     organism_sample_source.name AS source_name
# FROM organism_sample
# JOIN organism_sample_genome
# ON organism_sample.organism_sample_genome_id = organism_sample_genome.id
# JOIN organism_sample_source
# ON organism_sample.organism_sample_source_id = organism_sample_source.id;
class ViewOrganismSample(models.Model):

    class Meta:
        db_table = "view_organism_sample"
        managed = False

    id = models.IntegerField(primary_key=True)
    sample_id = models.UUIDField()
    na_sequence_hash = models.TextField()
    source_name = models.CharField(max_length=32)


# DROP VIEW IF EXISTS view_organism_sample_gisaid;
# CREATE VIEW view_organism_sample_gisaid AS
# SELECT
#     '1'::text AS id, -- django workaround
#     organism_sample.id AS sample_id,
#     organism_sample_genome.na_sequence_hash,
#     organism_sample_source.name AS source_name,
#     organism_sample_gisaid.id AS sample_gisaid_id,
#     organism_sample_gisaid.covv_accession_id,
#     organism_sample_gisaid.data_file_line_num,
#     organism_sample_gisaid.file_position,
#     file_system.filepath
# FROM organism_sample
# JOIN organism_sample_genome
# ON organism_sample.organism_sample_genome_id = organism_sample_genome.id
# JOIN organism_sample_source
# ON organism_sample.organism_sample_source_id = organism_sample_source.id
# JOIN organism_sample_gisaid
# ON organism_sample.id = organism_sample_gisaid.organism_sample_id
# JOIN file_system
# ON organism_sample_gisaid.provision_file_id = file_system.id;
class ViewOrganismSampleGisaid(models.Model):

    class Meta:
        db_table = "view_organism_sample_gisaid"
        managed = False

    id = models.IntegerField(primary_key=True)
    sample_id = models.UUIDField()
    na_sequence_hash = models.TextField()
    source_name = models.CharField(max_length=32)
    sample_gisaid_id = models.UUIDField()
    covv_accession_id = models.CharField(max_length=24)
    data_file_line_num = models.IntegerField()
    file_position = models.BigIntegerField()
    filepath = models.CharField(max_length=255)