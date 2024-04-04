
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

class GenomeBase():

    strand_values = {
        "reference": 0,
        "complement": 1,
        "exotic": 2
    }

class Genome(models.Model):

    class Meta:
        db_table = "genome"

    id = models.UUIDField(
        primary_key=True,
        default=uuid.uuid4,
    )
    dir_name = models.CharField(255)
    # Todo: move to a GenBank database table and rename to genbank_version?
    ncbi_version = models.CharField(max_length=16)
    taxonomy_id = models.UUIDField()

    columns = [
        "id",
        "ncbi_version",
        "dir_name",
        "taxonomy_id"
    ]

    @classmethod
    def insert(cls,
        fasta_filename = None,
        ncbi_version = None,
        taxonomy_id = None,
        pg_cur = None
    ):
        pg_con = None
        close_con = False
        if not pg_cur:
            pg_con, pg_cur = m.Pg.get_cursor()
            close_con = True

        if not fasta_filename:
            fasta_filename = ""
        if not ncbi_version:
            ncbi_version = ""
        if not taxonomy_id:
            taxonomy_id = "NULL"
            
        pg_cur.execute(
            "INSERT INTO genome ("
                "fasta_filename, "
                "ncbi_version, "
                "taxonomy_id"
            ") VALUES ("
                f"'{fasta_filename}', "
                f"'{ncbi_version.upper()}', "
                f"{taxonomy_id}"
            ") RETURNING id;"
        )
        genome_id = pg_cur.fetchone()[0]

        # Close DB connection if it was opened
        if close_con == True:
            pg_cur.close()
            pg_con.close()

        return genome_id


    @classmethod
    def get_ncbi_version(cls, genome_id, pg_cur = None):
        scope_pref = f"{cls.__name__}.get_ncbi_version"

        pg_con = None
        close_con = False
        if not pg_cur:
            pg_con, pg_cur = m.Pg.get_cursor()
            close_con = True

        pg_cur.execute(f"SELECT ncbi_version FROM genome WHERE id = {genome_id}")
        ncbi_version = pg_cur.fetchone()[0]
        
        if close_con:
            pg_cur.close()
            pg_con.close()

        return ncbi_version


    @classmethod
    def get_id(cls, genome, pg_cur):
        scope_pref = f"{cls.__name__}.get_id"
        # Args
        #  genome: can be genome filepath or NCBI accession and version
        # Returns
        #  A single ID of the genome in this table

        # Try to get by NCBI version
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
        scope_pref = f"{cls.__name__}.get_row"
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
    def load_gbk(cls,
        ncbi_version,
        genetic_element,
        genetic_element_type,
        lineage = None,
        feature_qty_limit = None
    ):

        fasta_filename = os.path.join(ncbi_version, ".fa")

        genome, genome_created = m.Genome.objects.get_or_create(
            dir_name = ncbi_version,
            ncbi_version = ncbi_version
        )

        genetic_element_id = m.GeneticElement.add(
            name = genetic_element,
            type = genetic_element_type,
            genome_id = genome.id
        )

        m.GenBank.load_gbk_file_into_db(
            ncbi_version=ncbi_version,
            genetic_element_id=genetic_element_id,
            genome_id = genome.id,
            feature_qty_limit = feature_qty_limit
        )

class GenomeToGeneticElement(models.Model):

    class Meta:
        db_table = "genome_to_genetic_element"

    id = models.UUIDField(
        primary_key=True,
        default=uuid.uuid4,
    )
    genome_id = models.UUIDField()
    genetic_element_id = models.UUIDField()

# Traverse tables: genome > m2m > genetic_element > m2m > gene
# DROP VIEW IF EXISTS view_genome_to_gene;
# CREATE VIEW view_genome_to_gene AS
# SELECT
#     '1'::text AS id, -- django workaround
#     genome.id AS genome_id,
#     genome.ncbi_version,
#     genetic_element.name AS genetic_element_name,
#     genetic_element.type AS genetic_element_type,
#     gene.id AS gene_id,
#     gene.alias,
#     gene.description,
#     gene.symbol,
#     gene.ensembl_gene_id,
#     gene.ncbi_gene_id
# FROM genome
# JOIN genome_to_genetic_element
# ON genome.id = genome_to_genetic_element.genome_id
# JOIN genetic_element
# ON genome_to_genetic_element.genetic_element_id = genetic_element.id
# JOIN genetic_element_to_gene
# ON genetic_element.id = genetic_element_to_gene.genetic_element_id
# JOIN gene
# ON genetic_element_to_gene.gene_id = gene.id;
class ViewGenomeToGene(models.Model):

    class Meta:
        db_table = "view_genome_to_gene"
        managed = False

    id = models.IntegerField(primary_key=True)
    genome_id = models.UUIDField()
    ncbi_version = models.CharField(max_length=255)
    genetic_element_name = models.CharField(max_length=255)
    genetic_element_type = models.CharField(max_length=255)
    gene_id = models.UUIDField()
    alias = models.CharField(max_length=4095)
    description = models.CharField(max_length=4095)
    symbol = models.CharField(max_length=255)
    ensembl_gene_id = models.CharField(max_length=255)
    ncbi_gene_id = models.BigIntegerField()
