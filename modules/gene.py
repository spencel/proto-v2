
import logging
import os
import uuid

from devtools import debug
from django.db import models


import modules as m
import util


# Configs





class Gene(models.Model):

    class Meta:
        db_table = "gene"

    id = models.UUIDField(        # UUID PRIMARY KEY DEFAULT uuid_generate_v4 (),
        primary_key=True,
        default=uuid.uuid4,
    )
    location = models.TextField()            # TEXT,
    locus_tag = models.CharField(max_length=255)           # VARCHAR (255),
    symbol = models.CharField(max_length=255)              # VARCHAR (255), -- same as gene qualifier in gbk files
    alias = models.CharField(max_length=4095)               # VARCHAR (4095),
    gene_synonym = models.TextField()        # TEXT,
    description = models.CharField(max_length=4095)         # VARCHAR (4095),
    ensembl_gene_id = models.CharField(max_length=255)     # VARCHAR (255),
    ncbi_gene_id = models.BigIntegerField()        # BIGINT,
    ncbi_db_xref = models.CharField(max_length=255)        # VARCHAR (255),
    genbank_id = models.CharField(max_length=255)          # VARCHAR (255)


    columns = [
        "id",
        "alias",
        "description",
        "symbol",
        "ensembl_gene_id",
        "ncbi_gene_id"
    ]

    # @classmethod
    # def insert_old(cls,
    #     data = None,
    #     source_type = None,
    #     genetic_element_id = None,
    #     pg_cur = None
    # ):
    #     scope_pref = f"{cls.__name__}.insert"
    # 
    # 
    #     ### Handle GenBank Data #############################################################################################
    # 
    #     if source_type == "GenBank":
    #         # logging.debug(f"{scope_pref}: data: {data}")
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
    #         locus_tag = ""
    #         if "locus_tag" in qualifiers:
    #             locus_tag = m.Pg.clean_text_value(qualifiers["locus_tag"])
    #         symbol = ""
    #         if "gene" in qualifiers:
    #             symbol = m.Pg.clean_text_value(qualifiers["gene"])
    #         alias = ""
    #         gene_synonym = ""
    #         if "gene_synonym" in qualifiers:
    #             gene_synonym = m.Pg.clean_text_value(qualifiers["gene_synonym"])
    #         description = ""
    #         enemble_gene_id = ""
    #         ncbi_gene_id = "NULL"
    #         ncbi_db_xref = ""
    #         if "db_xref" in qualifiers:
    #             ncbi_db_xref = m.Pg.clean_text_value(qualifiers["db_xref"])
    #         genbank_id = ""
    #         if "genbank_id" in qualifiers:
    #             genbank_id = m.Pg.clean_text_value(qualifiers["genbank_id"])
    # 
    # 
    #     ### Handle list/row Data ############################################################################################
    # 
    #     elif isinstance(data, list):
    #         # logging.debug(f"{scope_pref}: data: {data}")
    # 
    #         location = ""
    #         locus_tag = ""
    #         symbol = data[2]
    #         alias = data[0]
    #         gene_synonym = ""
    #         description = data[1]
    #         enemble_gene_id = data[3]
    #         ncbi_gene_id = data[4]
    #         if ncbi_gene_id == "":
    #             ncbi_gene_id = "NULL"
    #         ncbi_db_xref = ""
    #         genbank_id = ""
    # 
    #     else:
    #         logging.error(f"{scope_pref}: Missing code to handle this type of data.")
    #         raise Exception
    # 
    # 
    #     ### Insert ##########################################################################################################
    # 
    #     pg_cur.execute(
    #         f"INSERT INTO gene ("
    #             "location, "
    #             "locus_tag, "
    #             "symbol, "
    #             "alias, "
    #             "gene_synonym, "
    #             "description, "
    #             "ensembl_gene_id, "
    #             "ncbi_gene_id, "
    #             "ncbi_db_xref, "
    #             "genbank_id"
    #         ") VALUES ("
    #             f"NULLIF('{location}', ''), "
    #             f"NULLIF('{locus_tag}', ''), "
    #             f"NULLIF('{symbol.upper()}', ''), "
    #             f"NULLIF('{alias}', ''), "
    #             f"NULLIF('{gene_synonym}', ''), "
    #             f"NULLIF('{description}', ''), "
    #             f"NULLIF('{enemble_gene_id}', ''), "
    #             f"{ncbi_gene_id}, "
    #             f"NULLIF('{ncbi_db_xref}', ''), "
    #             f"NULLIF('{genbank_id}', '')"
    #         ") RETURNING id;"
    #     )
    #     gene_id =  pg_cur.fetchone()[0]
    # 
    #     # Relate genetic element to gene
    #     pg_cur.execute(
    #         "INSERT INTO genetic_element_to_gene ("
    #             "genetic_element_id, "
    #             "gene_id"
    #         ") VALUES ("
    #             f"{genetic_element_id}, "
    #             f"{gene_id}"
    #         ")"
    #     )
    # 
    #     return gene_id
    
    @classmethod
    def insert(cls,
        data = None,
        source_type = None,
        genetic_element_id = None
    ):
        scope_pref = f"{cls.__name__}.insert"


        ### Handle GenBank Data #############################################################################################

        if source_type == "GenBank":
            # logging.debug(f"{scope_pref}: data: {data}")
            feature = data

            location = None
            if "location" in vars(feature):
                location = util.json_to_str(feature.location)
            feature_id = None
            if "id" in vars(feature):
                feature_id = m.Pg.clean_text_value(feature.id)

            qualifiers = feature.qualifiers
            locus_tag = None
            if "locus_tag" in qualifiers:
                locus_tag = m.Pg.clean_text_value(qualifiers["locus_tag"])
            symbol = None
            if "gene" in qualifiers:
                symbol = m.Pg.clean_text_value(qualifiers["gene"]).upper()
            # logging.debug(f"{scope_pref}: symbol: {symbol}")
            alias = None
            gene_synonym = None
            if "gene_synonym" in qualifiers:
                gene_synonym = m.Pg.clean_text_value(qualifiers["gene_synonym"])
            description = None
            enembl_gene_id = None
            ncbi_gene_id = None
            ncbi_db_xref = None
            if "db_xref" in qualifiers:
                ncbi_db_xref = m.Pg.clean_text_value(qualifiers["db_xref"])
            genbank_id = None
            if "genbank_id" in qualifiers:
                genbank_id = m.Pg.clean_text_value(qualifiers["genbank_id"])


        ### Handle list/row Data ############################################################################################

        elif isinstance(data, list):
            # logging.debug(f"{scope_pref}: data: {data}")

            location = None
            locus_tag = None
            symbol = data[2]
            # logging.debug(f"{scope_pref}: symbol: {symbol}")
            alias = data[0]
            gene_synonym = None
            description = data[1]
            enembl_gene_id = data[3]
            ncbi_gene_id = data[4]
            ncbi_db_xref = None
            genbank_id = None

        else:
            logging.error(f"{scope_pref}: Missing code to handle this type of data.")
            raise Exception


        ### Insert ##########################################################################################################

        gene, created = cls.objects.get_or_create(
            location = location,
            locus_tag = locus_tag,
            symbol = symbol,
            alias = alias,
            gene_synonym = gene_synonym,
            description = description,
            ensembl_gene_id = enembl_gene_id,
            ncbi_gene_id = ncbi_gene_id,
            ncbi_db_xref = ncbi_db_xref,
            genbank_id = genbank_id
        )

        # Relate genetic element to gene
        m.GeneticElementToGene.objects.get_or_create(
            genetic_element_id = genetic_element_id,
            gene_id = gene.id
        )

        return gene.id


    @classmethod
    def get_gene_and_genome_id(cls,
        gene_symbol=None,
        ref_genome=None,
        pg_cur=None
    ):

        pg_cur.execute(
            "SELECT gene_id, genome_id "
            "FROM view_genome_to_gene "
            "WHERE "
                f"ncbi_version = '{ref_genome.upper()}' "
                f"AND symbol = '{gene_symbol.upper()}'"
        )
        return pg_cur.fetchone()


    @classmethod
    def get_isoform_id_v3(cls,
        gene = None,
        pg_cur = None
    ):

        pg_cur.execute(
            f"SELECT id FROM gene_isoform_v3 WHERE gene_name = '{gene}'"
        )
        return pg_cur.fetchone()[0]

    @classmethod
    def get_isoform_symbol_v3(cls,
        gene = None,
        pg_cur = None
    ):

        pg_cur.execute(
            f"SELECT gene_name FROM gene_isoform_v3 WHERE id = {gene}"
        )
        return pg_cur.fetchone()[0]


class GeneCds(models.Model):

    class Meta:
        db_table = "gene_cds"

    id = models.UUIDField(        # UUID PRIMARY KEY DEFAULT uuid_generate_v4 (),
        primary_key=True,
        default=uuid.uuid4,
    )
    location = models.TextField()
    gene = models.CharField(max_length=255)
    locus_tag = models.CharField(max_length=255)
    ribosomal_slippage = models.BooleanField()
    gene_synonym = models.CharField(max_length=255)
    note = models.TextField()
    codon_start = models.IntegerField()
    product = models.CharField(max_length=255)
    ncbi_protein_id = models.CharField(max_length=255)
    ncbi_db_xref = models.CharField(max_length=255)
    translation = models.TextField()
    gene_id = models.UUIDField()
    genbank_id = models.CharField(max_length=255)

    # @classmethod
    # def insert_cds_old(cls,
    #                data=None,
    #                source_type=None,
    #                gene_id=None,
    #                pg_cur=None
    #                ):
    #     scope_pref = f"{cls.__name__}.insert_cds"
    #
    #     # # Todo: handle lists with more than 1 item
    #     # ## Could use list to string utility
    #     # for key in data:
    #     #     if len(data[key]) > 1:
    #     #         logging.debug(f"{scope_pref}: key: {key}")
    #     #         logging.debug(f"{scope_pref}: value: {data[key]}")
    #     #         raise
    #     if source_type == "GenBank":
    #         feature = data
    #
    #         if not gene_id:
    #             gene_id = "NULL"
    #
    #         location = ""
    #         if "location" in vars(feature):
    #             location = util.json_to_str(feature.location)
    #         genbank_id = ""
    #         if "id" in vars(feature):
    #             genbank_id = m.Pg.clean_text_value(feature.id)
    #
    #         qualifiers = feature.qualifiers
    #         gene = ""
    #         if "gene" in qualifiers:
    #             gene = util.list_to_string(qualifiers["gene"])
    #         locus_tag = ""
    #         if "locus_tag" in qualifiers:
    #             locus_tag = util.list_to_string(qualifiers["locus_tag"])
    #         ribosomal_slippage = "FALSE"
    #         if "ribosomal_slippage" in qualifiers:
    #             ribosomal_slippage = "TRUE"
    #         gene_synonym = ""
    #         if "gene_synonym" in qualifiers:
    #             gene_synonym = util.list_to_string(qualifiers["gene_synonym"])
    #         note = ""
    #         if "note" in qualifiers:
    #             note = util.list_to_string(qualifiers["note"])
    #         codon_start = "NULL"
    #         if "codon_start" in qualifiers:
    #             codon_start = util.list_to_string(qualifiers["codon_start"])
    #         product = ""
    #         if "product" in qualifiers:
    #             product = util.list_to_string(qualifiers["product"])
    #         ncbi_protein_id = ""
    #         if "protein_id" in qualifiers:
    #             ncbi_protein_id = util.list_to_string(qualifiers["protein_id"])
    #         ncbi_db_xref = ""
    #         if "db_xref" in qualifiers:
    #             ncbi_db_xref = util.list_to_string(qualifiers["db_xref"])
    #         translation = ""
    #         if "translation" in qualifiers:
    #             translation = util.list_to_string(qualifiers["translation"])
    #
    #         pg_cur.execute(
    #             "INSERT INTO gene_cds ("
    #             "location, "
    #             "gene, "
    #             "locus_tag, "
    #             "ribosomal_slippage, "
    #             "gene_synonym, "
    #             "note, "
    #             "codon_start, "
    #             "product, "
    #             "ncbi_protein_id, "
    #             "ncbi_db_xref, "
    #             "translation, "
    #             "gene_id, "
    #             "genbank_id "
    #             ") VALUES ("
    #             f"NULLIF('{location}', ''), "
    #             f"NULLIF('{gene}', ''), "
    #             f"NULLIF('{locus_tag}', ''), "
    #             f"{ribosomal_slippage}, "
    #             f"NULLIF('{gene_synonym}', ''), "
    #             f"NULLIF('{note}', ''), "
    #             f"{codon_start}, "
    #             f"NULLIF('{product}', ''), "
    #             f"NULLIF('{ncbi_protein_id}', ''), "
    #             f"NULLIF('{ncbi_db_xref}', ''), "
    #             f"NULLIF('{translation}', ''), "
    #             f"{gene_id}, "
    #             f"NULLIF('{genbank_id}', '')"
    #             ") RETURNING id;"
    #         )
    #         gene_cds_id = pg_cur.fetchone()[0]
    #
    #         return gene_cds_id
    #
    #     else:
    #         logging.error(f"{scope_pref}: Missing code to handle this type of data.")
    #         raise Exception
    
    @classmethod
    def insert(cls,
               data=None,
               source_type=None,
               gene_id=None,
               pg_cur=None
    ):
        scope_pref = f"{cls.__name__}.insert_cds"

        # # Todo: handle lists with more than 1 item
        # ## Could use list to string utility
        # for key in data:
        #     if len(data[key]) > 1:
        #         logging.debug(f"{scope_pref}: key: {key}")
        #         logging.debug(f"{scope_pref}: value: {data[key]}")
        #         raise
        if source_type == "GenBank":
            feature = data

            if not gene_id:
                gene_id = None

            location = None
            if "location" in vars(feature):
                location = util.json_to_str(feature.location)
            genbank_id = None
            if "id" in vars(feature):
                genbank_id = m.Pg.clean_text_value(feature.id)

            qualifiers = feature.qualifiers
            gene = None
            if "gene" in qualifiers:
                gene = util.list_to_string(qualifiers["gene"])
            locus_tag = None
            if "locus_tag" in qualifiers:
                locus_tag = util.list_to_string(qualifiers["locus_tag"])
            ribosomal_slippage = False
            if "ribosomal_slippage" in qualifiers:
                ribosomal_slippage = True
            gene_synonym = None
            if "gene_synonym" in qualifiers:
                gene_synonym = util.list_to_string(qualifiers["gene_synonym"])
            note = None
            if "note" in qualifiers:
                note = util.list_to_string(qualifiers["note"])
            codon_start = None
            if "codon_start" in qualifiers:
                codon_start = util.list_to_string(qualifiers["codon_start"])
            product = None
            if "product" in qualifiers:
                product = util.list_to_string(qualifiers["product"])
            ncbi_protein_id = None
            if "protein_id" in qualifiers:
                ncbi_protein_id = util.list_to_string(qualifiers["protein_id"])
            ncbi_db_xref = None
            if "db_xref" in qualifiers:
                ncbi_db_xref = util.list_to_string(qualifiers["db_xref"])
            translation = None
            if "translation" in qualifiers:
                translation = util.list_to_string(qualifiers["translation"])

            # debug(qualifiers)
            gene_cds, created =  cls.objects.get_or_create(
                location = location,
                gene = gene,
                locus_tag = locus_tag,
                ribosomal_slippage = ribosomal_slippage,
                gene_synonym = gene_synonym,
                note = note,
                codon_start = codon_start,
                product = product,
                ncbi_protein_id = ncbi_protein_id,
                ncbi_db_xref = ncbi_db_xref,
                translation = translation,
                gene_id = gene_id,
                genbank_id = genbank_id
            )

            return gene_cds.id

        else:
            logging.error(f"{scope_pref}: Missing code to handle this type of data.")
            raise Exception


class GeneIsoform(models.Model):

    class Meta:
        db_table = "gene_isoform"

    id = models.UUIDField(
        primary_key=True,
        default=uuid.uuid4,
    )
    name = models.CharField(max_length=255)
    strand = models.IntegerField() # 1 = Reference (typically referred to as +)
                                   # 0 = Complement (typically referred to as -)
    start = models.IntegerField()
    stop = models.IntegerField()
    cds_start = models.IntegerField()
    cds_stop = models.IntegerField()
    # gene_id = models.UUIDField()
    gene = models.ForeignKey("Gene", on_delete=models.CASCADE)
    ensembl_transcript_id = models.CharField(max_length=255)
    ncbi_protein_version = models.CharField(max_length=255)

    @classmethod
    def create(cls, filepath):
        
        f = open(filepath, 'r')

        # Header row
        headers = ""
        for line in f:
            headers = line
            break
        # Data
        for line in f:
            # Convert line to list
            line = util.line_to_list(line)

            # Get the gene symbol
            actual_gene_symbol = line[9]

            # Get the genetic element ID
            genetic_element_name = line[10]
            # logging.debug(f"{scope_pref}: genetic_element_name: {genetic_element_name}")
            genetic_element_id = m.GeneticElement.objects.get(name=genetic_element_name).id

            # Get gene ID: traverse genome > m2m > genetic_element > m2m > gene tables using a view
            genetic_element_to_genes = m.GeneticElementToGene.objects.filter(
                genetic_element_id = genetic_element_id
            )
            # debug(genetic_element_to_genes)

            # genes = m.Gene.objects.filter(
            #     symbol = actual_gene_symbol
            # )

            # Get the gene
            # for genetic_element_to_gene in genetic_element_to_genes:
            #     # debug(genetic_element_to_gene.gene_id)
            #     gene = m.Gene.objects.get(
            #         symbol = actual_gene_symbol
            #     )
            #     debug(gene)

            # Get the gene ID by gene symbol
            gene_id = m.Gene.objects.get(
                symbol = actual_gene_symbol
            ).id
            # debug(gene_id)


            # Insert into gene isoform table
            cls.objects.get_or_create(
                strand = m.GenomeBase.strand_values[line[0]],
                name = line[1],
                start = line[2],
                stop = line[3],
                cds_start = line[4],
                cds_stop = line[5],
                ncbi_protein_version = line[6],
                ensembl_transcript_id = line[7],
                gene_id = gene_id
            )

        f.close()
        

class GeneIsoformV3(models.Model):

    class Meta:
        db_table = "gene_isoform_v3"

    id = models.UUIDField(
        primary_key=True,
        default=uuid.uuid4,
    )
    name = models.CharField(max_length=255)
    strand = models.IntegerField() # 1 = Reference (typically referred to as +)
                                   # 0 = Complement (typically referred to as -)
    start = models.IntegerField()
    stop = models.IntegerField()
    cds_start = models.IntegerField()
    cds_stop = models.IntegerField()
    gene_id = models.UUIDField()
    gene_name = models.CharField(max_length=255)
    actual_gene_symbol = models.CharField(max_length=255)
    ncbi_protein_version = models.CharField(max_length=255)
    ensembl_transcript_id = models.CharField(max_length=255)
    gene_id = models.UUIDField()
    ref_genome_id = models.UUIDField()

    @classmethod
    def create(cls,
        filepath = None
    ):

        if filepath:
            f = open(filepath, 'r')

            # Header row
            headers = ""
            for line in f:
                headers = line
                break
            # Data
            for line in f:
                # Convert line to list
                line = util.line_to_list(line)

                # Get reference genome ID
                ncbi_genome_version = line[10]
                genome_id = m.Genome.objects.get(ncbi_version=ncbi_genome_version).id

                # Get gene ID: traverse genome > m2m > genetic_element > m2m > gene tables using a view
                actual_gene_symbol = line[9]
                gene_id = m.ViewGenomeToGene.objects.get(
                    genome_id = genome_id,
                    symbol = actual_gene_symbol.upper()
                ).gene_id


                cls.objects.get_or_create(
                    strand = m.GenomeBase.strand_values[line[0]],
                    name = line[1],
                    start = line[2],
                    stop = line[3],
                    cds_start = line[4],
                    cds_stop = line[5],
                    gene_name=line[8].upper(),
                    actual_gene_symbol=line[9].upper(),
                    ncbi_protein_version = line[6],
                    ensembl_transcript_id = line[7],
                    gene_id = gene_id,
                    ref_genome_id = genome_id,
                )

            f.close()