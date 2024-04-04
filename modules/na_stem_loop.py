
import logging
import os
import uuid

from devtools import debug
from django.db import models


import modules as m
import util


# Configs





class NaStemLoop(models.Model):

    class Meta:
        db_table = "na_stem_loop"

    id = models.UUIDField(
        primary_key=True,
        default=uuid.uuid4,
    )
    location = models.TextField()
    gene = models.CharField(max_length=255)
    locus_tag = models.CharField(max_length=255)
    inference = models.TextField()
    note = models.TextField()
    function = models.TextField()
    gene_id = models.UUIDField()
    genetic_element_id = models.UUIDField()
    genbank_id = models.CharField(max_length=255)


    # @classmethod
    # def insert_old(cls,
    #     data = None,
    #     source_type = None,
    #     gene_id = None,
    #     genetic_element_id = None,
    #     pg_cur = None
    # ):
    #     scope_pref = f"{cls.__name__}.insert"
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
    #             gene = m.Pg.clean_text_value(qualifiers["gene"])
    #         locus_tag = ""
    #         if "locus_tag" in qualifiers:
    #             locus_tag = m.Pg.clean_text_value(qualifiers["locus_tag"])
    #         inference = ""
    #         if "inference" in qualifiers:
    #             inference = m.Pg.clean_text_value(qualifiers["inference"])
    #         note = ""
    #         if "note" in qualifiers:
    #             note = m.Pg.clean_text_value(qualifiers["note"])
    #         function = ""
    #         if "function" in qualifiers:
    #             function = m.Pg.clean_text_value(qualifiers["function"])
    #         if not gene_id:
    #             gene_id = "NULL"
    #         if not genetic_element_id:
    #             genetic_element_id = "NULL"
    #         pg_cur.execute(
    #             "INSERT INTO na_stem_loop ("
    #                 "location, "
    #                 "gene, "
    #                 "locus_tag, "
    #                 "inference, "
    #                 "note, "
    #                 "function, "
    #                 "gene_id, "
    #                 "genetic_element_id, "
    #                 "genbank_id"
    #             ") VALUES ("
    #                 f"NULLIF('{location}', ''), "
    #                 f"NULLIF('{gene}', ''), "
    #                 f"NULLIF('{locus_tag}', ''), "
    #                 f"NULLIF('{inference}', ''), "
    #                 f"NULLIF('{note}', ''), "
    #                 f"NULLIF('{function}', ''), "
    #                 f"{gene_id}, "
    #                 f"{genetic_element_id}, "
    #                 f"NULLIF('{genbank_id}', '')"
    #             ") RETURNING id;"
    #         )
    #
    #     else:
    #         logging.error(f"{scope_pref}: Missing code to handle this type of data.")
    #         raise Exception
    
    @classmethod
    def insert(cls,
        data = None,
        source_type = None,
        gene_id = None,
        genetic_element_id = None,
        pg_cur = None
    ):
        scope_pref = f"{cls.__name__}.insert"

        # # Todo: handle lists with more than 1 item
        # ## Could use list to string utility
        # for key in data:
        #     if len(data[key]) > 1:
        #         logging.debug(f"{scope_pref}: key: {key}")
        #         logging.debug(f"{scope_pref}: value: {data[key]}")
        #         raise
        if source_type == "GenBank":
            feature = data

            location = None
            if "location" in vars(feature):
                location = util.json_to_str(feature.location)
            genbank_id = None
            if "id" in vars(feature):
                genbank_id = m.Pg.clean_text_value(feature.id)

            qualifiers = feature.qualifiers
            gene = None
            if "gene" in qualifiers:
                gene = m.Pg.clean_text_value(qualifiers["gene"])
            locus_tag = None
            if "locus_tag" in qualifiers:
                locus_tag = m.Pg.clean_text_value(qualifiers["locus_tag"])
            inference = None
            if "inference" in qualifiers:
                inference = m.Pg.clean_text_value(qualifiers["inference"])
            note = None
            if "note" in qualifiers:
                note = m.Pg.clean_text_value(qualifiers["note"])
            function = None
            if "function" in qualifiers:
                function = m.Pg.clean_text_value(qualifiers["function"])
            if not gene_id:
                gene_id = None
            if not genetic_element_id:
                genetic_element_id = None

            cls.objects.get_or_create(
                location = location,
                gene = gene,
                locus_tag = locus_tag,
                inference = inference,
                note = note,
                function = function,
                gene_id = gene_id,
                genetic_element_id = genetic_element_id,
                genbank_id = genbank_id
            )

        else:
            logging.error(f"{scope_pref}: Missing code to handle this type of data.")
            raise Exception