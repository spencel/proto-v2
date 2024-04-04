
import logging
import os
import uuid

from devtools import debug
from django.db import models


import modules as m
import util


# Configs





class NaUtr(models.Model):

    class Meta:
        db_table = "na_utr"

    id = models.UUIDField(
        primary_key=True,
        default=uuid.uuid4,
    )
    type = models.CharField(max_length=255)
    location = models.TextField()
    genetic_element_id = models.UUIDField()
    genbank_id = models.CharField(max_length=255)


    @classmethod
    def insert(cls,
        data = None,
        source_type = None,
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

            type = m.Pg.clean_text_value(feature.type)

            location = None
            if "location" in vars(feature):
                location = util.json_to_str(feature.location)
            genbank_id = None
            if "id" in vars(feature):
                genbank_id = m.Pg.clean_text_value(feature.id)

            if not genetic_element_id:
                genetic_element_id = None

            # pg_cur.execute(
            #     "INSERT INTO na_utr ("
            #         "type, "
            #         "location, "
            #         "genetic_element_id, "
            #         "genbank_id"
            #     ") VALUES ("
            #         f"NULLIF('{type}', ''), "
            #         f"NULLIF('{location}', ''), "
            #         f"{genetic_element_id}, "
            #         f"NULLIF('{genbank_id}', '')"
            #     ") RETURNING id;"
            # )
            cls.objects.get_or_create(
                type = type,
                location = location,
                genetic_element_id = genetic_element_id,
                genbank_id = genbank_id
            )

        else:
            logging.error(f"{scope_pref}: Missing code to handle this type of data.")
            raise Exception