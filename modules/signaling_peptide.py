
import logging
import os
import uuid

from devtools import debug
from django.db import models


import modules as m
import util


# Configs





class SignalingPeptide(models.Model):

    class Meta:
        db_table = "signaling_peptide"

    id = models.UUIDField(
        primary_key=True,
        default=uuid.uuid4,
    )
    location = models.TextField()
    gene_symbol = models.CharField(255)
    genetic_element_id = models.UUIDField()

    @classmethod
    def insert(cls,
        data = None,
        source_type = None,
        genetic_element_id = None
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
            gene_symbol = None
            if "gene" in vars(feature):
                gene_symbol = m.Pg.clean_text_value(feature.gene)


            if not genetic_element_id:
                genetic_element_id = None

            cls.objects.get_or_create(
                location = location,
                gene_symbol = gene_symbol,
                genetic_element_id = genetic_element_id
            )

        else:
            logging.error(f"{scope_pref}: Missing code to handle this type of data.")
            raise Exception