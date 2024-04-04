
import logging
import os
import uuid

from devtools import debug
from django.db import models

import modules as m
import util


# Configs




class Exon(models.Model):

    id = models.UUIDField(
        primary_key=True,
        default=uuid.uuid4,
    )
    start = models.IntegerField()
    stop = models.IntegerField()
    name = models.CharField(max_length=255)
    gene_isoform = models.ForeignKey("GeneIsoform", on_delete=models.PROTECT)

    class Meta:
        db_table = "exon"

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

            start = line[0]
            stop = line[1]
            name = line[2]
            gene_isoform_name = line[3]

            # Get gene isoform ID
            gene_isoform_id = m.GeneIsoform.objects.get(
                name = gene_isoform_name
            ).id

            # Insert into exon table
            cls.objects.get_or_create(
                start = start,
                stop = stop,
                name = name,
                gene_isoform_id = gene_isoform_id
            )

        f.close()
