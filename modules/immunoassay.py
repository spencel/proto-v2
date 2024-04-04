
import logging
import os
import uuid

from devtools import debug
from django.db import models


import modules as m


# Configs




class Immunoassay(models.Model):

    id = models.UUIDField(
        primary_key=True,
        default=uuid.uuid4,
    )
    target_aa_start = models.IntegerField()
    target_aa_end = models.IntegerField()
    reference_genome_id = models.UUIDField()
    target_gene_id = models.UUIDField()

    class Meta:
        db_table = "immunoassay"
        constraints = [
            models.UniqueConstraint(
                name="columns",
                fields=[
                    "target_aa_start",
                    "target_aa_end",
                    "reference_genome_id", # Todo: remove? since the table should be specific to a reference genome anyway
                    "target_gene_id"
                ]
            )
        ]

    @classmethod
    def insert(cls, data, pg_cur):
        scope_pref = f"{cls.__name__}.insert"
        # logging.debug(f"{scope_pref}: data: {data}")

        # Handle list/row data

        if isinstance(data, list):
            target_aa_start = data[0]
            target_aa_end = data[1]
            reference_genome_id = data[2]
            target_gene_id = data[3]
            # logging.debug(f"{scope_pref}: data: {data}")
        else:
            target_aa_start = data["target_aa_start"]
            target_aa_end = data["target_aa_end"]
            reference_genome_id = data["reference_genome_id"]
            target_gene_id = data["target_gene_id"]

            # # Handle blank number values
            # if data[4] == "":
            #     data[4] = "NULL"

        # Return id if it exists
        pg_cur.execute(
            "SELECT id FROM immunoassay WHERE "
                f"target_aa_start = {target_aa_start} "
                f"AND target_aa_end = {target_aa_end} "
                f"AND target_gene_id = {target_gene_id};"
        )
        res = pg_cur.fetchone()
        # logging.debug(f"{scope_pref}: res: {res}")

        if res:
            return res[0]
        else:
            pg_cur.execute(
                f"INSERT INTO immunoassay ("
                    "target_aa_start, "
                    "target_aa_end, "
                    "reference_genome_id, "
                    "target_gene_id"
                ") VALUES ("
                    f"{target_aa_start}, "
                    f"{target_aa_end}, "
                    f"{reference_genome_id}, "
                    f"{target_gene_id}"
                ") RETURNING id;"
            )
            return pg_cur.fetchone()[0]


class ImmunoassayGroup(models.Model):

    id = models.UUIDField(
        primary_key=True,
        default=uuid.uuid4,
    )
    immunoassay_ids = models.TextField()

    class Meta:
        db_table = "immunoassay_group"
        constraints = [
            models.UniqueConstraint(name="columns", fields=["immunoassay_ids"])
        ]


class ImmunoassayToGroup(models.Model):

    id = models.UUIDField(
        primary_key=True,
        default=uuid.uuid4,
    )
    immunoassay_id = models.UUIDField()
    immunoassay_group_id = models.UUIDField()

    class Meta:
        db_table = "immunoassay_to_group"
        constraints = [
            models.UniqueConstraint(name="columns", fields=["immunoassay_id", "immunoassay_group_id"])
        ]


class ImmunoassayIncident(models.Model):

    id = models.UUIDField(
        primary_key=True,
        default=uuid.uuid4,
    )
    immunoassay_id = models.UUIDField()
    variant_segment_id = models.UUIDField()

    class Meta:
        db_table = "immunoassay_incident"
        constraints = [
            models.UniqueConstraint(name="columns", fields=["immunoassay_id", "variant_segment_id"])
        ]

# Immunoassay Incident Severity
statuses = models.IntegerChoices(
    "statuses",
    "SCORING_NOT_STARTED " 
    "SCORING_COMPLETE "
    "SCORING_STARTED"
)
# statuses.SCORING_NOT_STARTED = 1
# status.values = [1, 2, 3]
# status.choices = [
# 	(1, 'Scoring Not Started',),
#   (2, 'Scoring Complete',),
# 	(3, 'Scoring Started',),
# ]
class ImmunoassayIncidentSeverity(models.Model):

    id = models.UUIDField(
        primary_key=True,
        default=uuid.uuid4,
    )
    statuses = statuses
    status = models.IntegerField(choices=statuses.choices)
    severity = models.IntegerField()
    full_severity = models.IntegerField()
    variant_codes = models.TextField()
    immunoassay_id = models.UUIDField()
    organism_sample_id = models.UUIDField()

    class Meta:
        db_table = "immunoassay_incident_severity"
        constraints = [
            models.UniqueConstraint(name="columns", fields=["immunoassay_id","organism_sample_id"])
        ]


