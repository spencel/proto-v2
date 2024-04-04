

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


class Taxonomy(models.Model):

    class Meta:
        db_table = "taxonomy"

    id = models.UUIDField(
        primary_key=True,
        default=uuid.uuid4,
    )
    name = models.CharField(max_length=255)
    rank = models.CharField(max_length=255)
    common_name = models.CharField(max_length=255)
    parent_id = models.UUIDField()

    @classmethod
    def dump_taxon_tree(cls,
        tree_filepath = None,
        indent_size = 4,
    ):

        tree_file = open(tree_filepath, 'w')

        # First select all roots, ie, they have no parent ID
        root_taxons = cls.objects.filter(parent_id=None)

        # Recursive function
        def dump_taxon_tree_recurse(children, level=0):

            # Increment level
            this_level = level + 1

            # Loop through children
            for child in children:
                indent_str = ("|" + " "*(indent_size-1)) * this_level
                tree_file.write(f"{indent_str}{child.name}\n")

                sub_children = cls.objects.filter(parent_id=child.id)
                if sub_children:
                    dump_taxon_tree_recurse(sub_children, this_level)

        # Initialize recursion
        for root_taxon in root_taxons:
            # logging.debug(f"{scope_pref}: root_taxon: {root_taxon}")

            tree_file.write(f"{root_taxon.name}\n")

            # Get children
            children = cls.objects.filter(parent_id=root_taxon.id)

            if children:
                dump_taxon_tree_recurse(children)

        tree_file.close()



    # This function gets the taxonomy ID or creates the lineage
    @classmethod
    def resolve_lineage(cls,
        name = None,
        common_name = None,
        lineage = None
    ):
        
        # Flesh out lineage in case it doesn't exist
        prev_taxon = None
        for lineage_name in lineage:

            # logging.debug(f"{scope_pref}: lineage_name: {lineage_name}")

            parent_id = None
            if prev_taxon:
                parent_id = prev_taxon.id
            # logging.debug(f"{scope_pref}: parent_id: {parent_id}")


            this_taxon, created = cls.objects.get_or_create(
                parent_id = parent_id,
                name = lineage_name
            )

            prev_taxon = this_taxon

        parent_id = prev_taxon.id
        # logging.debug(f"{scope_pref}: name: {name}")
        # logging.debug(f"{scope_pref}: common_name: {common_name}")
        # logging.debug(f"{scope_pref}: parent_id: {parent_id}")
        # Add this taxon
        this_taxon, created = cls.objects.get_or_create(
            parent_id = parent_id,
            name = name,
            defaults = {
                "common_name": common_name
            }
        )

        return this_taxon
        
        
    @classmethod
    def load_data_from_file(cls,
        filepath = None
    ):

        with open(filepath, 'r') as f:
            # Skip header row
            for line in f:
                headers = line
                break
            # Data
            col_idx = {
                "name": 0,
                "common_name": 1,
                "lineage": 2
            }
            for line in f:
                # Convert line to list
                line = line.rstrip("\n").split("\t")
                name = line[col_idx["name"]]
                common_name = line[col_idx["common_name"]]
                lineage = line[col_idx["lineage"]].split("|")

                cls.resolve_lineage(
                    name = name,
                    common_name = common_name,
                    lineage = lineage
                )

                