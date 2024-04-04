
import logging
import os
import sys

from Bio import SeqIO
from devtools import debug
from django.db import models


import modules as m
import util


# Configs



# logging.debug(f"{scope_pref}: os.getcwd(): {os.getcwd()}")
root_dir = os.getcwd()
data_dir = os.path.join(root_dir, "data")
genomes_dir = os.path.join(data_dir, "genomes")


class GenBank(models.Model):


    @classmethod
    def load_gbk_file_into_db(cls,
        gbk_filepath = None,
        ncbi_version = None,
        genetic_element_id = None,
        genome_id = None,
        feature_qty_limit = None,
    ):
        # Current database destinations of feature types:
        ## <gbk-feature> -> <table-name>
        ## source -> no table
        ## 5'UTR -> na_utr
        ## gene -> gene
        ## CDS -> gene_cds
        ## mat_peptide -> protein
        ## stem_loop -> na_stem_loop
        ## 3'UTR -> na_utr

        scope_pref = f"{cls.__name__}.load_gbk_file_into_db"
        # logging.debug(f"{scope_pref}: os.getcwd(): {os.getcwd()}")

        ### Handle args #####################################################################################################

        gbk_file = None
        if gbk_filepath:
            gbk_file = SeqIO.read(gbk_filepath, "genbank")
        else:
            # Check if GenBank file exists, it may have a gb or gbk file extension
            gbk_filepath = os.path.join(genomes_dir, ncbi_version, ".gb")
            # logging.debug(f"{scope_pref}: gbk_filepath: {gbk_filepath}")
            if not os.path.exists(gbk_filepath):
                gbk_filepath = os.path.join(genomes_dir, ncbi_version, ".gbk")
            gbk_file = SeqIO.read(gbk_filepath, "genbank")

        ### Do load fbk file code ###########################################################################################

        # Add lineage data to database


        # Get lineage
        # util.list_members(gbk_file, "out.txt")
        # debug(gbk_file.annotations)
        # logging.debug(f"{scope_pref}: source: {gbk_file.annotations['source']}")
        # logging.debug(f"{scope_pref}: organism: {gbk_file.annotations['organism']}")
        # logging.debug(f"{scope_pref}: taxonomy: {gbk_file.annotations['taxonomy']}")
        organism_name = gbk_file.annotations['organism']
        lineage = gbk_file.annotations['taxonomy']
        # logging.debug(f"{scope_pref}: lineage: {lineage}")
        taxon = m.Taxonomy.resolve_lineage(
            name = organism_name,
            lineage = lineage
        )
        # logging.debug(f"{scope_pref}: taxon.id: {taxon.id}")
        # logging.debug(f"{scope_pref}: taxon.name: {taxon.name}")

        # Update genome table with lineage
        genome = m.Genome.objects.get(id=genome_id)
        genome.taxonomy_id = taxon.id
        genome.save()

        # Update the genetic element length
        ## Assume it's whole sequenceo of GenBank file for now
        genetic_element_length = len(gbk_file._seq)
        genetic_element = m.GeneticElement.objects.get(id=genetic_element_id)
        genetic_element.length_bp = genetic_element_length
        genetic_element.save()


        # Genbank data structure inspection
        expected_feature_types = (
            "3'UTR", "5'UTR", "CDS", "gene", "mat_peptide",
            "misc_feature", "primer_bind", "regulatory", "repeat_region",
            "sig_peptide", "source", "stem_loop"
        )
        feature_types = set()
        feature_location_attr = set()
        feature_location_operators = set()
        feature_location_parts = set()
        sig_peptide = set()
        sig_peptide_qualifiers = set()
        source_attr = set()
        source_qualifiers = set()
        utr5_attr = set()
        utr5_qualifiers = set()
        gene_attr = set()
        gene_qualifiers = set()
        cds_attr = set()
        cds_qualifiers = set()
        mat_peptide_attr = set()
        mat_peptide_qualifiers = set()
        stem_loop_attr = set()
        stem_loop_qualifiers = set()
        utr3_attr = set()
        utr3_qualifiers = set()


        # Todo: Loops through the file multiple times but should be made smarter to loop through only one time
        # Order of feature type handling seems to matter for some features
        logging.info(f"Organism: {organism_name}")
        for i, feature in enumerate(gbk_file.features):

            if feature_qty_limit:
                if i > feature_qty_limit:
                    break

            # debug(feature)

            feature_types.add(feature.type)
            # logging.debug(f"{scope_pref}: feature.type: {feature.type}")
            # logging.debug(f"{scope_pref}: feature.location: {feature.location}")
            if hasattr(feature.location, "operator"):
                # logging.debug(f"{scope_pref}: feature.location.operator: {feature.location.operator}")
                feature_location_operators.add(feature.location.operator)
                for part in feature.location.parts:
                    # logging.debug(f"{scope_pref}: feature_location_part: {part}")
                    for attr in vars(part):
                        # logging.debug(f"{scope_pref}: feature_location_part_attr: {attr}")
                        feature_location_parts.add(attr)
            elif hasattr(feature, "location"):
                for property in vars(feature.location):
                    feature_location_attr.add(property)
            # logging.debug(f"{scope_pref}: feature.qualifiers: {feature.qualifiers}")

            # The source feature type is likely always first
            if feature.type == "source":
                # logging.debug(f"{scope_pref}: feature.type: {feature.type}")
                for attr in vars(feature):
                    source_attr.add(attr)
                for qualifier in feature.qualifiers:
                    source_qualifiers.add(qualifier)

            elif feature.type == "5'UTR":
                # logging.debug(f"{scope_pref}: feature.type: {feature.type}")
                for attr in vars(feature):
                    utr5_attr.add(attr)
                for qualifier in feature.qualifiers:
                    utr5_qualifiers.add(qualifier)
                m.NaUtr.insert(
                    data=feature,
                    source_type="GenBank",
                    genetic_element_id = genetic_element_id
                )

            elif feature.type == "gene":
                gene_id = None
                # logging.debug(f"{scope_pref}: feature.type: {feature.type}")
                for attr in vars(feature):
                    gene_attr.add(attr)
                for qualifier in feature.qualifiers:
                    gene_qualifiers.add(qualifier)
                gene_id = m.Gene.insert(
                    data = feature,
                    source_type = "GenBank",
                    genetic_element_id = genetic_element_id
                )

            # Likely always follows a gene feature type
            elif feature.type == "CDS":
                gene_cds_id = None
                # logging.debug(f"{scope_pref}: feature.type: {feature.type}")
                for attr in vars(feature):
                    cds_attr.add(attr)
                for qualifier in feature.qualifiers:
                    cds_qualifiers.add(qualifier)
                gene_cds_id = m.GeneCds.insert(
                    data = feature,
                    source_type = "GenBank",
                    gene_id = gene_id
                )
                m.Protein.insert(
                    data=feature,
                    source_type="GenBank",
                    gene_cds_id=gene_cds_id
                )

            # Likely always follows a gene feature type
            elif feature.type == "mat_peptide":
                # logging.debug(f"{scope_pref}: feature.type: {feature.type}")
                for attr in vars(feature):
                    mat_peptide_attr.add(attr)
                for qualifier in feature.qualifiers:
                    mat_peptide_qualifiers.add(qualifier)
                m.Protein.insert(
                    data = feature,
                    source_type = "GenBank",
                    gene_cds_id = gene_cds_id
                )

            # Todo: The for loop causes stem loops not within genes to be assigned one anyway.
            #       This needs to be fixed. Try mapping gene loci and comparing to that.
            # Can be outside of a gene
            elif feature.type == "stem_loop":
                # logging.debug(f"{scope_pref}: feature.type: {feature.type}")
                for attr in vars(feature):
                    stem_loop_attr.add(attr)
                for qualifier in feature.qualifiers:
                    stem_loop_qualifiers.add(qualifier)
                m.NaStemLoop.insert(
                    data = feature,
                    source_type = "GenBank",
                    gene_id = gene_id,
                    genetic_element_id = genetic_element_id
                )

            elif feature.type == "3'UTR":
                # logging.debug(f"{scope_pref}: feature.type: {feature.type}")
                for attr in vars(feature):
                    utr3_attr.add(attr)
                for qualifier in feature.qualifiers:
                    utr3_qualifiers.add(qualifier)
                m.NaUtr.insert(
                    data = feature,
                    source_type = "GenBank",
                    genetic_element_id = genetic_element_id
                )

            elif feature.type == "sig_peptide":
                # logging.debug(f"{scope_pref}: feature.type: {feature.type}")
                for attr in vars(feature):
                    sig_peptide.add(attr)
                for qualifier in feature.qualifiers:
                    sig_peptide_qualifiers.add(qualifier)
                m.SignalingPeptide.insert(
                    data = feature,
                    source_type = "GenBank",
                    genetic_element_id = genetic_element_id
                )

            elif feature.type == "regulatory":
                # logging.debug(f"{scope_pref}: feature.type: {feature.type}")
                for attr in vars(feature):
                    sig_peptide.add(attr)
                for qualifier in feature.qualifiers:
                    sig_peptide_qualifiers.add(qualifier)
                m.RegulatoryNaSeq.insert(
                    data=feature,
                    source_type="GenBank",
                    genetic_element_id=genetic_element_id
                )
                
            sys.stdout.write(f"\r  Processed {i} features.")

        sys.stdout.write("\n")

        # # Genbank data structure inspection
        # logging.debug(f"{scope_pref}: feature_types: {feature_types}")
        # logging.debug(f"{scope_pref}: feature_location_attr: {feature_location_attr}")
        # logging.debug(f"{scope_pref}: feature_location_operators: {feature_location_operators}")
        # logging.debug(f"{scope_pref}: feature_location_parts: {feature_location_parts}")
        # logging.debug(f"{scope_pref}: source_attr: {source_attr}")
        # logging.debug(f"{scope_pref}: source_qualifiers: {source_qualifiers}")
        # logging.debug(f"{scope_pref}: utr5_attr: {utr5_attr}")
        # logging.debug(f"{scope_pref}: utr5_qualifiers: {utr5_qualifiers}")
        # logging.debug(f"{scope_pref}: gene_attr: {gene_attr}")
        # logging.debug(f"{scope_pref}: gene_qualifiers: {gene_qualifiers}")
        # logging.debug(f"{scope_pref}: cds_attr: {cds_attr}")
        # logging.debug(f"{scope_pref}: cds_qualifiers: {cds_qualifiers}")
        # logging.debug(f"{scope_pref}: mat_peptide_attr: {mat_peptide_attr}")
        # logging.debug(f"{scope_pref}: mat_peptide_qualifiers: {mat_peptide_qualifiers}")
        # logging.debug(f"{scope_pref}: stem_loop_attr: {stem_loop_attr}")
        # logging.debug(f"{scope_pref}: stem_loop_qualifiers: {stem_loop_qualifiers}")
        # logging.debug(f"{scope_pref}: utr3_attr: {utr3_attr}")
        # logging.debug(f"{scope_pref}: utr3_qualifiers: {utr3_qualifiers}")
        # Throw error if feature type wasn't handled by the above code
        for feature_type in feature_types:
            if feature_type not in expected_feature_types:
                logging.error(f"{scope_pref}: feature type '{feature_type}' was unexpected and is not handled. Write handling code.")
                raise
