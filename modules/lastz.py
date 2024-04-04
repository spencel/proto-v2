
import json
import logging
import lzma
import os
import subprocess
import uuid

from devtools import debug
from django.db import models

import modules as m



# Configs


root_dir = os.getcwd()
data_dir = os.path.join(root_dir, "data")
genomes_dir = os.path.join(data_dir, "genomes")
modules_dir = os.path.join(root_dir, "modules")
modules_data_dir = os.path.join(modules_dir, "data")
lastz_data_dir = os.path.join(modules_data_dir, "lastz")


class Lastz(models.Model):

    class Meta:
        db_table = "lastz"

    id = models.UUIDField(
        primary_key=True,
        default=uuid.uuid4,
    )
    status = models.IntegerField()
    status_values = {
        "NOT_STARTED": 0,
        "DONE": 1,
        "STARTING": 2
    }
    organism_sample_id = models.UUIDField()
    reference_genome_id = models.UUIDField()

    @classmethod
    def create_qry_file(cls, filepath, fa_defline, sequence):
        with open(filepath, 'w') as f:
            f.write(f'>{fa_defline}\n')
            f.write(sequence)

    @classmethod
    def run(cls, ref_fasta_filepath, query_filepath, output_filepath):
        subprocess.run([
            'lastz',
            ref_fasta_filepath,
            query_filepath,
            '--ambiguous=iupac',
            '--noytrim',
            '--format=differences',
            f"--output={output_filepath}"
        ])
        return

    @classmethod
    def process_sample_gisaid(cls,
        organism_sample_gisaid,
        ref_genome
    ):

        # Get GISAID provision filepath
        file_system = m.Os.objects.get(
            id = organism_sample_gisaid.provision_file_id
        )
        provision_filepath = file_system.filepath

        # Define lastz query filepath
        qry_filepath = os.path.join(lastz_data_dir, f"{organism_sample_gisaid.id}-query.fa")

        # Get the sequence from the provision json.xz archive
        with lzma.open(provision_filepath, 'r') as f:

            # Go to the location of this sample
            f.seek(organism_sample_gisaid.file_position)

            # Shortcut for decompressing a json line
            for line in f:

                # Convert line to json and create query file
                cls.create_qry_file(
                    filepath = qry_filepath,
                    fa_defline = organism_sample_gisaid.id,
                    sequence = json.loads(line)["sequence"]
                )

                # Only reading one line
                break

        # Declare reference and output filepaths
        ref_genome_dirname = ref_genome.dir_name
        ref_fasta_filepath = os.path.join(genomes_dir, ref_genome_dirname, ".fa")
        output_filepath = os.path.join(lastz_data_dir, f"{organism_sample_gisaid.id}.out")

        # Run LASTZ
        cls.run(
            ref_fasta_filepath,
            qry_filepath,
            output_filepath
        )

        # Upsert the variant segments
        m.VariantSegment.add_from_lastz_output(
            lastz_output_filepath = output_filepath,
            organism_sample_id = organism_sample_gisaid.organism_sample_id,
            reference_genome_id = ref_genome.id
        )

        # Remove the files from the data folder
        os.remove(qry_filepath)
        os.remove(output_filepath)

        return




    @classmethod
    def process_sample(cls,
        lastz_row
    ):

        # Get organism genome sample source and id
        view_organism_sample = m.ViewOrganismSample.objects.get(
            sample_id=lastz_row.organism_sample_id
        )

        if view_organism_sample:
            # debug(view_organism_sample)
            # logging.debug(f"{scope_pref}: organism_sample_id: {view_organism_sample.sample_id}")
            # logging.debug(f"{scope_pref}: source_name: {view_organism_sample.source_name}")

            # Get the reference genome
            ref_genome = m.Genome.objects.get(
                id = lastz_row.reference_genome_id
            )

            if view_organism_sample.source_name == "GISAID":

                organism_sample_gisaid = m.OrganismSampleGisaid.objects.get(
                    organism_sample_id = view_organism_sample.sample_id
                )

                if organism_sample_gisaid:

                    cls.process_sample_gisaid(
                        organism_sample_gisaid,
                        ref_genome
                    )

                else:

                    raise Exception

                # Switch the sample record in the lastz table to DONE
                lastz_row.status = cls.status_values["DONE"]
                lastz_row.save()

        else:

            raise Exception

    @classmethod
    def get_absent_organism_sample(cls):

        return m.ViewOrganismSample.objects.raw(
            "SELECT "
                "id, "
                "sample_id, "
                "na_sequence_hash, "
                "source_name "
            "FROM view_organism_sample "
            "WHERE NOT EXISTS ("
                "SELECT organism_sample_id "
                "FROM lastz "
                "WHERE view_organism_sample.sample_id = lastz.organism_sample_id"
            ") LIMIT 1; "
        )
    
    @classmethod
    def add_organism_sample(cls, reference_genome_id):

        res = cls.get_absent_organism_sample()

        if len(res) > 0:

            view_organism_sample = res[0]

            cls.objects.get_or_create(
                status = cls.status_values["NOT_STARTED"],
                organism_sample_id = view_organism_sample.sample_id,
                reference_genome_id = reference_genome_id
            )

            return True

        # There are no organism samples to add to the lastz table
        else:

            return False
    
    @classmethod
    def start_processing_a_sample(cls,
        reference_genome_id
    ):

        # Check for rows in lastz table that haven't started yet and set their status to started
        res = cls.objects.raw(
            "UPDATE lastz SET "
                f"status = {cls.status_values['STARTING']} "
            "WHERE id = ("
                "SELECT id "
                "FROM lastz "
                "WHERE "
                    "status = 0 "
                    f"AND reference_genome_id = '{reference_genome_id}' "
                "LIMIT 1"
            ")  "
            "RETURNING *;"
        )
        if len(res) > 0:
            lastz_row = res[0]
            cls.process_sample(lastz_row=lastz_row)
            return True

        else:
            return False
        







    

        