
import logging
import os
import subprocess
import pathlib
import uuid

from devtools import debug
from django.db import models

import config
import modules as m
import util

# Subclasses
from modules.subclasses.blast.db import Db as DbSubclass

# Configs



class Blast():

    # BLAST modes
    DEFAULT = "DEFAULT"
    DEFAULT_DATABASES = [
        "ref_euk_rep_genomes",
        "ref_prok_rep_genomes",
        "ref_viruses_rep_genomes"
    ]
    DEFAULT_DATABASE_SOURCE = "ncbi"
    CROSS_REACTIVITY_ANALYSIS = "CROSS_REACTIVITY_ANALYSIS"
    DEFAULT_DATABASES_FILENAME_PREFIX = config.blast.blastdb_prefix

    ncbi_blast_bin_dpath = config.blast.ncbi_blast_bin_dpath

    Db = DbSubclass

    def __init__(self, data_dpath, blastdb_fpath_prefixes=None):
        """
        :param data_dpath:
        :param blastdb_fpath_prefixes: a path or list of paths
        """
        self.data_dpath = data_dpath
        self.in_fpath = os.path.join(data_dpath, "blastn.in")
        self.blastdb_fpath_prefixes = blastdb_fpath_prefixes




    @classmethod
    def run_blastn_raw(cls,
        blastdb_fpath_prefix,
        in_fpath,
        out_fpath,
        mode=DEFAULT,
        use_debug=False
    ):
        # print(f"blastdb_fpath_prefix: {blastdb_fpath_prefix}")
        # Initialize command
        command = [
            "blastn",
            "-db", blastdb_fpath_prefix,
            "-query", in_fpath,  # change to not read from file, since it's not an R script?
        ]
        # Print to console or file
        if not use_debug:
            command.extend(["-out", out_fpath])
        # Common arguments
        command.extend([
            "-penalty", str(-1),
            "-task", "blastn-short"
        ])
        # Modes
        # Cross-reactivity Analysis
        # Needs defline (sseqid) and the number of identities (nident)
        if mode == cls.CROSS_REACTIVITY_ANALYSIS:
            command.extend([
                "-evalue", str(1_000_000_000),
                "-outfmt", "6 sseqid nident qseq sseq mismatch gapopen gaps"
            ])
        # Default mode
        # Needs alignment start position (sstart), e-value (evalue), and percent identities (pident)
        else:
            command.extend([
                "-evalue", str(0.1),
                "-outfmt", "6 sstart evalue pident"
            ])
        # print(f"command: {command}")

        # Run the command
        subprocess.run(command)


    def run_blastn(self,
        blastdb_fpath_prefixes=None,
        qry_seq=None,
        mode=DEFAULT,
        use_debug=False
    ):
        """

        :param blastdb_fpath_prefixes: a list
        :param qry_seq:
        :param mode:
        :param use_debug:
        :return:
        """
        # The BLAST database path may not be specified if for example it's iterating over multiple genomes
        if not blastdb_fpath_prefixes:
            blastdb_fpath_prefixes = self.blastdb_fpath_prefixes

        # Write the sequence to the input file if it hasn't yet
        if qry_seq:
            with open(self.in_fpath, 'w') as f:
                f.write(qry_seq)

        if blastdb_fpath_prefixes:
            print(f"sl: blastdb_fpath_prefixes: {blastdb_fpath_prefixes}")
            for blastdb_fpath_prefix in blastdb_fpath_prefixes:
                blastdb_name = blastdb_dname = util.split_path(blastdb_fpath_prefix)[-2]
                out_fpath = os.path.join(
                    self.data_dpath,
                    blastdb_name + "-blastn.out"
                )
                logging.info(f"BLAST DB: {blastdb_name}")
                self.run_blastn_raw(blastdb_fpath_prefix, self.in_fpath, out_fpath, mode, use_debug)