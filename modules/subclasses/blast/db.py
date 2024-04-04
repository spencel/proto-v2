
import os
import pathlib
import re
import subprocess

import config
import modules as m
import util


class Db():


    @staticmethod
    def get_dpath(fpath_prefix) -> str:
        """
        Gets the directory path of the BLAST database from the filename prefix path.
        :param fpath_prefix:
        :return:
        """
        # >>> existGDBPath = pathlib.Path(r'T:\Data\DBDesign\DBDesign_93_v141b.mdb')
        # >>> wkspFldr = existGDBPath.parent
        return str(pathlib.Path(fpath_prefix).parent)


    @staticmethod
    def make(
        in_fpath,
        dbtype=None,
        input_type=None,
        out_fpath_prefix=None,
        genomes=None,
        taxon_id_map_fpath=None,
        overwrite=False,
        remove_input_files=False
    ):
        """
        Argumenents for makeblastdb command:  https://www.ncbi.nlm.nih.gov/books/NBK279684/table/appendices.T.makeblastdb_application_opt/
        :param in_fpath:
        :param dbtype:
        :param input_type:
        :param out_fpath_prefix:
        :param genomes:
        :param taxon_id_map_fpath:
        :param overwrite:
        :param remove_input_files:
        :return:
        """

        if out_fpath_prefix:
            blastdb_dpath = os.path.split(out_fpath_prefix)[0]
        else:
            blastdb_dpath = m.File.get_fpath_without_extension(in_fpath)
            out_fpath_prefix = os.path.join(blastdb_dpath, "_")
        print(f"sl: out_fpath_prefix: {out_fpath_prefix}")


        # Make the BLAST database if it doesn't exist yet
        if m.Dir.make_if_not_exist(blastdb_dpath) == False:
            if overwrite:
                m.Dir.delete(blastdb_dpath)
            # Exit if the database already exists
            else:
                return

        if genomes:
            # Concatenate genome fasta files into a single file
            fasta_fpaths = []
            for genome in genomes:
                fasta_fpaths.append(genomes[genome]["fasta_fpath"])
            m.Fasta.combine(fasta_fpaths, in_fpath)

        if not dbtype:
            dbtype = "nucl"
        if not input_type:
            input_type = "fasta"

        print(f"sl: in_fpath: {in_fpath}")
        print(f"sl: is_file: {os.path.isfile(in_fpath)}")
        command = [
            f"{config.blast.ncbi_blast_bin_dpath}/makeblastdb",
            "-in", f"{in_fpath}",
            "-dbtype", dbtype,
            "-input_type", input_type,
            "-out", f"{out_fpath_prefix}"
        ]
        print(f"sl: command: {command}")
        if taxon_id_map_fpath:
            command.extend(["-taxid_map", taxon_id_map_fpath])

        subprocess.run(command)

        # Delete the combined FASTA file
        if remove_input_files:
            os.remove(in_fpath)


    @classmethod
    def count_species(cls, db_fpath_prefix):

        path_split = os.path.split(db_fpath_prefix)
        db_dpath = path_split[0]
        entries_fpath = os.path.join(db_dpath, "entries.tsv")
        if not os.path.isfile(entries_fpath):
            cls.get(db_fpath_prefix)

        with open(entries_fpath, 'r') as f:
            taxon_ids = set()
            for line in f:
                taxon_ids.add(line.split("\t")[7])
            print(f"sl: len(taxon_ids): {len(taxon_ids)}")


    @classmethod
    def get(cls,
        source=None,
        database=None,
        database_fname_prefix=None,
        download_destination=None,
        update=True,
        redownload=False
    ):

        if not source:
            source = cls.DEFAULT_DATABASE_SOURCE
        if not database:
            databases = cls.DEFAULT_DATABASES
        if not database_fname_prefix:
            database_fname_prefix = cls.DEFAULT_DATABASES_FILENAME_PREFIX

        # Save original working directory path
        original_working_dpath = os.getcwd()

        this_blast_db_dpath: str
        if not download_destination:
            # Get dpath of this BLAST database
            this_blast_db_dpath = os.path.join(
                config.paths.blastdb,
                database
            )
        else:
            this_blast_db_dpath = os.path.join(
                download_destination,
                database
            )

        # Force download even if BLAST database is up-to-date
        if redownload:
            util.delete_directory(this_blast_db_dpath)
        # Don't update BLAST database, but download it if it hasn't been already
        # Todo: for now it just exits before downloading
        elif update:
            return

        util.create_dpath_if_not_exist(this_blast_db_dpath)

        # Change directory to the BLAST database
        os.chdir(this_blast_db_dpath)

        update_blastdb_fpath = os.path.join(
            cls.ncbi_blast_bin_dpath,
            "update_blastdb.pl"
        )

        subprocess.run([
            "perl", update_blastdb_fpath,
            "--source", source,
            "--decompress", database
        ])

        # Rename files to the prefix
        # Todo: do this later

        # Change directory back to what it was originally
        os.chdir(original_working_dpath)


    @classmethod
    def export_entries(
        cls,
        db_fpath_prefix,
        out_fpath=None,
        overwrite=False,
        is_custom_defline=False
    ):

        #
        if not out_fpath:
            db_dpath = cls.get_dpath(db_fpath_prefix)
            out_fpath = os.path.join(
                os.getcwd(),
                db_dpath + "-entries.tsv"
            )

        # This command doesn't work unless executed from the same directory as the BLAST database
        path_split = os.path.split(db_fpath_prefix)
        # Get BLAST DB path
        db_dpath = path_split[0]
        # Get BLAST filename prefix
        db_fname_prefix = path_split[1]
        # Save original working directory path
        original_working_dpath = os.getcwd()
        # Change directory to the BLAST database
        os.chdir(db_dpath)

        if overwrite:
            m.File.delete(out_fpath)

        # The following was taken from the blastdbcmd -help output
        # %a = 0, accession (eg, NC_002645.1)
        # %g = 1, gi (eg, 12175745)
        # %o = 2, ordinal id (OID) (eg, 113)
        # %i = 3, sequence id (eg, ref|NC_002645.1|)
        # %t = 4, sequence title (eg, Human coronavirus 229E, complete genome)
        #         This can also be where custom entry data is stored from a FASTA defline.
        #         If the FASTA's defline was tab delimited, then the sequence title will be 3-space delimited: "   "
        # %l = 5, sequence length (eg, 27317)
        # %h = 6, sequence hash value (eg, 0X14677E4E)
        # %T = 7, taxid (eg, 11137)
        # %X = 8, leaf-node taxids (eg, 11137)
        # %e = 9, membership integer (eg, 4)
        # %L = 10, common taxonomic name (Human coronavirus 229E)
        # %C = 11, common taxonomic names for leaf-node taxids (eg, Human coronavirus 229E)
        # %S = 12, scientific name (eg, Human coronavirus 229E)
        # %N = 13, scientific names for leaf-node taxids (eg, Human coronavirus 229E)
        # %B = 14, BLAST name (eg, viruses)
        # %K = 15, taxonomic super kingdom (eg, Viruses)
        # %P = 16, PIG (eg, -1)
        # Don't update entries.tsv if it already
        outfmt: str
        if is_custom_defline:
            outfmt = "%t"
        else:
            outfmt = "%a\t%g\t%o\t%i\t%t\t%l\t%h\t%T\t%X\t%e\t%L\t%C\t%S\t%N\t%B\t%K\t%P"

        if not os.path.exists(out_fpath):
            subprocess.run([
                f"{config.blast.ncbi_blast_bin_dpath}/blastdbcmd",
                "-db", db_fname_prefix,
                "-entry", "all",
                "-outfmt", outfmt,
                "-out", out_fpath
            ])
        # Change directory back to what it was originally
        os.chdir(original_working_dpath)
        # print(f"get_entries: end")