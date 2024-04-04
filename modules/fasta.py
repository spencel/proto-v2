
import fileinput
import json
import logging
import os
import re
import sys

from Bio import Entrez

import config
import modules as m

# Configs



Entrez.email = config.ncbi.entrez_email


class Fasta():

  def __init__(self,
    fpath,
    valid_taxons_fpath = None
  ):
    self.fpath = fpath
    self.valid_taxons_fpath = valid_taxons_fpath

  # Combines fasta files into a single one.
  # Accepts names as a list of directory names in Genomes directory or a list of paths.
  @staticmethod
  def combine(fasta_fpaths, out_fpath):
    with open(out_fpath, 'w') as f_out:
      for fasta_fpath in fasta_fpaths:
        with open(fasta_fpath, 'r') as f_in:
          for line in f_in:
            f_out.write(line)

  # Modify defline and remove empty lines at end of file.
  @staticmethod
  def clean(genome, fpath):
    idx = 0
    for line in fileinput.input(fpath, inplace=True):
      # This is the defline, so modify it
      if line[0] == ">":
        # fileinput redirects print to file
        print(f">{genome}-{idx}\n", end="")
        idx += 1
      # This is not the defline, so just write it
      elif line != "" and line != "\n":
        # fileinput redirects print to file
        print(line, end="")
      # Dont write empty lines
  

  ##############################################################################
  # Deflines
        
  def get_deflines(self):
    deflines = []
    with open(self.fpath, 'r') as f:
      for line in f:
        if line.startswith(">"):
          defline = line[1:].strip("\n")
          deflines.append(defline)
    return deflines


  def deflines(self, out_fpath=None):
    deflines = self.get_deflines(self.fpath)
    if out_fpath:
      with open(out_fpath, 'w') as f:
        for defline in deflines:
          f.write(f"{defline}\n")
    return self.get_deflines(self.fpath)
  

  def export_deflines(self,
    fpath = None, # Option 1: Specify filename
    dpath = None  # Option 2: Specfiy directory, uses default filename
  ):
    
    if dpath:
      out_fname = m.File.get_fname_without_extension(
        self.fpath
      ) + "_deflines.tsv"
      out_fpath = os.path.join(dpath, out_fname)

    print("Scanning for deflines...")
    with open(self.fpath, 'r') as f_in, \
      open(out_fpath, 'w') as f_out:
      print(f" Reading {self.fpath}...")
      defline_qty = 0
      for i, line in enumerate(f_in):
        sys.stdout.write("\r")
        if line.startswith(">"):
          defline_qty += 1
          f_out.write(line)
        if (i + 1) % 1_000_000 == 0:
          sys.stdout.write(f"  {i+1} lines read")
      sys.stdout.write(f"  {i+1} lines read")
      print(f"\n  defline_qty: {defline_qty}")

  # This concates files, eg, created by export_deflines()
  @staticmethod
  def concat_defline_files(
    fpaths: list,
    out_fpath: str,
    keep_originals: True
  ):
    with open(out_fpath, 'w') as f_out:
      for fpath in fpaths:

        fa_fname = m.File.get_fname_without_extension(
          fpath
        ).replace("_deflines", ".fa")

        with open(fpath, 'r') as f_in:
          for line in f_in:
            f_out.write(f"{fa_fname}\t{line}")
        if not keep_originals:
          m.File.delete(fpath)


  ##############################################################################

  
  def lineages(self,
    accession_id_re_pattern = None,
    omit = None # omit these even if they satisfy the re pattern
  ):
    print(omit)
    # Get accession IDs from Fasta file
    accession_ids = set()
    with open(self.fpath, 'r') as f:
      for line in f:
        if line.startswith(">"):
          defline = line[1:].strip("\n")
          accession_id = None
          if accession_id_re_pattern:
            match = re.search(accession_id_re_pattern, defline)
            accession_id = match.group()
          if accession_id and accession_id not in omit:
            accession_ids.add(accession_id)
    if not len(accession_ids) > 0:
      return None
    

    webenv, query_key = m.Entrez.get_epost_webenv_and_query_key(
      db="nuccore",
      id=list(accession_ids)
    )

    i = 0
    retstart = 0
    retmax = 500
    while i < len(accession_ids):
      esummary_handle = m.Entrez.esummary({
        "db": "nuccore",
        "retstart": retstart,
        "retmax": retmax,
        "webenv": webenv,
        "query_key": query_key
      })
      esummary_json = json.loads(esummary_handle.readline().decode('utf-8'))
      return esummary_json
    
    # print(f"accession_id: {accession_id}")
    # handle = Entrez.esummary(
    #   db="nuccore",
    #   id=accession_id,
    #   retmode="xml"
    #   )
    # record = Entrez.read(handle)

    # lineages = {
    #   # accession_id: [
    #   #   Taxon ID,
    #   #   Lineage
    #   # ]
    # }

    # return record
  

  # def validate_deflines(self):