
import fileinput

import json
import re
import sys

from Bio import Entrez

import config as c
import modules as m


class Fasta():

  def __init__(self,
    fpath: str,
    valid_taxons_fpath: str|None = None
  ):
    self.fpath = fpath
    self.valid_taxons_fpath = valid_taxons_fpath

  
  @staticmethod
  def inspect_deflines(
    fpath: str,
    validate: list[str] = ['ncbi_accession_id'],
    out_fpath: str|None = None,
    inspection_fpath: str|None = None
  ) -> dict:
    """_summary_
    Example defline:
      >NZ_NIYS01000083.1 Shigella boydii strain ESBL-W3-2 NODE_67_length_21055_cov_46.3416_ID_133.ctg_1, whole genome shotgun sequence
    
    This currently doesn't account for GenBank accession prefixes such as: PAAA-PZZZ	NCBI	WGS	4+8 or more. Therefore, they will be identified as invalid even though they are likely valid.

    Args:
        fpath (str): _description_
        validate (list[str], optional): _description_. Defaults to 'ncbi_accession_id'.

    Returns:
        dict: A dictionary that includes refseq metrics, a list of invalid data, etc.
    """
    
    data = {
      'defline_qty': 0,
      'invalid_deflines': '',
      'refseq': {},
      'refseq_qty': 0,
      'non_refseq_qty': 0,
      'genbank': {},
      'no_genbank_prefix_qty': 0
    }

    for refseq_prefix in c.ncbi.refseq.accession_prefixes:
      data['refseq'][refseq_prefix] = {
        'qty': 0
      }

    if not out_fpath:
      out_fpath = m.file_sys.File.get_fpath_without_extension(
        fpath
      ) + '-accession-ids.tsv'

    with open(fpath, 'r') as f_in, \
      open(out_fpath, 'w') as f_out:

      for i, line in enumerate(f_in, start=1):

        # Skip non deflines
        if not line.startswith('>'):
          continue

        data['defline_qty'] += 1

        accession_id = line.split(' ')[0][1:]
        # print(f'accession_id: {accession_id}')
        refseq_prefix = accession_id[0:3] # NZ_ etc.


        # Check for RefSeq ID
        
        # If the first 3 characters is a RefSeq prefix
        if refseq_prefix in c.ncbi.refseq.accession_prefixes:
          data['refseq_qty'] += 1
          data['refseq'][refseq_prefix]['qty'] += 1
        
        # If the 3rd character is not an underscore, it should be in the list of
        # RefSeq prefixes.
        elif refseq_prefix[-1:] == '_':
          # It could be an invalid RefSeq prefix
          data['invalid_deflines'] += str(i) + ','

        # Otherwise, it's not a RefSeq prefix
        else:
          data['non_refseq_qty'] += 1
        

        # Check GenBank ID

        # Alias genbank_access_id 
        gb_acc_id = accession_id
        if '_' in accession_id:
          gb_acc_id = accession_id.split('_')[1]
        
        gb_pref = ''.join(filter(str.isalpha, gb_acc_id))
        
        # If there is no GenBank prefix
        if gb_pref == '':
          data['no_genbank_prefix_qty'] += 1
        
        elif gb_pref in c.ncbi.genbank.accession_prefixes:
          if gb_pref not in data['genbank']:
            data['genbank'][gb_pref] = {'qty': 1}
          else:  
            data['genbank'][gb_pref]['qty'] += 1

        else:
          data['invalid_deflines'] += str(i) + ','

        out_str = []
        if 'ncbi_accession_id' in validate:
          out_str.append(accession_id )
        f_out.write('\t'.join(out_str) + '\n')

    # Remove trailing ','
    data['invalid_deflines'] = data['invalid_deflines'][:-1]
    
    if not inspection_fpath:
      inspection_fpath = m.file_sys.File.get_fpath_without_extension(
        fpath
      ) + '-inspection.json'

    m.json.save_to_file(data, inspection_fpath)

    return data



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
  

  def export(self,
    out_dpath: str|None = None,
    seq_hash_type: str|None = m.hash.SHA256
  ) -> dict:
    """_summary_

    Args:
        is_deflines (bool, optional): _description_. Defaults to True.
        seq_hash_type (str | None, optional): _description_. Defaults to None.

    Returns:
        dict: _description_
    """

    if not out_dpath:
      export_fpath = m.file_sys.File.get_fpath_without_extension(
        self.fpath
      ) + '-export.tsv'
      meta_fpath = m.file_sys.File.get_fpath_without_extension(
        export_fpath
      ) + '-meta.json'
    
    meta = m.dict.create_dotdict()
    meta.fasta_fpath = self.fpath
    meta.export_fpath = export_fpath
    meta.meta_fpath = meta_fpath
    meta.defline_qty = 0
    meta.seq_qty = 0
    meta.seq_hash_type = seq_hash_type
    meta.total_bp = 0

    block = m.dict.create_dotdict()
    block.out_line = ''
    block.defline = ''
    block.seq = []
    block.seq_length = 0

    def write_line(f_out, block, meta):
      
      # Concat sequence
      block.seq = ''.join(block.seq)

      # Get sequence length
      block.seq_length = len(block.seq)

      # Handle meta
      meta.defline_qty += 1
      meta.seq_qty += 1
      meta.total_bp += block.seq_length

      # Make hash
      hash_str = m.hash.get_hash(block.seq, meta.seq_hash_type)
      
      # Build output line
      block.out_line += f'\t{block.seq_length}\t{hash_str}\n'
      
      # Write output line
      f_out.write(block.out_line)
      
      # Reset block data
      block.out_line = ''
      block.defline = ''
      block.seq = []
      block.seq_length = 0

      return block, meta

    # Export data
    with open(self.fpath, 'r') as f_in, \
      open(export_fpath, 'w') as f_out:
      print(f"Reading {self.fpath}...")

      is_seq = False
      for i, line in enumerate(f_in):
        sys.stdout.write("\r")

        if line.startswith(">"):
          
          # Write line & reset block
          if is_seq:
            is_seq = False
            # Write output line
            block, meta = write_line(f_out, block, meta)

          # Get defline
          block.defline = line.rstrip('\n')
          
          # Add defline to output line
          block.out_line += block.defline

          # Next line(s) are sequences
          is_seq = True
        
        # Build sequence string
        elif is_seq:
          block.seq.append(line.rstrip('\n'))

        if (i + 1) % 1_000_000 == 0:
          sys.stdout.write(f" {i+1} lines read")
      
      # Write last line
      block, meta = write_line(f_out, block, meta)

      sys.stdout.write(f" {i+1} lines read")
      print(f"\n  defline_qty: {meta.defline_qty}")

    # Save meta data
    m.json.save_to_file(meta, meta_fpath)

    return meta

  # This concates files, eg, created by export_deflines()
  @staticmethod
  def concat_defline_files(
    fpaths: list,
    out_fpath: str,
    keep_originals: True
  ):
    with open(out_fpath, 'w') as f_out:
      for fpath in fpaths:

        fa_fname = m.file_sys.File.get_fname_without_extension(
          fpath
        ).replace("_deflines", ".fa")

        with open(fpath, 'r') as f_in:
          for line in f_in:
            f_out.write(f"{fa_fname}\t{line}")
        if not keep_originals:
          m.file_sys.File.delete(fpath)


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
    

    webenv, query_key = m.ncbi.Entrez.get_epost_webenv_and_query_key(
      db="nuccore",
      id=list(accession_ids)
    )

    i = 0
    retstart = 0
    retmax = 500
    while i < len(accession_ids):
      esummary_handle = m.ncbi.Entrez.esummary({
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