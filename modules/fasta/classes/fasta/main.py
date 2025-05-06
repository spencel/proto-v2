
import fileinput

import json
import re
import sys

import config as c
from modules import file_sys, hash
from modules import list as m_list
from modules import json as m_json
from modules import dict as m_dict
from modules import ncbi


def print_file_read_status(line_idx:int, interval:int=1_000_000) -> None:
  if (line_idx + 1) % interval == 0:
    sys.stdout.write(f"\033[F{line_idx + 1} lines read\n")


    
class DEFLINE_SELECTION:
  ALL = 'ALL'
  MITOCHONDRIA = 'MITOCHONDRIA'
  NON_GENBANK = 'NON_GENBANK'
  NON_MITOCHONDRIA = 'NON_MITOCHONDRIA'
  NON_REFSEQ = 'NON_REFSEQ'

class DefLine():

  MITOCHONDRIA_KEY_WORDS = {
    'mitochondria',
    'mitochondrial scaffold',
    'mitochondrion'
  }

  @classmethod
  def get_accession_id(cls, defline):
    return defline.split(' ')[0][1:]

class Fasta():

  def __init__(self,
    fpath: str = '',
    valid_taxons_fpath: str = ''
  ):
    self.fpath = fpath
    self.valid_taxons_fpath = valid_taxons_fpath
    self.deflines: list = []
    self.accession_ids: list = []

  def inspect_deflines(self,
    validate: list[str] = ['ncbi_accession_id'],
    inspection_fpath: str = ''
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
      'refseq': {
        'qty': 0,
        'non_qty': 0,
        'invalid_qty': 0,
        'items': {}
      },
      'genbank': {
        'qty': 0,
        'non_qty': 0,
        'items': {}
      }
    }

    for refseq_pref in c.ncbi.refseq.accession_prefixes:
      data['refseq']['items'][refseq_pref] = {
        'qty': 0
      }

    with open(self.fpath, 'r') as f:
      

      for i, line in enumerate(f, start=1):

        # Skip non deflines
        if not line.startswith('>'):
          print_file_read_status(i, 1_000_000)
          continue

        data['defline_qty'] += 1

        accession_id = DefLine.get_accession_id(line)
        
        # Check for RefSeq ID
        if '_' in accession_id:
          refseq_pref = accession_id.split('_')[0] # NZ, etc.
          if refseq_pref in c.ncbi.refseq.accession_prefixes:
            data['refseq']['items'][refseq_pref]['qty'] += 1
            data['refseq']['qty'] += 1
          else:
            data['refseq']['invalid_qty'] += 1
        else:
          data['refseq']['non_qty'] += 1
        
        # Check GenBank ID
        # Alias genbank_access_id 
        gb_acc_id = accession_id
        if '_' in accession_id:
          gb_acc_id = accession_id.split('_')[1]
        
        gb_pref = ''.join(filter(str.isalpha, gb_acc_id))
        
        # If there is no GenBank prefix
        if gb_pref:
          if gb_pref not in data['genbank']['items']:
            data['genbank']['items'][gb_pref] = {'qty': 1}
          else:
            data['genbank']['items'][gb_pref]['qty'] += 1
          data['genbank']['qty'] += 1
        else:
          data['genbank']['non_qty'] += 1

        print_file_read_status(i, 1_000_000)

    return data

  def export_defline_inspection(self,
    out_fpath: str = '',
    validate: list[str] = ['ncbi_accession_id'],
    inspection_fpath: str = ''
  ):
    
    if not out_fpath:
      out_fpath = file_sys.File.get_fpath_without_extension(
        self.fpath
      ) + '-inspection.json'
    
    data = self.inspect_deflines(
      validate = validate,
      inspection_fpath = inspection_fpath
    )
    
    # with open(out_fpath, 'w') as f:
    #   out_str = []
    #   if 'ncbi_accession_id' in validate:
    #     out_str.append(accession_id)
    #   f_out.write('\t'.join(out_str) + '\n')

    json.dump(data, open(out_fpath, 'w'), indent=2)

  # Combines fasta files into a single one.
  # Accepts names as a list of directory names in Genomes directory or a list of paths.
  @staticmethod
  def combine(fasta_fpaths, out_fpath):
    with open(out_fpath, 'w') as f_out:
      for fasta_fpath in fasta_fpaths:
        with open(fasta_fpath, 'r') as f_in:
          for line in f_in:
            f_out.write(line)

  @staticmethod
  def separate(fasta_fpath: str, limit: int = 0):
    if limit: 
      print(True)
    else:
      print(False)

    """_summary_

    Args:
        fasta_fpath (str): _description_
        limit (int, optional): Limit number of files to generate. Defaults to 0.
    """
    i_out_file = 0
    with open(fasta_fpath) as f_in:
      out_file = None
      for line in f_in:
        if line.startswith('>'):
          if limit and i_out_file >= limit:
            break
          if out_file:
            out_file.close()
          out_fpath = str(
            f'{file_sys.remove_extension(fasta_fpath)}'
            f'-{i_out_file}'
            '.fa'
          )
          print(f'Writing to file No. {i_out_file}...')
          out_file = open(out_fpath, 'w')
          i_out_file += 1
        out_file.write(line)
      if out_file:
        out_file.close()

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
  
  def copy(self,
    ncbi_accession_ids: set[str]
  ):
    
    with open(self.fpath) as f_in:
      f_out = None
      ncbi_accession_id = ''
      has_ncbi_accession_id = False
      for line in f_in:
        if line.startswith('>'):
          if f_out:
              f_out.close()
              f_out = None
          for ncbi_accession_id in ncbi_accession_ids:
            if ncbi_accession_id in line:
              has_ncbi_accession_id = True
              break
            has_ncbi_accession_id = False
          if has_ncbi_accession_id:
            out_fpath = f'{ncbi_accession_id}.fa'
            f_out = open(out_fpath, 'w')
        if f_out:
          f_out.write(line)
      if f_out:
        f_out.close()
  
  def copy_sequence(self,
    start_locus: int,
    stop_locus: int
  ):
    
    cur_locus = 1 # first base is 1 for sam file loci
    with open(self.fpath) as f_in, \
    open(self.fpath + '.sub_seq', 'w') as f_out:
      defline = next(f_in).strip('\n')
      f_out.write(f'{defline}, loci {start_locus} to {stop_locus}\n')
      for line in f_in:
        line = line.rstrip('\n')
        for char in line:
          if start_locus <= cur_locus and cur_locus <= stop_locus:
            f_out.write(char)
          cur_locus += 1        
    

################################################################################
# Deflines
        
  def load_deflines(
    self,
    selections: list[str] = [DEFLINE_SELECTION.ALL]
  ):
    with open(self.fpath, 'r') as f:
      for i, line in enumerate(f):
        if line.startswith(">"):
          defline = line.strip("\n")
          
          is_load_defline = True

          if DEFLINE_SELECTION.NON_MITOCHONDRIA in selections:
            for keyword in DefLine.MITOCHONDRIA_KEY_WORDS:
              if keyword in defline.lower():
                is_load_defline = False
          
          if DEFLINE_SELECTION.MITOCHONDRIA in selections:
            is_load_defline = False
            for keyword in DefLine.MITOCHONDRIA_KEY_WORDS:
              if keyword in defline.lower():
                is_load_defline = True
          
          if is_load_defline:
            self.deflines.append(defline)

        print_file_read_status(i)

  def get_accession_ids(
    self,
    selections: list[str] = [DEFLINE_SELECTION.ALL]
  ):
    if not self.deflines:
      self.load_deflines(
        selections = selections
      )
    for defline in self.deflines:
      self.accession_ids.append(DefLine.get_accession_id(defline))

  def export_deflines(self,
    out_fpath: str = '',
    selections: list[str] = [DEFLINE_SELECTION.ALL]
  ) -> None:
    
    if not out_fpath:
      out_fpath = file_sys.File.get_fpath_without_extension(
        self.fpath
      ) + '-deflines.txt'
    if not self.deflines:
      self.load_deflines(
        selections = selections
      )

    if DEFLINE_SELECTION.ALL in selections :
      m_list.export(self.deflines, out_fpath)
      return 
    
    selected_deflines = []

    if DEFLINE_SELECTION.NON_GENBANK in selections:
      for defline in self.deflines:
        accession_id = DefLine.get_accession_id(defline)
        # Remove refseq prefix if present
        if '_' in accession_id:
          accession_id = accession_id.split('_')[-1]
        gb_pref = ''.join(filter(str.isalpha, accession_id))
        if not gb_pref:
          selected_deflines.append(defline)

    if DEFLINE_SELECTION.NON_REFSEQ in selections:
      for defline in self.deflines:
        accession_id = DefLine.get_accession_id(defline)
        refseq_pref: str = ''
        if '_' in accession_id:
          refseq_pref = accession_id.split('_')[0]
        if refseq_pref and refseq_pref not in c.ncbi.refseq.accession_prefixes:
          if defline not in selected_deflines:
            selected_deflines.append(defline)
        if not refseq_pref and defline not in selected_deflines:
          selected_deflines.append(defline)

    if DEFLINE_SELECTION.MITOCHONDRIA in selections:
      for defline in self.deflines:
        is_export_defline = False
        for keyword in DefLine.MITOCHONDRIA_KEY_WORDS:
          if keyword in defline.lower():
            is_export_defline = True
        if is_export_defline:
          selected_deflines.append(defline)
    
    if DEFLINE_SELECTION.NON_MITOCHONDRIA in selections:
      for defline in self.deflines:
        is_export_defline = True
        for keyword in DefLine.MITOCHONDRIA_KEY_WORDS:
          if keyword in defline.lower():
            is_export_defline = False
        if is_export_defline:
          selected_deflines.append(defline)

    m_list.export(selected_deflines, out_fpath)

  def export_defline_and_seq_hash(self,
    out_dpath: str|None = None,
    seq_hash_type: str|None = hash.SHA256
  ) -> dict:
    """_summary_

    Args:
        is_deflines (bool, optional): _description_. Defaults to True.
        seq_hash_type (str | None, optional): _description_. Defaults to None.

    Returns:
        dict: _description_
    """

    if not out_dpath:
      export_fpath = file_sys.File.get_fpath_without_extension(
        self.fpath
      ) + '-export.tsv'
      meta_fpath = file_sys.File.get_fpath_without_extension(
        export_fpath
      ) + '-meta.json'
    
    meta = m_dict.create_dotdict()
    meta.fasta_fpath = self.fpath
    meta.export_fpath = export_fpath
    meta.meta_fpath = meta_fpath
    meta.defline_qty = 0
    meta.seq_qty = 0
    meta.seq_hash_type = seq_hash_type
    meta.total_bp = 0

    block = m_dict.create_dotdict()
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
      hash_str = hash.get_hash(block.seq, meta.seq_hash_type)
      
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
    m_json.save_to_file(meta, meta_fpath)

    return meta

  def export_accession_ids(self,
    out_fpath: str = '',
    selections: list[str] = [DEFLINE_SELECTION.ALL]
  ) -> None:

    if not out_fpath:
      out_fpath = file_sys.File.get_fpath_without_extension(
        self.fpath
      ) + '-accession-ids.txt'
    if not self.accession_ids:
      self.get_accession_ids(
        selections = selections
      )
    m_list.export(self.accession_ids, out_fpath)

  # This concates files, eg, created by export_deflines()
  @staticmethod
  def concat_defline_files(
    fpaths: list,
    out_fpath: str,
    keep_originals: True
  ):
    with open(out_fpath, 'w') as f_out:
      for fpath in fpaths:

        fa_fname = file_sys.File.get_fname_without_extension(
          fpath
        ).replace("_deflines", ".fa")

        with open(fpath, 'r') as f_in:
          for line in f_in:
            f_out.write(f"{fa_fname}\t{line}")
        if not keep_originals:
          file_sys.File.delete(fpath)

# End Deflines
################################################################################

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
    

    webenv, query_key = ncbi.Entrez.get_epost_webenv_and_query_key(
      db="nuccore",
      id=list(accession_ids)
    )

    i = 0
    retstart = 0
    retmax = 500
    while i < len(accession_ids):
      esummary_handle = ncbi.Entrez.esummary({
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