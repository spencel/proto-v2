
import os

from modules import file_sys


DB_UPDATE_SPECIES_TAXON_IDS_FPATH = os.path.join(
	'data', 'data/db-update-species-taxon-ids.tsv'
)

def get_taxon_ids_from_file(
	deflines_fpaths: list[str],
  taxon_ids_fpath: str = DB_UPDATE_SPECIES_TAXON_IDS_FPATH,
  out_fpath: str|None = None
):
  """This script is unfinished. It's uncertain if the defline can be used to with running into inadvertent duplicates.

  Args:
      deflines_fpaths (list[str]): _description_
      taxon_ids_fpath (str, optional): _description_. Defaults to DB_UPDATE_SPECIES_TAXON_IDS_FPATH.
      out_fpath (str | None, optional): _description_. Defaults to None.
  """
	
  if not out_fpath:
    out_fpath = file_sys.File.get_fpath_without_extension(taxon_ids_fpath) + 'deflines.tsv'

  species_X_taxon_ids: dict = dict()
  with open(taxon_ids_fpath) as f:
    # Skip first line
    next(f)
    # Continue
    for line in f:
      ar_line = line.rstrip('\n').split('\t')
      species = ar_line[0]
      taxon_id = ar_line[1]
      species_X_taxon_ids[species] = taxon_id