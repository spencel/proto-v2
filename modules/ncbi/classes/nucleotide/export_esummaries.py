
import os

from modules import file_sys as m_file_sys

from .get_summary_batch import *


def export_esummaries(cls,
  accession_ids_fpath: str,
  col_idx: int = 0,
  export_items: set = {'taxid', 'gi', 'accessionversion', 'organism'},
  out_fpath: str|None = None
) -> int:
  """_summary_

  Args:
      accession_ids_fpath (str): _description_
      export_items (set | None, optional): _description_. Defaults to None.
      out_fpath (str | None, optional): _description_. Defaults to None.
  """
  if not out_fpath:
    out_fpath = os.path.join(
      m_file_sys.File.get_fpath_without_extension(accession_ids_fpath)
      + '-esummaries.tsv'
    )
  
  esummaries = get_esummary_batch(cls,
    accession_ids_fpath = accession_ids_fpath
  )
  print(esummaries)

  esummary_qty = 0
  with open(out_fpath, 'w') as f:
    # Make column name row
    f.write('\t'.join(export_items) + '\n')
    for gi_id, esummary in esummaries.items():
      esummary_qty += 1
      out_line = []
      for export_item in export_items:
        if export_item in esummary:
          out_line.append(str(esummary[export_item]))
        else:
          raise Exception(f'Error: {export_item} missing from esummary: {esummary}')
      f.write('\t'.join(out_line) + '\n')  

  return esummary_qty