
import logging
import os
import time

from modules import file_sys as m_file_sys

from .get_esummaries import *


log = logging.getLogger(__name__)
debug = log.debug


def export_esummaries(cls,
  accession_ids_fpath: str,
  col_idx: int = 0,
  export_items: list = ['taxid', 'gi', 'accessionversion', 'organism'],
  out_fpath: str|None = None,
  batch_qty: int = 200 # 200 appears to be stable for esummary
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
  elapsed_times = []
  esummary_added_qtys = []
  
  # Resume progress
  saved_accession_ids = set()
  if os.path.isfile(out_fpath):
    with open(out_fpath) as f:
      # Skip first line
      for line in f:
        break
      # Continue
      for line in f:
        saved_accession_id = line.rstrip('\n').split('\t')[2]
        saved_accession_ids.add(saved_accession_id)
  else:
    # Make column name row
    with open(out_fpath, 'w') as f:
      f.write('\t'.join(export_items) + '\n')
  debug(f'len(saved_accession_ids): {len(saved_accession_ids)}')
  
  esummary_qty = len(saved_accession_ids)
  debug(f'esummary_qty: {esummary_qty}')

  accession_ids = []
  if accession_ids_fpath:
    with open(accession_ids_fpath, 'r') as f:
      for line in f:
        accession_id = line.rstrip('\n')
        if accession_id not in saved_accession_ids:
          ar_line = line.strip('\n').split('\t')
          accession_id = ar_line[col_idx]
          accession_ids.append(accession_id)
  
  del saved_accession_ids

  accession_id_qty = len(accession_ids)
  debug(f'accession_id_qty: {accession_id_qty}')
  remaining_accession_id_qty = accession_id_qty
  
  retmax = 500
  # The max qty of characters a query term can be
  # max_term_length = 4041
  max_term_length = 2000
  id_idx = 0
  while id_idx < accession_id_qty:
    print(f'id_idx: {id_idx}')
    start_time = time.time()
    
    next_length = 0
    query_ids = []
    
    while id_idx < accession_id_qty:
      accession_id = accession_ids[id_idx]
      print(f'accession_id: {accession_id}')
      next_length += len(accession_id) + 1
      
      if next_length > max_term_length \
        or len(query_ids) >= batch_qty:
        break
      
      query_ids.append(accession_id)
      # Add +1 to account for commas
      id_idx += 1
    debug(f'next_length:{next_length}')
    debug(f'id_idx:{id_idx}')

    debug(f'len(query_ids): {len(query_ids)}')
    term = ','.join(query_ids)
    # debug(f'term: {term}')

    # Initiliaze amount returned to start the loop
    esearch_retcount = retmax + 1

    while esearch_retcount > retmax:
      print(f'esearch_retcount: {esearch_retcount}')

      # Download GI IDs
      def entrez_esearch(args):
        try:
          return Entrez.esearch(args)
        except:
          debug(f'Entrez.search(): Error')
          return None

      params = {
        "db": "nuccore",
        "retmax": retmax,
        "term": term,
        "idtype": "gi"
      }
      esearch_response = None
      while esearch_response == None:
        esearch_response = entrez_esearch(params)
        print(f'esearch_response: {esearch_response}')
      
      esearch_result = esearch_response
      
      esearch_retcount = int(esearch_result['Count'])
      debug(f'esearch_retcount: {esearch_retcount}')
      
      # Try again if return count is zero
      if esearch_retcount == 0:
        esearch_retcount = retmax + 1
        continue

      # Get the Accession IDs from the results
      gi_ids = esearch_result["IdList"]
      gi_id_qty = len(gi_ids)
      debug(f'gi_id_qty: {gi_id_qty}')

      # Upload a list of GI IDs to Entrez for the batch request
      epost_results = Entrez.epost({
        "db": "nuccore",
        "id": gi_ids
      })

      # Get the identifiers of the batch request
      webenv = epost_results["WebEnv"]
      query_key = epost_results["QueryKey"]
      
      # Initiliaze amount returned to start the loop
      uid_qty = 0
      esummary_retstart = 0
      while esummary_retstart < gi_id_qty:
        print(f'esummary_retstart:{esummary_retstart}')

        # Download ESumaries
        def entrez_esummary(args):
          try:
            return json.loads(Entrez.esummary(args).read())
          except:
            debug(f'Entrez.esummary(): Error')
            return None

        params = {
          "db": "nuccore",
          "retmax": retmax,
          "retmode": 'json',
          "retstart": esummary_retstart,
          "webenv": webenv,
          "query_key": query_key
        }
        esummary_response = None
        while esummary_response == None:
          esummary_response = entrez_esummary(params)
          print(f'esummary_response: {esummary_response}')

        esummary_result = esummary_response['result']
        # debug(f'esummary_result: {esummary_result}')

        if 'uids' in esummary_result:
          del esummary_result['uids']

        uid_qty += len(esummary_response['result'])
        debug(f'uid_qty: {uid_qty}')
          
        # Try again if return count is zero
        if uid_qty == 0:
          continue

        esummary_retstart += uid_qty
        debug(f'esummary_retstart: {esummary_retstart}')

        # Write ESummaries
        with open(out_fpath, 'a') as f:
          for gi_id, esummary in esummary_result.items():
            print(f'gi_id: {gi_id}')
            esummary_qty += 1
            out_line = []
            for export_item in export_items:
              if export_item in esummary:
                out_line.append(str(esummary[export_item]))
              else:
                raise Exception(f'Error: {export_item} missing from esummary: {esummary}')
            f.write('\t'.join(out_line) + '\n')
        debug(f'esummary_qty: {esummary_qty}')
        
        # Metrics
        esummary_added_qtys.append(uid_qty)
        avg_esummary_added_qty = sum(esummary_added_qtys) / len(esummary_added_qtys)
        print(f' avg_esummary_added_qty: {avg_esummary_added_qty}')
        elapsed_times.append(time.time() - start_time)
        avg_elapsed_time = sum(elapsed_times) / len(elapsed_times)
        print(f' avg_elapsed_time: {avg_elapsed_time}') 
        remaining_accession_id_qty -= uid_qty
        print(f' remaining_accession_id_qty: {remaining_accession_id_qty}')
        time_remaining = (remaining_accession_id_qty / avg_esummary_added_qty * avg_elapsed_time) / 60 / 60
        print(f' time_remaining: {time_remaining}') 

  return esummary_qty