
import datetime
import json
import logging
from typing import Union, List

from .. import Entrez

def _get_esearch_res(
  ids,
  id_idx,
  id_qty,
  esearch_term_max_len
) -> Union[int, List[str], str]:
  
  # Build search term for ESearch
  query_ids = []
  # -1 because last character would be a comma, which is removed
  esearch_term_len = -1 
  while id_idx < id_qty:
    id = ids[id_idx]
    esearch_term_len += len(f'{id},')

    # Finish adding accession IDs to query if max query length would be
    #   exceeded.
    if esearch_term_len > esearch_term_max_len:
      break
    
    # Otherwise, add this ID to the query list
    query_ids.append(id)
    id_idx += 1

  esearch_term = ','.join(query_ids)

  esearch_res = Entrez.esearch({
    'db': 'nuccore',
    'term': esearch_term
  })

  return id_idx, query_ids, esearch_res


def _get_epost_res(
  uids
):

  epost_uids = ','.join(uids)

  epost_res = Entrez.epost({
    'db': 'nuccore',
    'id': epost_uids
  })

  return epost_res


# Alias:
#   acc = accession_id, eg: 'JBAJFE010000013.1'
#   len = length
def get_esummaries(
  gi_ids: list[str] = None,
  accession_ids: list[str] = None,
  gi_and_accession_ids: list[str] = None
) -> dict:
  """_summary_

  Args:
    accession_ids (list[str], optional): _description_. Defaults to None.
      Example: ['JBAJFE010000013.1', ...]
  
  Returns:
    _type_: _description_
  """
  esummaries = []

  # Handle ID input
  epost_max_query_ids = 1_000_000_000_000
  all_ids = []
  if gi_ids:
    all_ids.extend(gi_ids)
  if accession_ids:
    all_ids.extend(accession_ids)
    epost_max_query_ids = 500
  if gi_and_accession_ids:
    all_ids.extend(gi_and_accession_ids)
    epost_max_query_ids = 500
  
  id_idx = 0
  id_qty = len(all_ids)

  while id_idx < id_qty:

################################################################################
## Get UIDs using ESearch

    uids = []

    while len(uids) == 0:
      id_idx, query_ids, esearch_res = _get_esearch_res(
        all_ids,
        id_idx,
        id_qty,
        esearch_term_max_len = 4000
      )
      uids = esearch_res['IdList']
      if len(uids) == 0:
        # Rollback id_idx
        id_idx -= len(query_ids)
        print(f'Failed  to get UIDs. Retrying ESearch...')

################################################################################
## EPost IDs and get WebEnv/Query
  
    # Get EPost term
    epost_res = _get_epost_res(
      uids
    )

    webenv = epost_res["WebEnv"]
    query_key = epost_res["QueryKey"]

################################################################################
## Get ESummaries

    esummary_res = Entrez.esummary(
      esummary_params = {
        'db': 'nuccore',
        'webenv': webenv,
        'query_key': query_key
      },
      as_json = True
    )

    for uid, items in esummary_res['result'].items():
      # Skip UID list
      if uid == 'uids':
        continue
      # Otherwise, add essumary
      esummaries.append(items)
    
    print(f'id_idx: {id_idx}')

  return esummaries

################################################################################


def get_esummaries_old(
  accession_ids: list[str] = None,
  taxon_id=None,
  clade_name=None,
  query=None,
  retmax=9999,
  retmode="json",
  time_frame=None,
  sort_by=None,
  get_complete_genome=True
) -> dict:
  """_summary_

  Args:
    taxon_id (_type_, optional): _description_. Defaults to None.
    clade_name (_type_, optional): _description_. Defaults to None.
    query (_type_, optional): _description_. Defaults to None.
    retstart (int, optional): _description_. Defaults to 0.
    retmax (int, optional): _description_. Defaults to 9999.
    retmode (str, optional): _description_. Defaults to "json".
    time_frame (_type_, optional): _description_. Defaults to None.
    sort_by (_type_, optional): _description_. Defaults to None.
    get_complete_genome (bool, optional): _description_. Defaults to True.

  Returns:
    _type_: _description_
  """
  esearch_results: list = []
  esummary_results: dict = dict()

  if retmode == "json":
    retmax = 500

  term: str
  
  if accession_ids:

    # The max qty of characters a query term can be
    max_term_length = 4041
    accession_id_qty = len(accession_ids)
    id_idx = 0
    while id_idx < accession_id_qty:
      
      next_length = 0
      query_ids = []
      
      while id_idx < accession_id_qty:
        accession_id = accession_ids[id_idx]
        next_length += len(accession_id) + 1
        
        if next_length > max_term_length:
          break
        
        query_ids.append(accession_id)
        # Add +1 to account for commas
        id_idx += 1

      term = ','.join(query_ids)

      # Initiliaze amount returned to start the loop
      retcount = retmax + 1

      while retcount > retmax:

        def esearch(args):
          return Entrez.esearch(args)

        esearch_result = esearch({
          "db": "nuccore",
          "retmax": retmax,
          "term": term,
          "idtype": "gi",
          "sort": sort_by
        })
        retcount = int(esearch_result['Count'])
        
        # Try again if return count is zero
        if retcount == 0:
          retcount = retmax + 1
          continue

        esearch_results.append(esearch_result)

    # Exit function if there are no results
    if esearch_results[0]["Count"] == "0":
      return None

    # Get the Accession IDs from the results
    gi_ids = []
    for esearch_result in esearch_results:
      gi_ids.extend(esearch_result["IdList"])
    gi_id_qty = len(gi_ids)

    if len(gi_ids) == 0:
      return "NO ID LIST RETURNED BY QUERY"

    # Upload a list of GI IDs to Entrez for the batch request
    epost_results = Entrez.epost({
      "db": "nuccore",
      "id": gi_ids
    })

    # Get the identifiers of the batch request
    webenv = epost_results["WebEnv"]
    query_key = epost_results["QueryKey"]
    
    # Initiliaze amount returned to start the loop
    esummary_qty = 0
    uid_qty = 0
    retstart = 0
    
    while esummary_qty < gi_id_qty:

      def esummary(args):
        return json.loads(
          Entrez.esummary(args).read()
        )
    
      esummary_response = esummary({
          "db": "nuccore",
          "retmax": retmax,
          "retmode": retmode,
          "retstart": retstart,
          "webenv": webenv,
          "query_key": query_key
      })
      uid_qty = len(esummary_response['result']['uids'])
        
      # Try again if return count is zero
      if uid_qty == 0:
        continue
      
      esummary_result = esummary_response['result']
      del esummary_result['uids']

      for gi_id, result in esummary_result.items():
        esummary_results[gi_id] = result

      esummary_qty += uid_qty
      retstart += uid_qty 

    return esummary_results

  else:
    term = ''
    # Special handling of queries
    if query:
      term = query
    else:
      if taxon_id:
        term = f"txid{taxon_id}[Organism]"
      elif clade_name:
        term = f"{clade_name}[Organism]"

    if not query and get_complete_genome:
      term += " AND \"complete genome\"[Title]"

    query_time_frame: str = None
    if time_frame:
      # Eg, past 6 months
      # ("2022/10/19"[PDAT] : "2023/04/19"[PDAT])
      if time_frame.startswith("past"):
        end_date = datetime.datetime.today()
        time_frame = time_frame.split(" ")
        units = time_frame[2]
        if units in ("monnth", "months"):
          time_offset = 365 / 12 * int(time_frame[1])
          start_date = end_date - datetime.timedelta(days=time_offset)

      query_time_frame = str(
        f"(\"{start_date.strftime('%Y/%m/%d')}\"[PDAT] : \"{end_date.strftime('%Y/%m/%d')}\"[PDAT])"
      )
    if query_time_frame:
      term += " AND " + query_time_frame
    
    # Get a list of the most recent records
    esearch_results = Entrez.esearch({
      "db": "nuccore",
      "retmax": retmax,
      "retstart": retstart,
      "term": term,
      "idtype": "gi",
      "sort": sort_by
    })
    json.dump(esearch_results, open('logs/esearch_results.log', 'w'), indent=2)

    # Exit function if there are no results
    if esearch_results["Count"] == "0":
      return None

    # Get the Accession IDs from the results
    gi_ids = esearch_results["IdList"]

    if len(gi_ids) == 0:
      return "NO ID LIST RETURNED BY QUERY"

    # Upload a list of GI IDs to Entrez for the batch request
    epost_results = Entrez.epost({
      "db": "nuccore",
      "id": gi_ids
    })

    # Get the identifiers of the batch request
    webenv = epost_results["WebEnv"]
    query_key = epost_results["QueryKey"]

    http_res_json = json.loads(
        Entrez.esummary({
          "db": "nuccore",
          "retmode": retmode,
          "webenv": webenv,
          "query_key": query_key
      }).read()
    )

    esummaries_json = http_res_json['result']
    del esummaries_json['uids']

    return esummaries_json
