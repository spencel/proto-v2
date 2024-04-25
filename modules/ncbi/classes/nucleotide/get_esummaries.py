
import datetime
import json
import logging

from .. import Entrez


log = logging.getLogger(__name__)
debug = log.debug


def get_esummaries(cls,
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
    debug(f'len(accession_ids): {len(accession_ids)}')

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
      debug(f'next_length:{next_length}')
      debug(f'id_idx:{id_idx}')

      debug(f'len(query_ids): {len(query_ids)}')
      term = ','.join(query_ids)
      debug(f'term: {term}')

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
        debug(f'retcount: {retcount}')
        
        # Try again if return count is zero
        if retcount == 0:
          retcount = retmax + 1
          continue

        debug(f'esearch_result: {esearch_result}')
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
      debug(f'esummary_response: {esummary_response}')
      uid_qty = len(esummary_response['result']['uids'])
      debug(f'uid_qty: {uid_qty}')
        
      # Try again if return count is zero
      if uid_qty == 0:
        continue
      
      esummary_result = esummary_response['result']
      del esummary_result['uids']
      debug(f'esummary_result: {esummary_result}')

      for gi_id, result in esummary_result.items():
        esummary_results[gi_id] = result

      esummary_qty += uid_qty
      debug(f'esummary_qty: {esummary_qty}')
      retstart += uid_qty
      debug(f'retstart: {retstart}')      

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
    
    debug(f'retmax: {retmax}')
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
