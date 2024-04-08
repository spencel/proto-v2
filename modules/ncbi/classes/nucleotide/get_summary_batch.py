
import datetime
import json

from .. import Entrez


def get_esummary_batch(cls,
  accession_ids_fpath: str = None,
  taxon_id=None,
  clade_name=None,
  query=None,
  retstart=0,
  retmax=9999,
  retmode="json",
  time_frame=None,
  sort_by=None,
  get_complete_genome=True
):
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

  if retmode == "json":
    retmax = 500

  term: str = None

  accession_ids = []
  if accession_ids_fpath:
    with open(accession_ids_fpath, 'r') as f:
      for line in f:
        accession_ids.append(line[:-1])
    
    term = ','.join(accession_ids)
    del accession_ids
  
  else:
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

  with open('log.log', 'w') as f:
    json.dump(esummaries_json, f, indent=2)
  print(esummaries_json)

  return esummaries_json
