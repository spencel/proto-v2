
import json


class Taxonomy():


  @staticmethod
  def get_sub_taxon_ids(taxon_name, retstart=0):
    """
    Returns a separate set of all taxon IDs for each given parent taxon ID.
    :param taxon_ids:
    :return:
    """
    term = str(taxon_name) + "[Subtree]"

    # print(f"sl: term: {term}")

    esearch_json = m.ncbi.Entrez.esearch(
      {
        "db": "taxonomy",
        "retstart": retstart,
        "term": term
      }
    )

    return esearch_json


  @classmethod
  def get_species_of_genus(cls, taxon_name):

    species = dict()
    # i_record = 0

    # Therefore, assume there is at least 1 record at first.
    esearch_retstart = 0
    esearch_record_qty = 1
    while esearch_retstart < esearch_record_qty:
      esearch_json = cls.get_sub_taxon_ids(taxon_name, retstart=esearch_retstart)
      esearch_record_qty = len(esearch_json["IdList"])
      # print(f"sl: esearch_record_qty: {esearch_record_qty}")

      epost_json = m.ncbi.Entrez.epost({
        "db": "Taxonomy",
        "id": esearch_json["IdList"]
      })
      webenv = epost_json["WebEnv"]
      query_key = epost_json["QueryKey"]

      esummary_retstart = 0
      while esummary_retstart < esearch_record_qty:
        esummary_json = m.ncbi.Entrez.esummary({
          "db": "Taxonomy",
          "retstart": esummary_retstart,
          "webenv": webenv,
          "query_key": query_key
        }, as_json=True)
        esummary_record_qty = len(esummary_json["result"]["uids"])
        # print(f"sl: esummary_record_qty: {esummary_record_qty}")
        m.Log.write(f"esummary_json: {json.dumps(esummary_json, indent=2)}")

        for key in esummary_json["result"]:
          if key == "uids":
            continue
          # i_record += 1
          data = esummary_json["result"][key]
          if data["rank"] == "species":
            species[data["scientificname"]] = data["taxid"]

        esummary_retstart = esummary_retstart + esummary_record_qty

      esearch_retstart = esearch_retstart + esearch_record_qty

    # print(f"sl: i_record: {i_record}")

    return species


  @staticmethod
  def get_esummary(taxon_ids):

    webenv, query_key = m.ncbi.Entrez.get_epost_webenv_and_query_key(
      db="Taxonomy",
      id=taxon_ids
    )

    i_taxon = 0
    retstart = 0
    retmax = 500
    while i_taxon < len(taxon_ids):
      esummary_handle = m.ncbi.Entrez.esummary({
        "db": "taxonomy",
        "retstart": retstart,
        "retmax": retmax,
        "webenv": webenv,
        "query_key": query_key
      })
      esummary_json = json.loads(esummary_handle.readline().decode('utf-8'))

      # efetch_handle = m.ncbi.Entrez.efetch(
      #     {
      #         "db": "taxonomy",
      #         "retstart": retstart,
      #         "retmax": retmax,
      #         "webenv": webenv,
      #         "query_key": query_key
      #     }
      # )
      # for line in efetch_handle:
      #     cmd.write_to_log(line)

      i_taxon += len(esummary_json["result"]) - 1  # first result is just the submitted UIDs
  

  @staticmethod
  def get_taxon_id(query):

    esearch_json = m.ncbi.Entrez.esearch(
      esearch_params = {
        "db": "taxonomy",
        "term": query
      }
    )
    '''
    Example JSON returned:
    <DictionaryElement({
        'Count': '1',
        'RetMax': '1',
        'RetStart': '0',
        'IdList': ['102862'],
        'TranslationSet': [],
        'TranslationStack': [
            <DictionaryElement({
                'Term': 'Proteus penneri[All Names]',
                'Field': 'All Names',
                'Count': '1',
                'Explode': 'N',
            })>,
            'GROUP',
        ],
        'QueryTranslation': 'Proteus penneri[All Names]',
    })> (DictionaryElement) len=7
    '''

    taxon_id_list = esearch_json['IdList']

    if len(taxon_id_list) > 1:
      raise Exception(
         "Expected at most one taxon ID."
        f" query: {query}"
        f" taxon_id_list: {taxon_id_list}"
      )
    
    if len(taxon_id_list) < 1:
      raise Exception(
         "No taxon ID returned."
        f" query: {query}"
      )

    taxon_id = taxon_id_list[0]
    
    return taxon_id
