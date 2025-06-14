
import json
import os
import xml.etree.ElementTree as ET

from modules import file_sys as m_file_sys
from modules import json as m_json

from .. import Entrez


class Taxonomy():


  @staticmethod
  def get_sub_taxon_ids(taxon_name, retstart=0):
    """
    Returns a separate set of all taxon IDs for each given parent taxon ID.
    :param taxon_ids:
    :return:
    """
    term = str(taxon_name) + "[Subtree]"

    esearch_json = Entrez.esearch(
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

      epost_json = Entrez.epost({
        "db": "Taxonomy",
        "id": esearch_json["IdList"]
      })
      webenv = epost_json["WebEnv"]
      query_key = epost_json["QueryKey"]

      esummary_retstart = 0
      while esummary_retstart < esearch_record_qty:
        esummary_json = Entrez.esummary({
          "db": "Taxonomy",
          "retstart": esummary_retstart,
          "webenv": webenv,
          "query_key": query_key
        }, as_json=True)
        esummary_record_qty = len(esummary_json["result"]["uids"])
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

    return species


  @staticmethod
  def get_esummary(taxon_ids):

    webenv, query_key = Entrez.get_epost_webenv_and_query_key(
      db="Taxonomy",
      id=taxon_ids
    )

    i_taxon = 0
    retstart = 0
    retmax = 500
    while i_taxon < len(taxon_ids):
      esummary_handle = Entrez.esummary({
        "db": "taxonomy",
        "retstart": retstart,
        "retmax": retmax,
        "webenv": webenv,
        "query_key": query_key
      })
      esummary_json = json.loads(esummary_handle.readline().decode('utf-8'))

      # efetch_handle = Entrez.efetch(
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

    esearch_json = Entrez.esearch(
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
  
  @staticmethod
  def get_taxon_ids(taxon_names=None, fpath=None):

    if fpath:
      taxon_names = set()
      with open(fpath) as f:
        for line in f:
          taxon_names.add(line.rstrip("\n"))
    
    taxon_ids = {}

    for taxon_name in taxon_names:
      taxon_id = Taxonomy.get_taxon_id(taxon_name)
      taxon_ids[taxon_name] = taxon_id

    return taxon_ids

  @staticmethod
  def export_lineages(
    taxon_ids_fpath: str|None = None,
    taxon_ids: set|None = None,
    col_idx: int = 0,
    has_header_row: bool = False,
    export_items: set = set(),
    out_fpath: str|None = None
  ) -> dict:
    """_summary_

    Args:
        taxon_ids_fpath (str): _description_
        col_idx (int, optional): Column index of the taxon IDs in the input. Defaults to 0.
        export_items (set, optional): _description_. Defaults to [].
        out_fpath (str | None, optional): _description_. Defaults to None.

    Returns:
        dict: Metrics...
          {
            'taxon_id_qty': (int)
            'unique_taxon_qty': (int)
            'unique_species_qty': (int)
          }
    """
    metrics = {
      'taxon_id_qty': 0,
      'unique_taxon_qty': 0,
      'unique_species_qty': 0
    }

    if not out_fpath:
      out_fpath = m_file_sys.File.get_fpath_without_extension(taxon_ids_fpath) + '-lineages.tsv'
    metrics_out_fpath = m_file_sys.File.get_fpath_without_extension(out_fpath) + '-metric.json'
    
    # Make list of unique taxon IDs
    if taxon_ids_fpath and not taxon_ids:
      taxon_ids = set()
      with open(taxon_ids_fpath, 'r') as f:
        
        # Skip column names row
        if has_header_row:
          for line in f:
            break
        
        for line in f:
          metrics['taxon_id_qty'] += 1
          taxon_id = line.rstrip('\n').split('\t')[col_idx]
          taxon_ids.add(taxon_id)

    # Post Taxon IDs
    epost_res = Entrez.epost({
      'db': 'taxonomy',
      'id': ','.join(taxon_ids)
    })

    # Get the identifiers of the batch request
    webenv = epost_res["WebEnv"]
    query_key = epost_res["QueryKey"]

    http_res_xml = Entrez.efetch({
      "db": "taxonomy",
      "webenv": webenv,
      "query_key": query_key
    }).read().decode('utf-8')

    tax_tree = dict()

    root = ET.fromstring(http_res_xml)
    taxons = root.findall('Taxon')
    for taxon in taxons:
      
      taxon_id = taxon.find('TaxId').text
      name = taxon.find('ScientificName').text
      rank = taxon.find('Rank').text

      lineage_taxons = taxon.findall('LineageEx/Taxon')
      previous_lineage_taxon_id = None
      for lineage_taxon in lineage_taxons:
        lineage_taxon_id = lineage_taxon.find('TaxId').text
        lineage_name = lineage_taxon.find('ScientificName').text
        lineage_rank = lineage_taxon.find('Rank').text

        if lineage_taxon_id not in tax_tree:
          metrics['unique_taxon_qty'] += 1

          tax_tree[lineage_taxon_id] = {
            'name': lineage_name,
            'rank': lineage_rank,
            'parent_id': previous_lineage_taxon_id
          }

        previous_lineage_taxon_id = lineage_taxon_id
      
      if taxon_id not in tax_tree:
        metrics['unique_species_qty'] += 1
        metrics['unique_taxon_qty'] += 1

        tax_tree[taxon_id] = {
          'name': name,
          'rank': rank,
          'parent_id': previous_lineage_taxon_id
        }

    table_lines = m_json.to_table(
      tax_tree,
      {'key_col_name': 'taxon_id'}
    )
    with open(out_fpath, 'w') as f:
      for line in table_lines:
        clean_line = [str(item) if item is not None else 'None' for item in line]
        f.write('\t'.join(clean_line) + '\n')
    
    m_json.save_to_file(metrics, metrics_out_fpath)

    # return metrics