
import datetime
import json
import logging
import os
import subprocess

import Bio
import Bio.Entrez

import config as c
from ..assembly import Assembly
from ..biosample import BioSample
from ..nucleotide import Nucleotide
from ..taxonomy import Taxonomy

# Configs


Bio.Entrez.email = c.ncbi.entrez_email


class Entrez():

  Assembly = Assembly
  BioSample = BioSample
  Nucleotide = Nucleotide
  Taxonomy = Taxonomy

  # List of official species names, ie, only genus and species
  CORRECT_SPECIES_NAMES = {
    "anguilla japonica",
    "anser albifrons",
    "apis cerana",
    "apis mellifera",
    "aselliscus stoliczkanus",
    "bison bison",
    "camelus dromedarius",
    "chaerephon plicata",
    "chlorocebus aethiops",
    "chlorocebus sabaeus",
    "elephas maximus",
    "equus caballus",
    "gallus gallus",
    "homo sapiens",
    "macaca mulatta",
    "mandrillus sphinx",
    "meleagris gallopavo",
    "mus musculus",
    "nocardia brasiliensis",
    "rhinolophus affinis",
    "rhinolophus cornutus",
    "rhinolophus ferrumequinum",
    "rhinolophus pusillus",
    "rhinolophus sinicus",
    "rhinolophus sp.",
    "rhinolophus stheno",
    "seriola quinqueradiata",
    "vespertilio murinus",
    "zalophus californianus"
  }

  # List of names that aren't official species names, ie, something else besides genus and species
  SPECIAL_SPECIES_NAMES = {
    "bat": {},
    "beer": {},
    "bovine": {},
    "camel": {},
    "cow": {},
    "fish": {},
    "guinea pig": {},
    "homo sapiens 29 year old woman": {
      "species_name": "Homo sapiens"
    },
    "homo sapiens 40 year old man": {
      "species_name": "Homo sapiens"
    },
    "human (child)": {
      "species_name": "Homo sapiens"
    },
    "kefir grains": {},
    "mice": {},
    "mouse": {},
    "mus musculus subsp. domesticus": {
      "comment": "This is a subspecies."
    },
    "pasture gramineae": {},
    "pig": {
      "comment": "Check taxon ID 9822."
    },
    "plateau pika": {},
    "protozoa": {},
    "rat": {},
    "rhinilophus ferrumequinum": {
      "species_name": "rhinolophus ferrumequinum"
    },
    "swine": {}
  }

  """
  The BioPython Entrez Subpackage
  Source: https://biopython.org/docs/1.76/api/Bio.Entrez.html
  Bio.Entrez
    .efetch(db, **keywords)
    .esearch(db, term, **keywds)
    .elink(**keywds)
    .einfo(**keywds)
      Returns a summary of Entrez databases or if a database name is provided then it provides metadata on that data base.
    .esummary(**keywds)
    .egquery(**keywds)
    .espell(**keywds)
    .ecitmatch(**keywds)
    .read(handle, validate=True, escape=False)
    .parse(handle, validate=True, escape=False)
      Parses an XML file.
  """



### BioPython Entrez Wrapper Methods ########################################################################################

  @classmethod
  def einfo(cls, einfo_params=None, return_handle=False):
    """

    :param einfo_params: {
      "db": None | str, eg, "nuccore"
        If no database is provided, then returns meta data on available databases.
      "version": "2.0", only "2.0" is available
      "retmode": "xml" | "json"
    }
    :return:
    """

    # Handle einfo arguments

    # Defaults
    db = None
    version = "2.0"
    retmode = "json"

    einfo_handle = None
    if einfo_params:
      if "db" in einfo_params:
        db = einfo_params["db"]
      if "version" in einfo_params:
        version = einfo_params["version"]
      if "retmode" in einfo_params:
        retmode = einfo_params["retmode"]

      if db:
        # Get metadata on a database
        einfo_handle = Bio.Entrez.einfo(
          db=db,
          version=version,
          retmode=retmode
        )
      else:
        # Get a list of available databses
        einfo_handle = Bio.Entrez.einfo(
          retmode=retmode
        )

    # If no einfo_params were passed
    else:
      # Get a list of available databses
      einfo_handle = Bio.Entrez.einfo(
        retmode=retmode
      )

    if return_handle:
      return einfo_handle

    einfo_results = None
    if retmode == "json":
      einfo_results = json.load(einfo_handle)
    else:
      einfo_results = Bio.Entrez.read(einfo_handle)

    return einfo_results


    # einfo_results = Bio.Entrez.read(einfo_handle)


  @classmethod
  def esearch(cls, esearch_params=None, return_handle=False):
    """
    ESearch Documentation
      https://www.ncbi.nlm.nih.gov/books/NBK25499/#_chapter4_ESearch_

    Nucleotide/Nuccore, Protein, EST, GSS Field Discriptions
      https://www.ncbi.nlm.nih.gov/books/NBK49540/

    :param esearch_params:
      {
        "db": str, <database-name>, eg, "nuccore"
        "retmax": int
        "retstart": int
        "term": str, eg, "txid10515[Organism] AND (complete genome[Title])"
        "idtype": str, eg, "acc" - accession
        "sort": str, eg, "default"
      }
    :return: esearch_results = Bio.Entrez.read(esearch_response) = {
      "Count": str,
      "RetMax": str,
      "RetStart": str,
      "IdList": [
        str,
        ...
      ],
      "TranslationSet": [???],
      "TranslationStack": [
        {
          "Term": str, eg, "txid11137[Organism]",
          "Field": str, eg, "Organism",
          "Count": str, eg, "814",
          "Explode": str, eg, "Y"
        },
        {
          "Term": str, eg, "complete genome[Title]",
          "Field": str, eg, "Title",
          "Count": str, eg, "2173442",
          "Explode": str, eg, "N"
        },
        str, eg, "AND",
        ...?
      ],
      "QueryTranslation": str, eg, "txid11137[Organism] AND complete genome[Title]"
    }
    """

    # Handle parameters

    db = "nuccore"
    if "db" in esearch_params:
      db = esearch_params["db"]

    retmax = 9_999
    if "retmax" in esearch_params:
      retmax = esearch_params["retmax"]

    retstart = 0
    if "retstart" in esearch_params:
      retstart = esearch_params["retstart"]

    term = ""
    if "term" in esearch_params:
      term = esearch_params["term"]

    idtype = None
    if "idtype" in esearch_params:
      idtype = esearch_params["idtype"]

    esearch_handle = Bio.Entrez.esearch(
      db=db,
      retmax=retmax,
      retstart=retstart,
      term=term,
      idtype=idtype
    )

    if return_handle:
      return esearch_handle


    esearch_json = Bio.Entrez.read(esearch_handle)
    esearch_handle.close()

    return esearch_json
    # {
    #   "Count": str, # Count of results
    #   "RetMax": str,
    #   "RetStart": str,
    #   "IdList": [
    #       str,
    #       ...
    #   ]
    # }


  @classmethod
  def elink(cls, elink_params=None, as_json=True):
    """
    https://www.ncbi.nlm.nih.gov/books/NBK25499/#chapter4.ELink
    :param elink_params:
      {
        "db": str, Retrieves UIDs from this database.
        "dbfrom": str, Database containing the input UIDs.
        "cmd": str, The mode of the link.
          "neighbor" - (default)
          "neighbor_score" - similarity scores
          "neighbor_history" - returns query_key and WebEnv
          "acheck" - checks available elink lists
          "ncheck" - checks for links within the input database
          "lcheck" - checks for external links (LinkOuts)
          "llinks" - checks for LinkOuts that are not libraries
          "llinkslib" - check for LinkOUts including libraries
          "prlinks" - lists primary LinkOut provider

      }
    :return:
    """

    db = "nuccore"
    if "db" in elink_params:
      db = elink_params["db"]

    dbfrom = "taxonomy"
    if "dbfrom" in elink_params:
      dbfrom = elink_params["dbfrom"]

    cmd = "neighbor"
    #   neighbor (default)
    #   neighbor_score
    #   neighbor_history
    #   acheck
    #   ncheck
    #   lcheck
    #   llinks
    #   llinkslib
    #   prlinks
    if "cmd" in elink_params:
      cmd = elink_params["cmd"]

    id = None
    if "id" in elink_params:
      id = ",".join(elink_params["id"])

    query_key = None
    if "query_key" in elink_params:
      query_key = elink_params["query_key"]
    # print(f"sl: query_key: {query_key}")

    webenv = None
    if "webenv" in elink_params:
      webenv = elink_params["webenv"]
    # print(f"sl: webenv: {webenv}")

    retmode: str = "json"
    #   xml (default)
    #   json
    if "retmode" in elink_params:
      retmode = elink_params["retmode"]
    if retmode == "xml":
      as_json = False

    idtype: str = None
    #   GI (default)
    #   acc - accession ID
    if "idtype" in elink_params:
      idtype = elink_params["idtype"]

    linkname: str = None
    # Eg: taxonomy_nucleotide_exp
    if "linkname" in elink_params:
      linkname = elink_params["linkname"]

    term = None
    if "term" in elink_params:
      term = elink_params["term"]

    holding = None
    if "holding" in elink_params:
      holding = elink_params["holding"]

    elink_handle = Bio.Entrez.elink(
      db=db,
      dbfrom=dbfrom,
      cmd=cmd,
      id=id,
      query_key=query_key,
      webenv=webenv,
      retmode=retmode,
      idtype=idtype,
      linkname=linkname,
      term=term,
      holding=holding
    )

    if not as_json:
      return elink_handle

    elink_json = json.load(elink_handle)

    return elink_json


  @classmethod
  def epost(cls, epost_params=None, return_handle=False):
    """
    https://www.ncbi.nlm.nih.gov/books/NBK25499/#_chapter4_EPost_
    :param epost_params: {
      "db": str, <database-name>, eg, "nuccore"
      "id": list, eg, [0, 1, ...], the method will covert the list to a string.
    }
    :param return_handle:
    :return:
    """

    db = "nuccore"
    if "db" in epost_params:
      db = epost_params["db"]

    # For sequence databases (nuccore, popset, protein), the UID list may be a mixed list of GI numbers and
    # accession.version identifiers. Note: When using accession.version identifiers, there is a conversion step that
    # takes place that causes large lists of identifiers to time out, even when using POST. Therefore, we recommend
    # batching these types of requests in sizes of about 500 UIDs or less, to avoid retrieving only a partial amount
    # of records from your original POST input list.
    id = "0"
    if "id" in epost_params:
      id = epost_params["id"]
      if not isinstance(id[0], str):
        id = ",".join(map(str, id))
      else:
        id = ",".join(id)

    epost_handle = None
    try:
      epost_handle = Bio.Entrez.epost(
        db=db,
        id=id
      )
    except Exception as e:
      with open("log.err", 'w') as f:
        f.write(f"Bio.Entrez.epost(): id = \"{id}\"")
      raise e

    if return_handle:
      return epost_handle

    try:
      epost_json = Bio.Entrez.read(epost_handle)
      epost_handle.close()

      return epost_json

    except Exception as e:
      with open("log.err", 'w') as f:
        f.write(f"Bio.Entrez.read(): epost params: id = \"{id}\"")
      raise e


  @classmethod
  def efetch(cls, efetch_params=None):
    """
    https://www.ncbi.nlm.nih.gov/books/NBK25499/#_chapter4_EFetch_
    :param efetch_params: {
      "db": str, <database-name>, eg, "nuccore"
      "rettype": str,
      "retmode": str,
      "webenv": webenv
      "query_key": query_key
    }

    rettype and retmode combinations:
    source: https://www.ncbi.nlm.nih.gov/books/NBK25499/table/chapter4.T._valid_values_of__retmode_and/
    <rettype>
      <retmode>
      ...
    "docum"
      "xml" (default) - SEQUENCES
      "text" - SEQUENCES
    null - binary ASN.1 or text ASN.1
      "text (default) - SEQUENCES - seems to refer to text
      "asn.1" - seems to refer to binary
    "native" - full record in XML
      "xml (default) - SEQUENCES
      "text" - SEQUENCES
    "acc" - accession number(s)
      "text"
    "fasta" - FASTAs - SEQUENCES
      "text" - SEQUENCES
    "seqid" - SeqID string
      "text"
    "gb"
      "txt" - GenBank flat file
      "xml" - GBSeq XML
    "gbc"
      "xml" - INSDSeq XML
    "ft" - feature table
      "text"
    "gbwithparts" - GenBank flat file with full sequence (contigs)
      "text"
    "fasta_cds_na" - CDS nucleotide FASTA
      "text"
    "fasta_cds_aa" - CDS protein FASTA
      "text"

    :return:
    """

    db = "nuccore"
    if "db" in efetch_params:
      db = efetch_params["db"]

    id = None
    if "id" in efetch_params:
      id = efetch_params["id"]

    rettype = "docsum"
    if "rettype" in efetch_params:
      rettype = efetch_params["rettype"]

    retmode = "xml"
    if "retmode" in efetch_params:
      retmode = efetch_params["retmode"]

    webenv = None
    if "webenv" in efetch_params:
      webenv = efetch_params["webenv"]

    query_key = None
    if "query_key" in efetch_params:
      query_key = efetch_params["query_key"]

    # Get the records identified by the batch request
    efetch_handle = Bio.Entrez.efetch(
      db=db,
      id=id,
      rettype=rettype,
      retmode=retmode,
      webenv=webenv,
      query_key=query_key
    )

    return efetch_handle


  @classmethod
  def esummary(cls, esummary_params=None, as_json=False):
    """If set to xml, not as much information is returned when compared to a json request.

    Args:
        esummary_params (_type_, optional): _description_. Defaults to None.
        as_json (bool, optional): _description_. Defaults to False.

    Returns:
        _type_: _description_
        
        Example JSON response from Bio.Entrez.esummary():
        "2125450602": {
            "uid": "2125450602",
            "term": "2125450602",
            "caption": "OK073084",
            "title": "Human coronavirus 229E strain P42, complete genome",
            "extra": "gi|2125450602|gb|OK073084.1|",
            "gi": 2125450602,
            "createdate": "2022/10/13",
            "updatedate": "2022/10/13",
            "flags": "",
            "taxid": 11137,
            "slen": 27239,
            "biomol": "genomic",
            "moltype": "rna",
            "topology": "linear",
            "sourcedb": "insd",
            "segsetsize": "",
            "projectid": "0",
            "genome": "genomic",
            "subtype": "strain|host|country|isolation_source|collection_date",
            "subname": "P42|Homo sapiens|China|nasopharyngeal aspirate|2017",
            "assemblygi": "",
            "assemblyacc": "",
            "tech": "",
            "completeness": "complete",
            "geneticcode": "1",
            "strand": "",
            "organism": "Human coronavirus 229E",
            "strain": "P42",
            "biosample": "",
            "statistics": [
            {
                "type": "Length",
                "count": 27239
            },
            {
                "type": "Length",
                "subtype": "literal",
                "count": 27239
            },
            {
                "type": "all",
                "count": 15
            },
            {
                "type": "blob_size",
                "count": 217441
            },
            {
                "type": "cdregion",
                "count": 7
            },
            {
                "type": "cdregion",
                "subtype": "CDS",
                "count": 7
            },
            {
                "type": "gene",
                "count": 7
            },
            {
                "type": "gene",
                "subtype": "Gene",
                "count": 7
            },
            {
                "type": "org",
                "count": 1
            },
            {
                "type": "pub",
                "count": 2
            },
            {
                "type": "pub",
                "subtype": "unpublished",
                "count": 1
            },
            {
                "source": "CDS",
                "type": "all",
                "count": 7
            },
            {
                "source": "CDS",
                "type": "prot",
                "count": 7
            },
            {
                "source": "CDS",
                "type": "prot",
                "subtype": "Prot",
                "count": 7
            },
            {
                "source": "all",
                "type": "Length",
                "count": 27239
            },
            {
                "source": "all",
                "type": "all",
                "count": 22
            },
            {
                "source": "all",
                "type": "blob_size",
                "count": 217441
            },
            {
                "source": "all",
                "type": "cdregion",
                "count": 7
            },
            {
                "source": "all",
                "type": "gene",
                "count": 7
            },
            {
                "source": "all",
                "type": "org",
                "count": 1
            },
            {
                "source": "all",
                "type": "prot",
                "count": 7
            },
            {
                "source": "all",
                "type": "pub",
                "count": 2
            }
            ],
            "properties": {
            "na": "1",
            "value": "1"
            },
            "oslt": {
            "indexed": true,
            "value": "OK073084.1"
            },
            "accessionversion": "OK073084.1"
        }
      
      Example XML response from Bio.Entrez.esummary():
        <DocSum>
        <Id>2125450602</Id>
        <Item Name="Caption" Type="String">OK073084</Item>
        <Item Name="Title" Type="String">Human coronavirus 229E strain P42, complete genome</Item>
        <Item Name="Extra" Type="String">gi|2125450602|gb|OK073084.1|[2125450602]</Item>
        <Item Name="Gi" Type="Integer">2125450602</Item>
        <Item Name="CreateDate" Type="String">2022/10/13</Item>
        <Item Name="UpdateDate" Type="String">2022/10/13</Item>
        <Item Name="Flags" Type="Integer">0</Item>
        <Item Name="TaxId" Type="Integer">11137</Item>
        <Item Name="Length" Type="Integer">27239</Item>
        <Item Name="Status" Type="String">live</Item>
        <Item Name="ReplacedBy" Type="String"/>
        <Item Name="Comment" Type="String"> </Item>
        <Item Name="AccessionVersion" Type="String">OK073084.1</Item>
        </DocSum>
    """

    db = "nuccore"
    if "db" in esummary_params:
      db = esummary_params["db"]

    id = None
    if "id" in esummary_params:
      id = ",".join(esummary_params["id"])

    retstart = 0
    if "retstart" in esummary_params:
      retstart = esummary_params["retstart"]

    retmode = "json"
    if "retmode" in esummary_params:
      retmode = esummary_params["retmode"]

    retmax = 9999
    if "retmax" in esummary_params:
      retmax = esummary_params["retmax"]
    # The max for json format is 500 records
    if retmode == "json":
      retmax = 500

    webenv = None
    if "webenv" in esummary_params:
      webenv = esummary_params["webenv"]

    query_key = None
    if "query_key" in esummary_params:
      query_key = esummary_params["query_key"]

    esummary_handle = Bio.Entrez.esummary(
      db=db,
      id=id,
      retstart=retstart,
      retmax=retmax,
      retmode=retmode,
      webenv=webenv,
      query_key=query_key
    )

    if not as_json:
      return esummary_handle

    esummary_json = json.load(esummary_handle)

    return esummary_json

### END BioPython Entrez Wrapper Methods ####################################################################################


  @staticmethod
  def get_taxon_id_from_fpath(fpath) -> str:
    return os.path.basename(fpath).split("-")[0]


  @classmethod
  def get_einfo(cls, einfo_params=None):

    return cls.einfo(einfo_params)


  @classmethod
  def get_genbank_records(cls, taxon_id, at_most_record_qty=None, cut_off_date=None):
    """
    Retrieves at most of a given quantity of records for a given taxon ID. The GenBank files can be queried by taxon ID
    using the following pattern: txid<taxon-id>[Organism]. As an example, for the taxon ID 11137, the search term is
    txid11137[Organism].

    :param taxon_id:
    :param at_most_record_qty:
    :param cut_off_date: Keep in mind that available samples may be before the cut off date but will be included anyway
    if that is the case.
    """

    if not at_most_record_qty:
      at_most_record_qty = 20

    if not cut_off_date:
      cut_off_date = datetime.datetime(1900, 1, 1)

    # Special handling of queries
    term = f"txid{taxon_id}[Organism]"
    match taxon_id:

      # Influenza A & B	211044	Influenza A virus (A/Puerto Rico/8/1934(H1N1))	NC_002016.1	No complete genome, had records such as 'complete cds'.
      case 211044:
        pass

      # Influenza A & B	518987	Influenza B virus (B/Lee/1940)	NC_002205.1	No complete genome, had records such as 'complete cds'.
      case 518987:
        pass

      case _:
        term += " AND (complete genome[Title])"

    # Get a list of the most recent records
    esearch_results = cls.esearch({
      "db": "nucleotide",
      "retmax": at_most_record_qty,
      "retstart": 0,
      "term": term,
      "idtype": "gi"
    })

    # print(json.dumps(esearch_results, indent=2))

    # Get the Accession IDs from the results
    accession_ids = esearch_results["IdList"]

    if len(accession_ids) == 0:
      return "NO ID LIST RETURNED BY QUERY"

    # print(json.dumps(accession_ids, indent=2))

    # Upload a list of GI IDs to Entrez for the batch request
    try:
      results = cls.epost({
        "db": "nucleotide",
        "id": accession_ids
      })
    except:
      with open("debug.log", 'a') as f:
        f.write(json.dumps(esearch_results, indent=2)+"\n")
      raise Exception


    # Get the identifiers of the batch request
    webenv = results["WebEnv"]
    query_key = results["QueryKey"]

    # Get the records identified by the batch request
    handle = cls.efetch({
      "db": "nucleotide",
      # "rettype": "gb", # this might exclude sequence as opposed to gbwithparts
      "rettype": "gbwithparts",
      "retmode": "text",
      "webenv": webenv,
      "query_key": query_key
    })

    return {
      "handle": handle,
      "esearch_results": esearch_results,
      "accession_ids": accession_ids
    }


  @classmethod
  def get_biosample(
      cls,
      taxon_id,
      at_most_record_qty=None,
      cut_off_date=None,
      esearch_retmax=9_999,
      esearch_retstart=0
  ):
    """
    Retrieves at most of a given quantity of records for a given taxon ID. The GenBank files can be queried by taxon ID
    using the following pattern: txid<taxon-id>[Organism]. As an example, for the taxon ID 11137, the search term is
    txid11137[Organism].

    :param taxon_id:
    :param at_most_record_qty:
    :param cut_off_date: Keep in mind that available samples may be before the cut off date but will be included anyway
    if that is the case.
    """

    # print("sl: m.ncbi.Entrez.get_biosample: start")

    if not at_most_record_qty:
      at_most_record_qty = 100

    # Default is one year ago
    if not cut_off_date:
      cut_off_date = m.Time.get_one_year_ago()

    # Query building and special handling of queries
    term = f"txid{taxon_id}[Organism]"

    # Get a list of the most recent records
    esearch_results = cls.esearch({
      "db": "biosample",
      "retmax": esearch_retmax,
      "retstart": esearch_retstart,
      "term": term,
      "idtype": "bios"
    })

    # Exit function if there are no results
    if esearch_results["Count"] == "0":
      # print(f"sl: esearch_results[\"Count\"]: {esearch_results['Count']}")
      return None

    # print(f"sl: esearch_results B: {json.dumps(esearch_results, indent=2)}")

    # Get the Accession IDs from the results
    biosample_ids = esearch_results["IdList"]

    if len(biosample_ids) == 0:
      return "NO ID LIST RETURNED BY QUERY"

    # print(f"sl: biosample_ids: {json.dumps(biosample_ids, indent=2)}")

    # Upload a list of GI IDs to Entrez for the batch request
    try:
      results = cls.epost({
        "db": "biosample",
        "id": biosample_ids
      })
    except:
      with open("debug.log", 'a') as f:
        f.write(json.dumps(esearch_results, indent=2) + "\n")
      raise Exception

    # Get the identifiers of the batch request
    webenv = results["WebEnv"]
    # print(f"sl: webenv: {webenv}")
    query_key = results["QueryKey"]
    # print(f"sl: query_key: {query_key}")

    # Get the records identified by the batch request
    handle = cls.efetch({
      "db": "biosample",
      # "rettype":"gb", # this might exclude sequence as opposed to gbwithparts
      # "rettype":"docsum",
      "rettype": "full",
      # "retmode":"text",
      "webenv": webenv,
      "query_key": query_key
    })
    # print(f"sl: handle: {json.dumps(str(handle), indent=2)}")

    return {
      "handle": handle,
      "esearch_results": esearch_results,
      "biosample_ids": biosample_ids
    }


# OLD #######################################################################################################################

  @classmethod
  def get_esummary_old2(cls, taxon_id, at_most_record_qty=None, cut_off_date=None):
    """
    Retrieves at most of a given quantity of records for a given taxon ID. The GenBank files can be queried by taxon ID
    using the following pattern: txid<taxon-id>[Organism]. As an example, for the taxon ID 11137, the search term is
    txid11137[Organism].

    :param taxon_id:
    :param at_most_record_qty:
    :param cut_off_date: Keep in mind that available samples may be before the cut off date but will be included anyway
    if that is the case.
    """

    if not at_most_record_qty:
      at_most_record_qty = 100

    # Default is one year ago
    if not cut_off_date:
      cut_off_date = m.Time.get_one_year_ago()

    # Query building and special handling of queries
    term = f"txid{taxon_id}[Organism]"
    match taxon_id:

      # Influenza A & B	211044	Influenza A virus (A/Puerto Rico/8/1934(H1N1))	NC_002016.1	No complete genome, had records such as 'complete cds'.
      case 211044:
        pass

      # Influenza A & B	518987	Influenza B virus (B/Lee/1940)	NC_002205.1	No complete genome, had records such as 'complete cds'.
      case 518987:
        pass

      case _:
        term += " AND (complete genome[Title])"

    # Get a list of the most recent records
    esearch_results = cls.esearch({
      "db": "nucleotide",
      "retmax": 10_000,
      "retstart": 0,
      "term": term,
      "idtype": "gi"
    })
    # print(f"sl: esearch_results A: {esearch_results}")

    # print(json.dumps(esearch_results, indent=2))

    # Get the Accession IDs from the results
    accession_ids = esearch_results["IdList"]

    if len(accession_ids) == 0:
      return "NO ID LIST RETURNED BY QUERY"

    # print(json.dumps(accession_ids, indent=2))

    # Upload a list of GI IDs to Entrez for the batch request
    try:
      results = cls.epost({
        "db": "nucleotide",
        "id": accession_ids
      })
    except:
      with open("debug.log", 'a') as f:
        f.write(json.dumps(esearch_results, indent=2) + "\n")
      raise Exception

    # Get the identifiers of the batch request
    webenv = results["WebEnv"]
    query_key = results["QueryKey"]

    # Get the records identified by the batch request
    handle = cls.efetch({
      "db": "nucleotide",
      # "rettype"="gb", # this might exclude sequence as opposed to gbwithparts
      "rettype": "gbwithparts",
      "retmode": "text",
      "webenv": webenv,
      "query_key": query_key
    })

    return {
      "handle": handle,
      "esearch_results": esearch_results,
      "accession_ids": accession_ids
    }


  @classmethod
  def get_esummary_old1(cls, taxon_id, at_most_record_qty=None, cut_off_date=None):

    if not at_most_record_qty:
      at_most_record_qty = 20

    if not cut_off_date:
      cut_off_date = datetime.datetime(1900, 1, 1)

    term = f"txid{taxon_id}[Organism]"

    # Get a list of the most recent records
    esearch_results = cls.esearch({
      "db": "nucleotide",
      "retmax": at_most_record_qty,
      "retstart": 0,
      "term": term,
      "idtype": "gi"
    })

    # Get the Accession IDs from the results
    accession_ids = esearch_results["IdList"]

    if len(accession_ids) == 0:
      return "NO ID LIST RETURNED BY QUERY"

    # print(json.dumps(accession_ids, indent=2))

    # Upload a list of GI IDs to Entrez for the batch request
    results = cls.epost({
      "db": "nucleotide",
      "id": accession_ids
    })

    # Get the identifiers of the batch request
    webenv = results["WebEnv"]
    query_key = results["QueryKey"]

    # Get the records identified by the batch request
    handle = cls.esummary({
      "db": "nucleotide",
      # "rettype": "gb", # this might exclude sequence as opposed to gbwithparts
      "retmode": "xml",
      "webenv": webenv,
      "query_key": query_key
    })

    return {
      "handle": handle,
      "esearch_results": esearch_results,
      "accession_ids": accession_ids
    }


  @staticmethod
  def retrieve_fastas_old(genomes):
    genbank_prefixes = ("AF", "AP", "CP", "MN", "NC_", "NZ_", "OM", "ON", "OP")
    assembly_prefixes = ("GCA_",)

    # Download the fasta or assembly if it hasn't yet
    # Todo: this doesn't always work, the E-utilities API is convoluted
    for genome in genomes:

      if not genome.startswith(genbank_prefixes + assembly_prefixes):
        raise Exception(f"Unexpected genome ID encountered: {genome}")

      dpath = os.path.join(c.paths.genomes, genome)
      m.Dir.make_if_not_exist(dpath)

      # Download fasta
      if genome.startswith(genbank_prefixes):
        fpath = os.path.join(dpath, ".fa")
        if not os.path.exists(fpath):
          logging.info(f"Downloading {genome} fasta...")
          subprocess.run(
            f"esearch -db nucleotide -query \"{genome}\" | efetch -format fasta > {fpath}",
            shell=True
          )

      # Download assembly
      elif genome.startswith(assembly_prefixes):
        fpath = os.path.join(dpath, ".fna")
        if not os.path.exists(fpath):
          logging.info(f"Downloading {genome} assembly...")
          subprocess.run(
            f"esearch -db assembly -query {genome} | elink -target nucleotide -name assembly_nuccore_insdc | efetch -format fasta > {fpath}",
            shell=True
          )
      # Ensure defline is correct
      m.Fasta.clean(genome, fpath)
      # Add fasta file path to genomes data
      genomes[genome]["fasta_fpath"] = fpath

    return genomes


  @classmethod
  def get_epost_webenv_and_query_key(cls, db, id):
    epost_results = cls.epost({
      "db": db,
      "id": id,
      "return_handle": False
    })
    webenv = epost_results["WebEnv"]
    query_key = epost_results["QueryKey"]
    return webenv, query_key

  @classmethod
  def get_elink_webenv_and_query_key(cls, elink_params):
    elink_params["cmd"] = "neighbor_history"
    elink_json = cls.elink(elink_params, as_json=True)
    # {
    #   "header": {
    #     "type": "elink",
    #     "version": "0.3"
    #   },
    #   "linksets": [
    #     {
    #       "dbfrom": "taxonomy",
    #       "ids": [
    #         "2776901",
    #         "1383066",
    #         "1383065",
    #         "1383064",
    #         "1383063",
    #         "1382367",
    #         "1370127",
    #         "1312904",
    #         "1236535"
    #       ],
    #       "linksetdbhistories": [
    #         {
    #           "dbto": "nuccore",
    #           "linkname": "taxonomy_nucleotide_exp",
    #           "querykey": "1"
    #         }
    #       ],
    #       "webenv": "MCID_641d468a11143c4a080f4f54"
    #     }
    #   ]
    # }
    webenv = elink_json["linksets"][0]["webenv"]
    query_key = elink_json["linksets"][0]["linksetdbhistories"][0]["querykey"]
    return webenv, query_key