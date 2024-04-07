
import json
import os
import re

import modules as m


class Assembly():

    BASE_FNAME_SUFFIX = "-assembly"
    RAW_ESUMMARY_FNAME_SUFFIX = f"{BASE_FNAME_SUFFIX}-esummary-raw.tsv"
    NORMALIZED_ESUMMARY_FPATH_SUFFIX = f"{BASE_FNAME_SUFFIX}-esummary-normalized.tsv"
    ESUMMARY_ANALYTICS_FPATH_SUFFIX = f"{BASE_FNAME_SUFFIX}-esummary-metrics.tsv"


    @staticmethod
    def count(taxon_id):

        esearch_results = m.ncbi.Entrez.esearch({
            "db": "Assembly",
            "retmax": 0,
            "term": f"txid{taxon_id}[Organism]"
        })

        return esearch_results["Count"]


    @classmethod
    def get_esummary_batch(
            cls,
            taxon_id=None,
            clade_name=None,
            query=None,
            retstart=0,
            retmax=9999,
            retmode="json"
    ):

        if retmode == "json":
            retmax = 500

        # Special handling of queries
        term: str = None
        if query:
            term = query
        else:
            if taxon_id:
                term = f"txid{taxon_id}[Organism]"
            elif clade_name:
                term = f"{clade_name}[Organism]"

        # Get a list of the most recent records
        print(term)
        esearch_results = m.ncbi.Entrez.esearch({
            "db": "assembly",
            "retmax": retmax,
            "retstart": retstart,
            "term": term,
            "idtype": "gi"
        })
        # print(f"sl: esearch_results: {json.dumps(esearch_results, indent=2)}")

        # Exit function if there are no results
        if esearch_results["Count"] == "0":
            # print(f"sl: esearch_results[\"Count\"]: {esearch_results['Count']}")
            return None

        # print(f"sl: esearch_results B: {json.dumps(esearch_results, indent=2)}")

        # Get the Accession IDs from the results
        assembly_ids = esearch_results["IdList"]
        # print(f"sl: nuccore_ids: {nuccore_ids}")

        if len(assembly_ids) == 0:
            return "NO ID LIST RETURNED BY QUERY"

        # print(f"sl: biosample_ids: {json.dumps(biosample_ids, indent=2)}")

        # Upload a list of GI IDs to Entrez for the batch request
        epost_results = m.ncbi.Entrez.epost({
            "db": "assembly",
            "id": assembly_ids
        })

        # Get the identifiers of the batch request
        webenv = epost_results["WebEnv"]
        query_key = epost_results["QueryKey"]

        handle = m.ncbi.Entrez.esummary({
            "db": "assembly",
            "retmode": retmode,
            "webenv": webenv,
            "query_key": query_key
        })

        return {
            "handle": handle,
            "esearch_results": esearch_results,
            "nuccore_ids": assembly_ids
        }


    @staticmethod
    def get_fasta_handle(id_list, id_type="GI"):
        """

        :param id_list:
        :param id_type: "GI" | "ACCESSION" | "MIXED"
        :return:
        """

        for id in id_list:
            if not re.search("^\d+$", id) and id_type == "GI":
                raise Exception(f"This doesn't seem to be a GI: {id}")

        # Note: When using accession.version identifiers, there is a conversion step that takes place that causes large lists of
        # identifiers to time out, even when using POST. Therefore, we recommend batching these types of requests in sizes of
        # about 500 UIDs or less, to avoid retrieving only a partial amount of records from your original POST input list.
        # https://www.ncbi.nlm.nih.gov/books/NBK25499/
        if id_type != "GI" and len(id_list) > 500:
            raise Exception(f"Accession IDs batch size should not exceed 500, this batch size is: {len(id_list)}")

        # Get GIs from assembly accession IDs
        term = ""
        for id in id_list:
            term += id + "[ASAC] OR "
        term = term[:-4]

        esearch_results = m.ncbi.Entrez.esearch({
            "db": "assembly",
            "term": term,
            "idtype": "GI"
        })
        m.Log.write(f"esearch_results: {json.dumps(esearch_results, indent=2)}")
        gi_ids = esearch_results["IdList"]

        # Get Nucleotide GIs from Assembly GIs using ELink
        elink_json = m.ncbi.Entrez.elink(
            {
                "db": "nuccore",
                "dbfrom": "assembly",
                "cmd": "neighbor",
                "id": gi_ids,
                "retmode": "json",
                "idtype": "gi"
            },
            as_json=True
        )
        m.Log.write(f"elink_json: {json.dumps(elink_json, indent=2)}")
        # Just get the ones from GenBank (assembly_nuccore_insdc), because they should be identical to RefSeq
        quit()

        epost_results = m.ncbi.Entrez.epost({
            "db": "assembly",
            "id": gi_ids
        })

        # Get the identifiers of the batch request
        webenv = epost_results["WebEnv"]
        query_key = epost_results["QueryKey"]

        handle = m.ncbi.Entrez.efetch({
            "db": "assembly",
            "rettype": "docsum",
            "retmode": "xml",
            "webenv": webenv,
            "query_key": query_key
        })
        # print(f"sl: handle: {handle}")

        return handle