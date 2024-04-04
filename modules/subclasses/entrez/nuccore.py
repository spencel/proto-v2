
import datetime
import json
import os
import re

import modules as m


class NucCore():

    BASE_FNAME_SUFFIX = "-nuccore"
    RAW_ESUMMARY_FNAME_SUFFIX = f"{BASE_FNAME_SUFFIX}-esummary-raw.tsv"
    NORMALIZED_ESUMMARY_FPATH_SUFFIX = f"{BASE_FNAME_SUFFIX}-esummary-normalized.tsv"
    ESUMMARY_ANALYTICS_FPATH_SUFFIX = f"{BASE_FNAME_SUFFIX}-esummary-metrics.tsv"

    CURATED_DATA = {
        # <GI>: { "host": ..., ... }
        "1828053884": {
            "host": "Homo sapiens"
        },
        "1827015957": {
            "host": "Homo sapiens"
        }
    }


    @staticmethod
    def count(taxon_id=None, term=None):

        if not term:
            term = f"txid{taxon_id}[Organism]"

        esearch_results = m.Entrez.esearch({
            "db": "Nucleotide",
            "retmax": 0,
            "term": term
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
        retmode="json",
        time_frame=None,
        sort_by=None,
        get_complete_genome=True
    ):
        # print(f"sl: sort_by: {sort_by}")

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
                # print(f"sl: start_date: {start_date}")
                # print(f"sl: end_date: {end_date}")

            query_time_frame = str(
                f"(\"{start_date.strftime('%Y/%m/%d')}\"[PDAT] : \"{end_date.strftime('%Y/%m/%d')}\"[PDAT])"
            )
        if query_time_frame:
            term += " AND " + query_time_frame
        # print(f"sl: term: {term}")

        # Get a list of the most recent records
        esearch_results = m.Entrez.esearch({
            "db": "nuccore",
            "retmax": retmax,
            "retstart": retstart,
            "term": term,
            "idtype": "gi",
            "sort": sort_by
        })
        # print(f"sl: esearch_results: {json.dumps(esearch_results, indent=2)}")

        # Exit function if there are no results
        if esearch_results["Count"] == "0":
            # print(f"sl: esearch_results[\"Count\"]: {esearch_results['Count']}")
            return None

        # print(f"sl: esearch_results B: {json.dumps(esearch_results, indent=2)}")

        # Get the Accession IDs from the results
        nuccore_ids = esearch_results["IdList"]
        # print(f"sl: nuccore_ids: {nuccore_ids}")

        if len(nuccore_ids) == 0:
            return "NO ID LIST RETURNED BY QUERY"

        # print(f"sl: biosample_ids: {json.dumps(biosample_ids, indent=2)}")

        # Upload a list of GI IDs to Entrez for the batch request
        epost_results = m.Entrez.epost({
            "db": "nuccore",
            "id": nuccore_ids
        })

        # Get the identifiers of the batch request
        webenv = epost_results["WebEnv"]
        query_key = epost_results["QueryKey"]

        handle = m.Entrez.esummary({
            "db": "nuccore",
            "retmode": retmode,
            "webenv": webenv,
            "query_key": query_key
        })

        return {
            "handle": handle,
            "esearch_results": esearch_results,
            "nuccore_ids": nuccore_ids
        }


    @classmethod
    def prepare_esummary_data(
        cls,
        f_in_fpath,
        f_out_fpath,
        remove_input_files=False
    ):
        col_names = m.File.get_last_line(f_in_fpath).strip("\n")
        col_qty = len(col_names.split("\t"))

        # Get index of collection date column
        collection_date_col_idx: int = None
        for i, col_name in enumerate(col_names.split("\t")):
            if col_name == "collection_date":
                collection_date_col_idx = i

        with open(f_in_fpath, 'r') as f_in, \
                open(f_out_fpath, 'w') as f_out:

            f_out.write(f"{col_names}\tNormalized Collection Date\n")

            for line in f_in:

                # Clean line
                line = line.strip("\n")
                # Add "\t" to line if needed, because total columns wasn't known earlier
                row_length = len(line.split("\t"))
                while row_length < col_qty:
                    line += "\t"
                    row_length += 1

                if collection_date_col_idx:
                    collection_date = line.split("\t")[collection_date_col_idx]
                else:
                    collection_date = "0000-00-00"

                # Skip row that has column names
                if line.startswith("GI\tAccession ID\tTitle\tStrain\tCompleteness"):
                    continue

                if not re.search("^\d{4}-\d{2}-\d{2}$", collection_date):
                    try:
                        normalized_date = m.DataNormalizer.normalize_date(collection_date)
                    except:
                        raise Exception(
                            f"File: {f_in_fpath}\n"
                            f"Line: {line}"
                        )
                else:
                    normalized_date = collection_date

                f_out.write(f"{line}\t{normalized_date}\n")

        if remove_input_files:
            m.File.delete(f_in_fpath)

    @classmethod
    def prepare_esummary_data_from_taxon_ids(
            cls,
            taxon_ids: list[str] = None,
            dpath: str = None
    ):

        for taxon_id in taxon_ids:

            f_in_fpath = os.path.join(dpath, f"{taxon_id}{cls.RAW_ESUMMARY_FNAME_SUFFIX}")
            f_out_fpath = os.path.join(dpath, f"{taxon_id}{cls.NORMALIZED_ESUMMARY_FPATH_SUFFIX}")

            cls.prepare_esummary_data(f_in_fpath, f_out_fpath)


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

        epost_results = m.Entrez.epost({
            "db": "nuccore",
            "id": id_list
        })

        # Get the identifiers of the batch request
        webenv = epost_results["WebEnv"]
        query_key = epost_results["QueryKey"]

        handle = m.Entrez.efetch({
            "db": "nuccore",
            "rettype": "fasta",
            "retmode": "text",
            "webenv": webenv,
            "query_key": query_key
        })
        # print(f"sl: handle: {handle}")

        return handle