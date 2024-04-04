
import json
import logging
import os


from modules.subclasses.ncbi_datasets.virus_lab import VirusLab as VirusLabSubclass
from modules.subclasses.ncbi_datasets.genome import Genome as GenomeSubclass


# Configs




class NcbiDatasets():

    VirusLab = VirusLabSubclass
    Genome = GenomeSubclass

    @staticmethod
    def test(data_dpath=None):

        out_fpath = os.path.join(
            data_dpath, "test.out"
        )
        collection_dates_fpath = os.path.join(
            data_dpath, "collection-dates.out"
        )
        host_names_fpath = os.path.join(
            data_dpath, "host-names.out"
        )

        # with open(out_fpath, 'w') as f:
        #     subprocess.run([
        #         "datasets",
        #         "summary", "virus", "genome",
        #         "taxon", "1335626",
        #         # "--host", "Homo sapiens"
        #         "--as-json-lines",
        #         "--complete-only",
        #         "--api-key", "85ef7e99b4fbac1481e6d3ba8eb2f9130c09"
        #     ], stdout=f)
        with_host_qty = 0
        with_collect_date_qty = 0
        try:
            with open(out_fpath, 'r') as f_in, \
            open(collection_dates_fpath, 'w') as f_dates, \
            open(host_names_fpath, 'w') as f_names:
                for line in f_in:
                    line_json = json.loads(line)
                    # Data structure:
                    # ...
                    # "host": {
                    #   "lineage": [...],
                    #   "organism_name":"Homo sapiens",
                    #   "tax_id":9606
                    # },
                    # ...
                    # "isolate": {
                    #   "collection_date":"2019-04-03",
                    #   ...
                    # },
                    # ...
                    if "host" in line_json:
                        f_names.write(line_json["host"]["organism_name"] + "\n")
                        with_host_qty += 1
                    if "isolate" in line_json:
                        if "collection_date" in line_json["isolate"]:
                            f_dates.write(line_json["isolate"]["collection_date"] + "\n")
                            with_collect_date_qty += 1
        except:
            with open("err.log", 'w') as f:
                f.write(json.dumps(line_json, indent=2))
                raise Exception

        print(f"sl: with_host_qty: {with_host_qty}")
        print(f"sl: with_collect_date_qty: {with_collect_date_qty}")

        # # Retrieve and print genomic metadata for assemblies belonging to the specified taxon
        # for assembly in ncbi_datasets_metadata_genome.get_assembly_metadata_by_taxon("11137"):
        #     # ncbi_datasets_metadata_genome.print_assembly_metadata_by_fields(
        #     #     assembly, ["assembly_accession", "assembly_level", "seq_length"]
        #     # )
        #     ncbi_datasets_metadata_genome.print_assembly_metadata_by_fields(assembly, ["host", "collection_date"])

