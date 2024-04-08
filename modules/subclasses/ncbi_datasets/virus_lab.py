
import json
import os

import pandas as pd

import modules as m


class VirusLab():

    @staticmethod
    def ingest_influenza_data(fpath=None, taxon_id=None):
        """
        The data can be found navigating from this page: https://www.ncbi.nlm.nih.gov/datasets/
        This is an example filter: https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/virus?SeqType_s=Nucleotide&VirusLineage_ss=Influenza%20A%20virus,%20taxid:11320&utm_source=data-hub
        :param fpath:
        :return:
        """

        # Get the data into a dataframe
        df = m.DataFrame.load_tsv(
            fpath,
            {
                "usecols": [
                    "Accession",
                    "Genotype",
                    "Segment",
                    "GenBank_Title",
                    "Collection_Date"
                ]
            }
        )

        # Sort so that the earliest samples are first and the chromosomes are grouped consecutively
        # The segment number in the GenBank title refers to which chromosome the record represents
        df = m.Pandas.DataFrame.sort_values(
            df,
            {
                "by": ["Collection_Date", "GenBank_Title"],
                "ascending": [False, False]
            }
        )

        # Going to group the chromosomes for each sample using a dictionary
        samples = dict()
        # {
        #   "Influenza A virus (A/Michigan/UOM10045650349/2022(H3N2))": [
        #       <Accession>,
        #       <Genotype>,
        #       <Segment>,
        #       <Collection_Date>
        #   ],
        #   ...
        # }

        max_samples: int
        if "partial-genome" in fpath:
            max_samples = 100
        elif "complete-genome" in fpath:
            max_samples = 3
        else:
            raise Exception("Error: add 'partial-genome' or 'complete-genome' to filename.")

        for i, row in df.iterrows():

            # It's difficult to parse the GenBank title, there may be more spaces than expected.
            # Use "segment" to separate the GenBank title into something more parsable.
            # Ignore any rows that have no segment in the title
            genbank_title = row.at["GenBank_Title"]
            index = genbank_title.find(" segment ")
            if index:
                modified_title = genbank_title[:index]
                if modified_title not in samples and len(samples.keys()) < max_samples:
                    samples[modified_title] = [[
                        row.at["Genotype"],
                        row.at["Segment"],
                        row.at["Accession"]
                    ]]
                elif modified_title in samples:
                    samples[modified_title].append([
                        row.at["Genotype"],
                        row.at["Segment"],
                        row.at["Accession"]
                    ])
        m.Log.write(json.dumps(samples, indent=2))

        # # Get the GIs for the selected records
        # Actually, just get FASTAs by Accession
        # accessions = []
        # for modified_title in samples:
        #     records = samples[modified_title]
        #     for record in records:
        #         accessions.append(record[2])
        # esearch_json = m.ncbi.Entrez.esearch({
        #     "db": "nuccore",
        #     # Term is a comma separated list of accessions
        #     "term": ",".join(accessions)
        # })

        # Write the samples to a file
        out_fpath = os.path.join(
            m.file_sys.File.get_fpath_without_extension(fpath) + "-prepared.tsv"
        )
        with open(out_fpath, 'w') as f:
            f.write("Modified Title\tGenotype\tSegment\tAccession\n")
            for modified_title in samples:
                for item in samples[modified_title]:
                    f.write(
                        f"{modified_title}\t"
                        f"{item[0]}\t"
                        f"{item[1]}\t"
                        f"{item[2]}\n"
                    )



    @classmethod
    def ingest_data(cls, fpath=None, taxon_id=None):

        # Influenza A & B
        if taxon_id in [11320, 11520]:
            cls.ingest_influenza_data(fpath, taxon_id)
