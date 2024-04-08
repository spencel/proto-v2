
import os

import pandas as pd

import modules as m


class Genome():


    @classmethod
    def ingest_data(cls, fpath=None, taxon_id=None):

        in_df = m.Pandas.DataFrame.load_tsv(
            fpath,
            {
                "usecols": [
                    "Assembly Accession",
                    "Organism Infraspecific Names Strain",
                    "Assembly Submission Date"
                ]
            }
        )


        in_df = m.Pandas.DataFrame.sort_values(
            in_df,
            {
                "by": "Assembly Submission Date",
                "ascending": False
            }
        )

        max_samples = 100

        out_df = pd.DataFrame(
            columns=[
                "Assembly Accession",
                "Organism Infraspecific Names Strain",
                "Assembly Submission Date"
            ]
        )

        for i, row in in_df.iterrows():
            if i >= max_samples:
                break

            out_df = out_df._append({
                "Assembly Accession": row.at["Assembly Accession"],
                "Organism Infraspecific Names Strain": row.at["Organism Infraspecific Names Strain"],
                "Assembly Submission Date": row.at["Assembly Submission Date"]
            },
                ignore_index=True
            )


        out_fpath = os.path.join(
            m.file_sys.File.get_fpath_without_extension(fpath) +
            "-prepared.tsv"
        )

        out_df.to_csv(
            out_fpath,
            sep="\t",
            index=False
        )