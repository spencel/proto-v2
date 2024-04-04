
import json
import os
import re

import modules as m


class BioSample():


    @staticmethod
    def count(taxon_id):

        esearch_results = m.Entrez.esearch({
            "db": "BioSample",
            "retmax": 0,
            "term": f"txid{taxon_id}[Organism]"
        })

        return esearch_results["Count"]
