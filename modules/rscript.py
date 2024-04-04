

import logging
import os
import subprocess


# Configs




class RScript():

    @classmethod
    def run(cls, r_script_path):
        scope_pref = f"{cls.__name__}.run"
        logging.info(f"{scope_pref}: R script: {r_script_path}")

        subprocess.run(
            ["Rscript", r_script_path],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.STDOUT
        )
        # subprocess.run(
        #     ["Rscript", r_script_path]
        # )