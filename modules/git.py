
import os
import subprocess

import modules as m


class Git():

    @staticmethod
    def clone(
        host=None,
        account=None,
        repo_names: list[str] = None,
        out_dpath=None,
        overwrite=False
    ):

        for repo_name in repo_names:

            repo_dpath = os.path.join(out_dpath, repo_name)

            if overwrite:
                m.Dir.delete(repo_dpath)

            url = f"{host}/{account}/{repo_name}.git"

            original_working_directory = os.getcwd()
            print(f"original_working_directory: {original_working_directory}")

            os.chdir(out_dpath)
            print(f"os.getcwd() A: {os.getcwd()}")

            subprocess.run([
                "git", "clone",
                url,
                repo_name
            ])

            os.chdir(original_working_directory)
            print(f"os.getcwd() B: {os.getcwd()}")