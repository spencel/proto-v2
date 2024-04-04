

import subprocess


class VarAaSeg():


    @classmethod
    def create_qry_file(filepath, fa_defline, sequence):
        with open(filepath, 'w') as lastz_query_file:
            lastz_query_file.write(f'>{fa_defline}\n')
            lastz_query_file.write(sequence)
