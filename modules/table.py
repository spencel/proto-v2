import json

import modules as m


class Table():


    def __init__(self, col_names, rows=[]):
        """

        :param col_names: a dictionary with column names as keys and indices as values.
        :param rows:
        """
        self.col_names = col_names
        self.rows = rows

    @classmethod
    def load_file(
        cls,
        fpath,
        delimiter="\t",
        skip_col_names=False
    ):

        col_names = dict()
        rows = []
        with open(fpath, 'r') as f:

            if not skip_col_names:
                col_idx = 0
                for col_name in f.readline().strip("\n").split(delimiter):
                    if col_name in col_names:
                        raise Exception(f"Unexpected behavior: duplicate column names: {col_name}")
                    col_names[col_name] = col_idx
                    col_idx += 1

            for line in f:
                rows.append(line.strip("\n").split(delimiter))

        return cls(col_names, rows)


    def write_to_file(
        self,
        fpath,
        delimeter="\t",
        skip_col_names=False
    ):
        with open(fpath, 'w') as f:
            f.write(
                f"{delimeter}".join(self.col_names.keys()) + "\n"
            )
            for row in self.rows:
                for i, item in enumerate(row):
                    if not item:
                        row[i] = ""
                    elif not isinstance(item, str):
                        row[i] = str(item)
                f.write(f"{delimeter}".join(row) + "\n")


    def add_col(self, col_name):
        col_idx = len(self.col_names.keys())
        self.col_names[col_name] = col_idx
        for row in self.rows:
            row.append("")


    def add_row(self, data):
        """

        :param data: a list or dictionary with column names as keys with values...
        :return:
        """
        if isinstance(data, dict):
            self.rows.append([None] * len(self.col_names))
            i_row = len(self.rows) - 1
            for col_name in data:
                i_col = self.col_names[col_name]
                self.rows[i_row][i_col] = data[col_name]
        elif isinstance(data, list):
            self.rows.append(data)


    def sort(self, col_name=None, col_idx=None, in_place=False, order="ASCENDING"):
        """
        Sort by column name if provided or use column index.
        :param table:
        :param col_name:
        :param col_idx:
        :return:
        """
        if col_name:
            col_idx = self.col_names[col_name]

        if in_place:
            if order == "ASCENDING":
                self.rows.sort(key=lambda x: x[col_idx])
            else:
                self.rows.sort(key=lambda x: x[col_idx], reverse=True)

        else:
            if order == "ASCENDING":
                self.rows = sorted(self.rows, key=lambda x: x[col_idx])
            else:
                self.rows = sorted(self.rows, key=lambda x: x[col_idx], reverse=True)
