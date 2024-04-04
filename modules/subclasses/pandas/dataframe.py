
import pandas as pd

import modules as m


class DataFrame():

    @staticmethod
    def create(dataframe_params):

        index = None
        if "index" in dataframe_params:
            index = dataframe_params["index"]

        columns = None
        if "columns" in dataframe_params:
            columns = dataframe_params["columns"]

        return pd.DataFrame(
            columns=columns
        )


    @classmethod
    def load_tsv(cls, fpath=None, read_csv_params={}):
        """
        https://pandas.pydata.org/docs/reference/api/pandas.read_csv.html
        :param fpath:
        :param read_csv_params:
        :return:
        """

        index_col = None
        if "index_col" in read_csv_params:
            index_col = read_csv_params["index_col"]

        usecols = None
        if "usecols" in read_csv_params:
            usecols = read_csv_params["usecols"]

        return pd.read_csv(
            fpath,
            header=0,
            index_col=index_col,
            usecols=usecols,
            sep="\t",
            dtype=str,
            parse_dates=False,
            keep_default_na=False,
            true_values=["Y"],
            false_values=["N"]
        )


    @staticmethod
    def sort_values(df, sort_values_params={}):
        """
        https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.sort_values.html
        :param sort_values_params:
        :return:
        """

        by = None
        if "by" in sort_values_params:
            by = sort_values_params["by"]

        ascending = True
        if "ascending" in sort_values_params:
            ascending = sort_values_params["ascending"]

        ignore_index = True
        if "ignore_index" in sort_values_params:
            ignore_index = sort_values_params["ignore_index"]

        return df.sort_values(
            by=by,
            ascending=ascending,
            ignore_index=ignore_index
        )