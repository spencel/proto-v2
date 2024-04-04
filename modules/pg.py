
import datetime
import logging
import os

from devtools import debug
import psycopg2

import config
import util


# Configs



class Pg():

    schema_dirpath = config.pg_schema_dpath

    @classmethod
    def get_cursor(cls, db_name=None):

        # Get default database name if a name wasn't passed as an argument
        if not db_name:
            db_name = config.db.default

        pg_con = psycopg2.connect(
            host=config.db[db_name].host,
            database=config.db[db_name].database,
            user=config.db[db_name].user,
            password=config.db[db_name].password,
            port=config.db[db_name].port
        )
        pg_con.set_session(autocommit=True)
        pg_cur = pg_con.cursor()

        return pg_con, pg_cur

    @classmethod
    def escape_single_quote(cls, arg_string):
        arg_string = arg_string.replace("'", "''")
        return arg_string

    @classmethod
    def clean_text_value(cls,
        value
    ):
        # Convert value to string if it's a list
        if isinstance(value, list):
            value = util.list_to_string(value)

        # Handle characters
        value = cls.escape_single_quote(value)

        return value

    @classmethod
    def create_db(cls,
        schemas = [],
        skip_schemas = []
    ):
        scope_pref = f"{cls.__name__}.create_db"

        pg_con, pg_cur = cls.get_cursor()

        for schema in schemas:
            # logging.debug(f"{scope_pref}: schema: {schema}")
            # Don't process any schemas in the skipped schemas list
            if schema in skip_schemas:
                continue

            schema_filepath = os.path.join(cls.schema_dirpath, f"{schema}.sql")
            with open(schema_filepath, 'r') as f:
                sql_qry = "".join(f.readlines())
                # logging.debug(f"{scope_pref}: executing {schema_filename}.sql")
                pg_cur.execute(sql_qry)

        pg_cur.close()
        pg_con.close()