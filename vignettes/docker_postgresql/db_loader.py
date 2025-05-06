#!/usr/bin/python3

import argparse
import logging
import os

from psycopg import Psycopg2


logging.basicConfig(level=logging.DEBUG)


if __name__ == '__main__':

	Psycopg2.gen_db(
		new_db_name = os.getenv("POSTGRES_DB_NAME"),
		schema_fpath = os.getenv("DB_SCHEMA_FPATH"),
		tables_dpath = os.getenv("DB_TABLES_DPATH"),
		db_conn_config = {
			'dbname': 'postgres',
			'user': f'{os.getenv("POSTGRES_DB_USERNAME")}',
			'password': f'{os.getenv("POSTGRES_DB_PASWORD")}',
			'host': f'{os.getenv("POSTGRES_DB_HOST")}',
			'port': f'{os.getenv("POSTGRES_DB_PORT")}'
		},
		override_existing_db = False
	)