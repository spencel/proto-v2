
import datetime
import logging
import os
import sys

from devtools import debug
# import mariadb

import config
import util


# Configs



class MariaDb():

    schema_dirpath = config.pg_schema_dpath

    # @classmethod
    # def get_cursor(cls, db_name=None):
    #     # Connect to MariaDB Platform
    #     try:
    #         # This doesn't work
    #         conn = mariadb.connect(
    #             user="spencer",
    #             password="",
    #             host="manager",
    #             port=3306,
    #             database="onramp"
    #
    #         )
    #     except mariadb.Error as e:
    #         print(f"Error connecting to MariaDB Platform: {e}")
    #         sys.exit(1)
    #
    #     # Get Cursor
    #     return conn, conn.cursor()