
import logging
import os
import uuid

from devtools import debug
from django.db import models

import config
import modules as m
import util

# Configs





# class Command(models.Model):
class Gisaid():

    # id = models.UUIDField(
    #     primary_key=True,
    #     default=uuid.uuid4,
    # )
    # name = models.CharField(max_length=255)
    #
    # class Meta:
    #     db_table = "command"
    this_dpath = m.file_sys.File.get_dpath(__file__)
    data_dpath = os.path.join(
        config.paths.dir_names.data,
        config.paths.dir_names.gisaid_data
    )
    data_fpath = os.path.join(
        data_dpath,
        config.paths.default_filenames.gisaid_data
    )