
import logging
import os

from devtools import debug
from django.db import models

import config
import modules as m
import util

# Configs


# class Command(models.Model):
class Command():

  # id = models.UUIDField(
  #   primary_key=True,
  #   default=uuid.uuid4,
  # )
  # name = models.CharField(max_length=255)
  #
  # class Meta:
  #   db_table = "command"

  # Constructor
  def __init__(self,
    cmd_fpath = None # typically: __file__
  ):
    # Set paths
    self.dpath = m.File.get_dpath(cmd_fpath)
    self.data_dpath = os.path.join(self.dpath, config.paths.dir_names.command_data)
    self.rscript_path = os.path.join(self.dpath, config.paths.default_filenames.r_script)
    self.log_fpath = os.path.join(self.data_dpath, "command.log")

    # Reset log file, for now...
    with open(self.log_fpath, 'w') as f:
      f.write("")

    # Create relative data directory if it doesn't already exist
    util.create_dpath_if_not_exist(self.data_dpath)


  def write_to_log(self, data):

    if not isinstance(data, str):
      data = str(data)

    if not data.endswith("\n"):
      data += "\n"

    with open(self.log_fpath, 'a') as f:
      f.write(data)


  def reset_log_file(self):
    with open(self.log_fpath, 'w') as f:
      f.write("")