
from .subclasses import *


class File():

  # Subclasses
  Text = Text
  

import gzip
import os
from pathlib import Path


class File():


  Text = Text


  @staticmethod
  def get_size(fpath, units):
    size = os.path.getsize(fpath) # Bytes
    match units:
      case 'KB':
        size = size / (1000 ** 1)
      case 'MB':
        size = size / (1000 ** 2)
      case 'GB':
        size = size / (1000 ** 3)
      case 'TB':
        size = size / (1000 ** 4)
    
    return size

  @staticmethod
  def delete(fpath):
    if os.path.isfile(fpath):
      os.remove(fpath)

  @staticmethod
  def get_extension(fpath) -> str:
    """
    Removes the "." as well.
    :param fpath:
    :return:
    """
    return os.path.splitext(fpath)[1][1:]


  @staticmethod
  def get_fname(fpath) -> str:
    return os.path.basename(fpath)


  @staticmethod
  def get_fname_without_extension(fpath) -> str:
    return os.path.basename(
      os.path.splitext(fpath)[0]
    )


  @staticmethod
  def get_fpath_without_extension(fpath) -> str:
    return os.path.splitext(fpath)[0]


  @staticmethod
  def get_dpath(fpath) -> str:
    return Path(fpath).parent


  @classmethod
  def get_dname(cls, fpath) -> str:
    dpath = cls.get_dpath(fpath)
    return os.path.basename(dpath)


  @staticmethod
  def get_last_line(fpath) -> str|bytes:

    with open(fpath, "rb") as f:

      try:  # catch OSError in case of a one line file
        f.seek(-2, os.SEEK_END)
        while f.read(1) != b'\n':
          f.seek(-2, os.SEEK_CUR)
      except OSError:
        f.seek(0)
      last_line = f.readline().decode()

    return last_line
  

  @staticmethod
  def read_lines(
    f_in,
    f_out = None,
    mode = 'txt',
    last_line_number = None
  ) -> None:
    i = 1
    if f_out and mode == 'txt':
      for line in f_in:
        f_out.write(line)
        i += 1
        if last_line_number and i > last_line_number:
          break
    
    # For gzip for example
    elif f_out:
      for line in f_in:
        f_out.write(line.decode('utf-8'))
        i += 1
        if last_line_number and i > last_line_number:
          break
    
    else:
      for line in f_in:
        print(line)
        i += 1
        if last_line_number and i > last_line_number:
          break

  
  # python -m fire main.py modules File print_lines --fpath /home/sl/dev/partiseq-bp/pipeline/illumina_v3/sample/E150029294_L01_51_2.fq.gz --last-line-number 300 --out-fpath E150029294_L01_51_2-300-lines.fq
  @classmethod
  def print_lines(
    cls,
    fpath,
    out_fpath = None,
    last_line_number = None
  ) -> None:

    f_in = None
    f_out = None
    mode = 'txt'

    # Handle gz files
    if fpath.endswith(".gz"):
      f_in = gzip.open(fpath, "rb")
      mode = 'gz'
      
    # Treat any other file as a text file
    else:
      f_in = open(fpath, 'r')
    
    if out_fpath:
      f_out = open(out_fpath, 'w')
    
    cls.read_lines(f_in, f_out, mode, last_line_number)

    f_in.close()
    f_out.close()


  @staticmethod
  def append(in_fpath, out_fpath, remove_input_file=False):
    """
    Appends one file to another.
    :param in_fpath:
    :param out_fpath:
    :param remove_input_file:
    :return:
    """
    with open(in_fpath, 'r') as f_in, \
    open(out_fpath, 'a') as f_out:
      for line in f_in:
        # Prevent empty lines in concatenated fasta file
        if line == "\n":
          continue
        f_out.write(line)
    if remove_input_file:
      os.remove(in_fpath)
  
  @staticmethod
  def to_list(fpath) -> list:
    lines_list = []
    with open(fpath, 'r') as f:
      for line in f:
        lines_list.append(line.strip("\n"))
    return lines_list
  

  @staticmethod
  def concat_files(fpaths, out_fpath, has_header=False):
    already_handled_header = False
    with open(out_fpath, 'w') as f_out:
      for fpath in fpaths:
        with open(fpath) as f_in:
          if has_header:
            line = f_in.readline()
            if not already_handled_header:
              f_out.write(line)
              already_handled_header = True
          for line in f_in:
            f_out.write(line)

